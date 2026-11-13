from __future__ import annotations

import gzip
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import array_vcf_tools as avt  # noqa: E402
import genome_report_hotfix13 as gr13  # noqa: E402
import knowledge_enrichment as ke  # noqa: E402
import patient_folder_discovery as discovery  # noqa: E402
import patient_input_profile as pip  # noqa: E402
import patient_manager as pm  # noqa: E402


def write_array_vcf(path: Path) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##source=Genotek\n"
        "##reference=hg19\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0|0\n"
        "chr1\t200\trs2\tC\tT\t.\t.\t.\tGT\t0|1\n"
        "chr1\t300\trs3\tG\tA\t.\t.\t.\tGT\t./.\n"
        "chr2\t400\trs4\tA\tAT\t.\t.\t.\tGT\t1|1\n",
        encoding="utf-8",
    )


class PatientInputProfileTests(unittest.TestCase):
    def test_genotek_dense_vcf_is_array_not_wgs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            patient = Path(directory)
            inp = patient / "input"
            inp.mkdir()
            write_array_vcf(inp / "sample.vcf")
            profile = pip.detect(patient)
        self.assertEqual(profile["kind"], pip.ARRAY_KIND)
        self.assertEqual(profile["source"], "Genotek")
        self.assertEqual(profile["assembly"], "GRCh37")
        self.assertFalse(profile["capabilities"]["whole_genome_negative_evidence"])
        self.assertFalse(profile["capabilities"]["structural_variant_annotation"])
        self.assertTrue(profile["capabilities"]["direct_rsid_export"])

    def test_startup_discovery_persists_array_profile_for_existing_patient(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            patient = pm.ensure_patient(root, "Demo 2")
            write_array_vcf(patient / "input" / "y6tvv2.vcf")

            refreshed = discovery.refresh_input_profiles(root)
            saved = json.loads((patient / "input_profile.json").read_text(encoding="utf-8"))

        self.assertIn(("Demo 2", pip.ARRAY_KIND), refreshed)
        self.assertEqual(saved["kind"], pip.ARRAY_KIND)
        self.assertEqual(saved["source"], "Genotek")
        self.assertEqual(saved["assembly"], "GRCh37")
        self.assertFalse(saved["capabilities"]["whole_genome_negative_evidence"])

    def test_measured_homref_is_kept_only_in_direct_export(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "array.vcf"
            write_array_vcf(source)
            output = root / "Dante_23andMe_style.txt.gz"
            avt.export_23andme(source, output, root / "export.txt")
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                rows = [line.strip() for line in handle if not line.startswith("#")]
        self.assertIn("rs1\t1\t100\tAA", rows)
        self.assertIn("rs2\t1\t200\tCT", rows)
        self.assertFalse(any(row.startswith("rs3\t") for row in rows))
        self.assertFalse(any(row.startswith("rs4\t") for row in rows))

    def test_annotation_layers_include_alt_only_and_strip_unproven_phase(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "array.vcf"
            write_array_vcf(source)
            snp = root / "snp.vcf.gz"
            indel = root / "indel.vcf.gz"
            residual = root / "residual.vcf.gz"
            avt.prepare_layers(source, snp, indel, residual, root / "prep.txt")
            with gzip.open(snp, "rt", encoding="utf-8") as handle:
                snp_text = handle.read()
            with gzip.open(indel, "rt", encoding="utf-8") as handle:
                indel_text = handle.read()
        self.assertNotIn("rs1", snp_text)
        self.assertNotIn("rs3", snp_text)
        self.assertIn("rs2", snp_text)
        self.assertIn("\tGT\t0/1\n", snp_text)
        self.assertIn("DANTE_ARRAY_PHASE_STRIPPED", snp_text)
        self.assertIn("rs4", indel_text)
        self.assertIn("\tGT\t1/1\n", indel_text)
        self.assertIn("\tPASS\t", snp_text)

    def test_array_runtime_hides_wgs_only_menu_actions(self) -> None:
        source = (SRC / "array_profile_runtime.ps1").read_text(encoding="utf-8")
        menu = source.split("function Select-Mode", 1)[1]
        for item in ("3", "4", "5", "13", "14", "15", "16", "17", "18", "19", "27"):
            self.assertNotIn(f'Write-MenuItem "{item}"', menu)
        for item in ("1", "2", "6", "7", "8", "9", "10", "11", "12", "21", "22", "23", "24", "25", "26"):
            self.assertIn(f'Write-MenuItem "{item}"', menu)
        self.assertIn("interrogated array loci only", menu)
        self.assertIn('"2" { return "rsid" }', menu)

    def test_runtime_overlay_is_array_only_and_preserves_wgs_fallback(self) -> None:
        array_runtime = (SRC / "array_profile_runtime.ps1").read_text(encoding="utf-8-sig")
        self.assertIn('$script:DanteOriginalSelectMode', array_runtime)
        self.assertIn('if (-not (Test-IsArrayProfile))', array_runtime)
        self.assertIn('& $script:DanteOriginalRunRsidExport', array_runtime)
        self.assertIn('& $script:DanteOriginalPrepareNormalizedLayers', array_runtime)
        self.assertIn('& $script:DanteOriginalRunFullVepAnnotation', array_runtime)

    def test_array_stage24_confidence_is_not_read_backed_wgs(self) -> None:
        event = ke.PatientEvent("SNV_INDEL", "v1", {"CATALOG"}, "1", 100, 100, ref="A", alt="G")
        old = os.environ.get("DANTE_INPUT_PROFILE_KIND")
        try:
            os.environ["DANTE_INPUT_PROFILE_KIND"] = "GENOTYPING_ARRAY_VCF"
            self.assertTrue(ke.event_technical_confidence(event).startswith("ARRAY_VENDOR_CALL_NO_READ_LEVEL_QC"))
        finally:
            if old is None:
                os.environ.pop("DANTE_INPUT_PROFILE_KIND", None)
            else:
                os.environ["DANTE_INPUT_PROFILE_KIND"] = old

    def test_stage24_state_records_array_scope(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            patient = Path(directory)
            inp = patient / "input"
            inp.mkdir()
            write_array_vcf(inp / "sample.vcf")
            profile = pip.detect(patient)
            pip.save_atomic(patient / "input_profile.json", profile)
            out = patient / "Dante_results" / "knowledge_enrichment"
            out.mkdir(parents=True)
            report = out / "Dante_knowledge_enrichment_report.txt"
            report.write_text("ORIGINAL REPORT\n", encoding="utf-8")
            state = {
                "schema": "dante-knowledge-enrichment-v1",
                "status": "complete",
                "rules_schema": "dante-knowledge-rules-v8",
                "outputs": {
                    "report": {
                        "name": report.name,
                        "length": report.stat().st_size,
                        "sha256": pip.sha256_file(report),
                    }
                },
            }
            pip.save_atomic(out / "Dante_knowledge_enrichment.state.json", state)
            updated = pip.annotate_stage24(patient)
            text = report.read_text(encoding="utf-8")
        self.assertEqual(updated["input_profile"]["kind"], "GENOTYPING_ARRAY_VCF")
        self.assertFalse(updated["input_profile"]["whole_genome_negative_evidence"])
        self.assertIn("Structural variants: NOT TESTED", text)
        self.assertIn("Whole-genome negative evidence: NO", text)

    def test_array_pdf_scope_never_claims_wgs_coverage(self) -> None:
        self.assertIn("не отчет по полногеномному секвенированию (WGS)", gr13.ARRAY_COMMERCIAL_LEGAL[0])
        self.assertTrue(any("UNKNOWN" in line for line in gr13.ARRAY_COMMERCIAL_LEGAL))
        state = {"input_profile": {"kind": "GENOTYPING_ARRAY_VCF"}}
        self.assertTrue(gr13.is_array_state(state))
        self.assertNotIn(3, gr13.ARRAY_STAGE_ITEMS)
        self.assertNotIn(15, gr13.ARRAY_STAGE_ITEMS)
        self.assertIn(24, gr13.ARRAY_STAGE_ITEMS)

    def test_deployment_manifest_tracks_profile_runtime(self) -> None:
        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        for path in (
            "src/patient_input_profile.py",
            "src/array_vcf_tools.py",
            "src/array_profile_runtime.ps1",
            "src/report_runtime.ps1",
            "src/knowledge_enrichment_hotfix7.py",
            "src/genome_report_hotfix13.py",
            "tests/test_patient_input_profiles.py",
        ):
            self.assertIn(path, manifest["deployment_paths"])


if __name__ == "__main__":
    unittest.main()
