from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import vcf_23me_converter as conv  # noqa: E402


VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##reference={assembly}\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
)


class Vcf23meConverterTests(unittest.TestCase):
    def run_convert(self, body: str, assembly: str = "GRCh37", gz: bool = False):
        td = tempfile.TemporaryDirectory()
        root = Path(td.name)
        source = root / ("input.vcf.gz" if gz else "input.vcf")
        payload = VCF_HEADER.format(assembly=assembly) + body
        if gz:
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write(payload)
        else:
            source.write_text(payload, encoding="utf-8")
        output = root / "Dante_23andMe_v5.txt"
        audit = root / "Dante_23andMe_v5_unresolved.tsv"
        report = root / "Dante_23andMe_v5_report.txt"
        conv.convert(source, output, audit, report, resolve_rsids=False)
        return td, output, audit, report

    def test_called_snv_forms_are_preserved_without_rsid_loss(self) -> None:
        td, output, audit, report = self.run_convert(
            "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\n"
            "1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
            "1\t300\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1\n"
            "X\t400\trs4\tC\tT\t.\tPASS\t.\tGT\t1\n"
            "2\t500\trs5\tA\tC,G\t.\tPASS\t.\tGT\t1/2\n"
            "3\t600\t.\tT\tC\t.\tPASS\t.\tGT\t0|1\n"
        )
        try:
            text = output.read_text(encoding="utf-8")
            self.assertIn("rs1\t1\t100\tAA", text)
            self.assertIn("rs2\t1\t200\tCT", text)
            self.assertIn("rs3\t1\t300\tAA", text)
            self.assertIn("rs4\tX\t400\tT", text)
            self.assertIn("rs5\t2\t500\tCG", text)
            self.assertIn("DANTE_GRCh37_3_600_T_C\t3\t600\tTC", text)
            audit_text = audit.read_text(encoding="utf-8")
            self.assertIn("NO_RSID_EXPORTED_WITH_DANTE_MARKER", audit_text)
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("called_canonical_snv_records: 6", report_text)
            self.assertIn("exported_called_snv_records: 6", report_text)
            self.assertIn("meaningful_called_snv_loss_records: 0", report_text)
            self.assertIn("loss_accounting_invariant: PASS", report_text)
            self.assertIn("phased_gt_collapsed_to_unphased_columns: 1", report_text)
        finally:
            td.cleanup()

    def test_non_snv_missing_gt_and_polyploid_are_audited_not_silently_lost(self) -> None:
        td, output, audit, report = self.run_convert(
            "1\t100\trs1\tA\tAT\t.\tPASS\t.\tGT\t0/1\n"
            "1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t./.\n"
            "1\t300\trs3\tG\tA\t.\tPASS\t.\tGT\t0/1/1\n"
            "1\t400\trs4\tT\tC\t.\tq10\t.\tGT\t0/1\n"
        )
        try:
            text = output.read_text(encoding="utf-8")
            self.assertIn("rs4\t1\t400\tTC", text)
            audit_text = audit.read_text(encoding="utf-8")
            self.assertIn("NON_SNV", audit_text)
            self.assertIn("MISSING_OR_INVALID_GT", audit_text)
            self.assertIn("GT_NOT_REPRESENTABLE_IN_23ANDME_COLUMNS", audit_text)
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("source_filter_nonpass_exported_records: 1", report_text)
            self.assertIn("meaningful_called_snv_loss_records: 0", report_text)
        finally:
            td.cleanup()

    def test_grch38_is_explicitly_labelled_and_not_fake_lifted(self) -> None:
        td, output, _, report = self.run_convert(
            "chr1\t101\trs10\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            assembly="GRCh38",
            gz=True,
        )
        try:
            text = output.read_text(encoding="utf-8")
            self.assertIn("# reference_build=GRCh38", text)
            self.assertIn("no silent GRCh37 liftover", text)
            self.assertIn("rs10\t1\t101\tAG", text)
            self.assertIn("Assembly: GRCh38", report.read_text(encoding="utf-8"))
        finally:
            td.cleanup()

    def test_unknown_assembly_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.vcf"
            source.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"
                "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(RuntimeError, "assembly is unknown"):
                conv.convert(
                    source,
                    root / "out.txt",
                    root / "audit.tsv",
                    root / "report.txt",
                    resolve_rsids=False,
                )

    def test_runtime_overlay_preserves_fallback_chain(self) -> None:
        runtime = (SRC / "vcf_export_runtime.ps1").read_text(encoding="utf-8-sig")
        self.assertIn("$script:DantePreVcfRunRsidExport", runtime)
        self.assertIn("$script:DantePreVcfSelectMode", runtime)
        self.assertIn("if (-not (Test-IsGenericSingleVcfProfile))", runtime)
        self.assertIn("Dante_23andMe_v5.txt", runtime)
        self.assertIn("DANTE_VCF_23ME_NO_RSID_RESOLUTION", runtime)
        self.assertIn("Docker is not required", runtime)

    def test_release_manifest_deploys_converter_runtime_and_tests(self) -> None:
        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["pipeline_release"], "9.1.0")
        self.assertEqual(manifest["display_version"], "9.1.0")
        for path in (
            "src/vcf_23me_converter.py",
            "src/vcf_export_runtime.ps1",
            "tests/test_vcf_23me_converter.py",
            "docs/VCF_23ME_V5_CONVERTER_RU.md",
        ):
            self.assertIn(path, manifest["deployment_paths"])


if __name__ == "__main__":
    unittest.main()
