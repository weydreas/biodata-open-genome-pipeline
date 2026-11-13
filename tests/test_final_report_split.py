from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
import genome_report as gr  # noqa: E402
import genome_report_hotfix9 as gr_impl  # noqa: E402


class FinalReportSplitTests(unittest.TestCase):
    def test_commercial_usage_rejects_research_only_sources(self) -> None:
        usage = {
            "schema": "knowledge-usage-v1",
            "scope": "OPEN_PLUS_RESEARCH_OVERLAYS",
            "sources_used": [
                {"source": "PanelApp", "publication_scope": "RESEARCH_ONLY", "prohibited_uses": ["commercial"]}
            ],
        }
        with self.assertRaises(RuntimeError):
            gr.validate_usage("commercial", usage)

    def test_reason_grouping_deduplicates_events_but_keeps_all_relation_rows(self) -> None:
        rows = [
            {"finding_id": "OF1", "priority_rank": "2", "event_key": "GRCh37:1:1:A:G", "gene": "GENE1", "disease_name": "D1", "reason_code": "K3_RARE_MISSENSE_VALIDATED_GENE"},
            {"finding_id": "OF2", "priority_rank": "2", "event_key": "GRCh37:1:1:A:G", "gene": "GENE1", "disease_name": "D2", "reason_code": "K3_RARE_MISSENSE_VALIDATED_GENE"},
            {"finding_id": "OF3", "priority_rank": "4", "event_key": "GRCh37:2:2:C:T", "gene": "GENE2", "reason_code": "K3_RARE_MISSENSE_VALIDATED_GENE"},
        ]
        grouped = gr.group_reason_events(rows)
        self.assertEqual(len(grouped["K3_RARE_MISSENSE_VALIDATED_GENE"]), 2)
        self.assertEqual(len(grouped["K3_RARE_MISSENSE_VALIDATED_GENE"][0]["rows"]), 2)
        counts = gr.reason_counts(grouped, gr.OPEN_CATEGORY_ORDER)
        missense = next(x for x in counts if x["reason"] == "K3_RARE_MISSENSE_VALIDATED_GENE")
        self.assertEqual(missense["events"], 2)
        self.assertEqual(missense["rows"], 3)

    def test_common_major_filter_and_heading(self) -> None:
        rows = [
            {
                "finding_id": "OF1", "priority_rank": "4", "event_type": "SNV_INDEL",
                "event_key": "GRCh37:11:17435840:A:G", "chrom": "11", "position": "17435840",
                "ref": "A", "alt": "G", "genotype": "1/1", "zygosity": "HOM_ALT",
                "gene": "ABCC8", "max_af": "0.9834", "resolved_rsid": "rs4148632",
                "reason_code": "K1_EXACT_CLINVAR_CONFLICT",
            },
            {
                "finding_id": "OF2", "priority_rank": "4", "event_type": "SNV_INDEL",
                "event_key": "GRCh37:11:17435840:A:G", "chrom": "11", "position": "17435840",
                "ref": "A", "alt": "G", "genotype": "1/1", "zygosity": "HOM_ALT",
                "gene": "LOC110121471", "max_af": "0.9834", "resolved_rsid": "rs4148632",
                "reason_code": "K1_EXACT_CLINVAR_CONFLICT",
            },
            {
                "finding_id": "OF3", "priority_rank": "2", "event_type": "SNV_INDEL",
                "event_key": "GRCh37:1:10:C:T", "chrom": "1", "position": "10",
                "ref": "C", "alt": "T", "genotype": "0/1", "gene": "GENE1",
                "max_af": "0.001", "resolved_rsid": "rs10",
                "reason_code": "K3_RARE_MISSENSE_VALIDATED_GENE",
            },
        ]
        grouped = gr.group_reason_events(rows)
        primary, reference = gr.split_common_major_events(grouped)
        self.assertEqual(sum(len(v) for v in primary.values()), 1)
        self.assertEqual(sum(len(v) for v in reference.values()), 1)
        event = reference["K1_EXACT_CLINVAR_CONFLICT"][0]
        self.assertAlmostEqual(gr.event_max_af(event), 0.9834)
        heading = gr.event_heading_label(event)
        self.assertIn("rs4148632", heading)
        self.assertIn("(ABCC8, LOC110121471)", heading)
        self.assertIn("аллель пациента G/G", heading)

    def test_common_major_threshold_is_strict_and_unknown_af_stays_primary(self) -> None:
        rows = [
            {"finding_id": "A", "event_key": "E1", "event_type": "SNV_INDEL", "ref": "A", "alt": "G", "genotype": "0/1", "max_af": "0.70", "reason_code": "K1_EXACT_CLINVAR_ASSOCIATION"},
            {"finding_id": "B", "event_key": "E2", "event_type": "SNV_INDEL", "ref": "A", "alt": "G", "genotype": "0/1", "max_af": "0.70001", "reason_code": "K1_EXACT_CLINVAR_ASSOCIATION"},
            {"finding_id": "C", "event_key": "E3", "event_type": "SNV_INDEL", "ref": "A", "alt": "G", "genotype": "0/1", "max_af": ".", "reason_code": "K1_EXACT_CLINVAR_ASSOCIATION"},
        ]
        primary, reference = gr.split_common_major_events(gr.group_reason_events(rows))
        primary_keys = {event["event_key"] for events in primary.values() for event in events}
        reference_keys = {event["event_key"] for events in reference.values() for event in events}
        self.assertEqual(primary_keys, {"E1", "E3"})
        self.assertEqual(reference_keys, {"E2"})

    def test_research_filter_uses_open_event_af_context(self) -> None:
        overlay_rows = [
            {"overlay_id": "RO1", "event_key": "E1", "gene": "GENE1", "reason_code": "K8_CADD_SUPPORT_TARGET"},
            {"overlay_id": "RO2", "event_key": "E2", "gene": "GENE2", "reason_code": "K8_CADD_SUPPORT_TARGET"},
        ]
        open_by_event = {
            "E1": [{"event_key": "E1", "max_af": "0.95"}],
            "E2": [{"event_key": "E2", "max_af": "0.001"}],
        }
        primary, reference = gr.split_common_major_events(gr.group_reason_events(overlay_rows), open_by_event)
        self.assertEqual({e["event_key"] for es in primary.values() for e in es}, {"E2"})
        self.assertEqual({e["event_key"] for es in reference.values() for e in es}, {"E1"})

    def test_commercial_report_does_not_require_or_open_research_layer(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            patient = root / "patients" / "Demo 5"
            knowledge = patient / "Dante_results" / "knowledge_enrichment"
            knowledge.mkdir(parents=True)
            state = {
                "status": "complete",
                "schema": "dante-knowledge-enrichment-v1",
                "rules_schema": "dante-knowledge-rules-v8",
                "pipeline_release": "9.0.0",
                "assembly": "GRCh37",
                "scan_stats": {"candidate_events": 2, "open_findings": 2, "research_overlays": 999},
            }
            (knowledge / "Dante_knowledge_enrichment.state.json").write_text(json.dumps(state), encoding="utf-8")
            fields = [
                "finding_id", "priority_rank", "event_type", "event_key", "chrom", "position", "ref", "alt",
                "genotype", "zygosity", "gene", "consequence", "max_af", "resolved_rsid", "clinvar_significance",
                "clinvar_conditions", "disease_id", "disease_name", "gene_disease_classification", "inheritance",
                "technical_confidence", "evidence_sources", "reason_code",
            ]
            with (knowledge / "Dante_knowledge_findings_open.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow({
                    "finding_id": "OF_RARE", "priority_rank": "2", "event_type": "SNV_INDEL",
                    "event_key": "GRCh37:1:1:A:G", "chrom": "1", "position": "1", "ref": "A", "alt": "G",
                    "genotype": "0/1", "zygosity": "HET_ALT", "gene": "GENE1", "consequence": "missense_variant",
                    "max_af": "0.001", "resolved_rsid": "rs1", "clinvar_significance": "Pathogenic",
                    "disease_name": "D1", "gene_disease_classification": "STRONG", "inheritance": "AD",
                    "technical_confidence": "CATALOG_PASS", "evidence_sources": "ClinVar;VEP/CATALOG",
                    "reason_code": "K1_EXACT_CLINVAR_PLP",
                })
                writer.writerow({
                    "finding_id": "OF_COMMON", "priority_rank": "4", "event_type": "SNV_INDEL",
                    "event_key": "GRCh37:11:17435840:A:G", "chrom": "11", "position": "17435840", "ref": "A", "alt": "G",
                    "genotype": "1/1", "zygosity": "HOM_ALT", "gene": "ABCC8", "consequence": "intron_variant",
                    "max_af": "0.9834", "resolved_rsid": "rs4148632", "clinvar_significance": "Conflicting_classifications_of_pathogenicity",
                    "technical_confidence": "STANDARD", "evidence_sources": "ClinVar;VEP/CATALOG",
                    "reason_code": "K1_EXACT_CLINVAR_CONFLICT",
                })
            usage = {
                "schema": "knowledge-usage-v1",
                "scope": "OPEN_FINDINGS_ONLY",
                "sources_used": [{"source": "ClinVar", "publication_scope": "COMMERCIAL_ALLOWED", "prohibited_uses": []}],
            }
            (knowledge / "Dante_knowledge_usage_manifest.json").write_text(json.dumps(usage), encoding="utf-8")
            # Deliberately do NOT create research TSV or research usage manifest.
            def fake_pdf(output: Path, **kwargs) -> None:
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_bytes(b"%PDF-1.4\n% synthetic unit-test pdf\n")

            argv = ["genome_report.py", "commercial", "--root", str(root), "--patient-root", str(patient)]
            with mock.patch.object(sys, "argv", argv), mock.patch.object(gr_impl, "render_report_pdf", side_effect=fake_pdf):
                self.assertEqual(gr.main(), 0)
            report_dir = patient / "Dante_results" / "reports"
            manifest = json.loads((report_dir / "Dante_genome_report_COMMERCIAL.manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema"], "dante-genome-report-v2")
            self.assertEqual(manifest["generation_mode"], "DETERMINISTIC_STAGE24_ONLY")
            self.assertFalse(manifest["ai_used"])
            self.assertTrue(manifest["all_stage24_findings_included"])
            self.assertTrue(manifest["all_stage24_findings_accounted_for"])
            self.assertTrue(manifest["primary_pdf_filtered"])
            self.assertEqual(manifest["finding_rows"], 1)
            self.assertEqual(manifest["common_major_filter"]["reference_rows"], 1)
            self.assertEqual(manifest["common_major_filter"]["rule"], "max_af > 0.70")
            self.assertFalse(manifest["research_layer_opened"])
            self.assertNotIn("research_overlays", manifest["source_files"])
            self.assertTrue((report_dir / "Dante_genome_report_COMMERCIAL.pdf").is_file())
            self.assertTrue((report_dir / "Dante_genome_report_COMMERCIAL_COMMON_MAJOR_REFERENCE.pdf").is_file())

    def test_research_mode_has_explicit_restriction_acknowledgement_contract(self) -> None:
        source = (SRC / "genome_report_hotfix9.py").read_text(encoding="utf-8")
        self.assertIn("Research mode requires --accept-research-restrictions", source)
        self.assertIn("Type RESEARCH to continue", source)
        self.assertIn("RESEARCH ONLY", source)
        self.assertIn("non-commercial", source)

    def test_report_generator_has_no_openai_or_top_n_contract(self) -> None:
        source = (SRC / "genome_report_hotfix9.py").read_text(encoding="utf-8")
        entrypoint = (SRC / "genome_report.py").read_text(encoding="utf-8")
        for text in (source, entrypoint):
            self.assertNotIn("from openai import OpenAI", text)
            self.assertNotIn("client.responses.create", text)
            self.assertNotIn("OPENAI_MAX_VARIANTS", text)
            self.assertNotIn("--max-events", text)
        self.assertIn('"ai_used": False', source)
        self.assertIn('"all_stage24_findings_accounted_for": True', source)

    def test_materialized_menu_exposes_two_separate_no_ai_report_actions(self) -> None:
        runtime = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8-sig")
        self.assertIn('Write-MenuItem "25" "Generate CLIENT PDF report (OPEN / COMMERCIAL-SAFE / NO AI)"', runtime)
        self.assertIn('Write-MenuItem "26" "Generate RESEARCH ONLY PDF report (NON-COMMERCIAL / NO AI)"', runtime)
        self.assertIn('Write-MenuItem "27" "Launch local IGV.js browser', runtime)
        self.assertIn('"25" { return "commercial-report" }', runtime)
        self.assertIn('"26" { return "research-report" }', runtime)
        self.assertIn('"27" { return "genome-browser" }', runtime)
        self.assertIn('"commercial-report" { Run-GenomeReport "commercial" }', runtime)
        self.assertIn('"research-report" { Run-GenomeReport "research" }', runtime)
        self.assertIn('dante-genome-report-v2', runtime)
        report_block = runtime.split('function Run-GenomeReport', 1)[1].split('function Select-Mode', 1)[0]
        self.assertNotIn('Ensure-OpenAiSdk', report_block)
        self.assertNotIn('OPENAI_MAX_VARIANTS', report_block)
        self.assertNotIn('--max-events', report_block)

    def test_start_launcher_does_not_run_openai_preflight(self) -> None:
        launcher = (ROOT / "start.cmd").read_text(encoding="utf-8")
        self.assertNotIn("DANTE_OPENAI_PREFLIGHT", launcher)
        self.assertNotIn("openai_sdk_preflight.py", launcher)
        self.assertIn("deterministic Stage-24 PDF renderers", launcher)

    def test_deployment_manifest_keeps_hotfix9_report_files_in_release_910(self) -> None:
        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        self.assertNotIn("src/openai_sdk_preflight.py", manifest["deployment_paths"])
        self.assertIn("src/genome_report.py", manifest["deployment_paths"])
        self.assertIn("src/genome_report_base_v2.py", manifest["deployment_paths"])
        self.assertIn("src/genome_report_hotfix9.py", manifest["deployment_paths"])
        self.assertEqual(manifest["pipeline_release"], "9.1.0")
        self.assertEqual(manifest["hotfix"], 0)
        self.assertEqual(manifest["display_version"], "9.1.0")

    def test_report_generator_self_test(self) -> None:
        self.assertEqual(gr.self_test(), 0)


if __name__ == "__main__":
    unittest.main()
