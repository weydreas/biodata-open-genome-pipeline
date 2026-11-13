from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import clinical_risk_evidence_hotfix1 as risk  # noqa: E402
import genome_report_hotfix15 as gr15  # noqa: E402


class ClinicalRiskReportTests(unittest.TestCase):
    def setUp(self) -> None:
        gr15.risk = risk
        gr15._register_events(
            {
                "GRCh37:7:140453136:A:T": {
                    "event_key": "GRCh37:7:140453136:A:T",
                    "rsid": "rs123",
                    "genes": ["BRAF"],
                    "patient_alt_allele": "T",
                    "patient_genotype_display": "A/T",
                    "clinvar_classification": "Pathogenic",
                    "metrics": [
                        {
                            "risk_summary_ru": "Пенетрантность 80% к возрасту 40 лет (ClinGen Actionability, Adult; gene-condition)"
                        }
                    ],
                }
            }
        )

    def test_stable_entrypoint_routes_through_hotfix15(self) -> None:
        source = (SRC / "genome_report.py").read_text(encoding="utf-8")
        self.assertIn("genome_report_hotfix15", source)
        self.assertIn("genome_report_hotfix14", source)
        self.assertIn("genome_report_hotfix12", source)
        self.assertIn("clinical_risk_evidence_hotfix1", source)
        self.assertIn("_impl.PUBLICATION_CONTRACT = _impl.h10.PUBLICATION_CONTRACT", source)
        self.assertIn("_impl.h11.PUBLICATION_CONTRACT = _impl.h10.PUBLICATION_CONTRACT", source)

    def test_risk_layer_has_no_network_client_and_preserves_explicit_unavailable(self) -> None:
        source = (SRC / "clinical_risk_evidence.py").read_text(encoding="utf-8")
        self.assertNotIn("urllib.request", source)
        self.assertNotIn("requests.", source)
        self.assertIn("patient_identifiers_sent_external", source)
        self.assertIn("Количественная оценка риска не установлена", source)
        self.assertIn("Bare ClinGen", source)

    @unittest.skipUnless(importlib.util.find_spec("reportlab") is not None, "reportlab not installed in CI runtime")
    def test_short_variant_index_is_augmented_to_six_columns(self) -> None:
        from reportlab import platypus
        from reportlab.lib.styles import getSampleStyleSheet

        style = getSampleStyleSheet()["BodyText"]
        p = lambda text: platypus.Paragraph(text, style)
        source_rows = [
            [p("rsID / вариант"), p("Ген(ы)"), p("Связанные состояния")],
            [p("rs123"), p("BRAF"), p("Noonan syndrome")],
        ]
        with gr15.enriched_summary_tables():
            table = platypus.Table(source_rows, colWidths=[100, 100, 100])
        self.assertEqual(len(table._cellvalues[0]), 6)
        headers = [cell.getPlainText() for cell in table._cellvalues[0]]
        self.assertEqual(
            headers,
            [
                "rsID / вариант",
                "Аллель пациента",
                "Ген(ы)",
                "Классификация ClinVar",
                "Риск / пенетрантность",
                "Связанные состояния",
            ],
        )
        row = [cell.getPlainText() for cell in table._cellvalues[1]]
        self.assertEqual(row[1], "T; A/T")
        self.assertEqual(row[3], "Pathogenic")
        self.assertIn("80%", row[4])
        self.assertEqual(row[5], "Noonan syndrome")

    def test_manifest_wrapper_records_risk_artifacts_and_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            patient_results = Path(directory) / "Dante_results"
            report_dir = patient_results / "reports"
            risk_dir = patient_results / "clinical_risk_evidence"
            report_dir.mkdir(parents=True)
            risk_dir.mkdir(parents=True)
            evidence = risk_dir / "Dante_clinical_risk_evidence.json"
            evidence_tsv = risk_dir / "Dante_clinical_risk_evidence.tsv"
            state_path = risk_dir / "Dante_clinical_risk_evidence.state.json"
            evidence.write_text(json.dumps({"schema": risk.SCHEMA, "status": "complete"}), encoding="utf-8")
            evidence_tsv.write_text("event_key\n", encoding="utf-8")
            state_path.write_text(
                json.dumps(
                    {
                        "source_policy": risk.SOURCE_POLICY,
                        "network_used": False,
                        "patient_identifiers_sent_external": False,
                        "counts": {"events": 1},
                    }
                ),
                encoding="utf-8",
            )
            manifest_path = report_dir / "report.manifest.json"

            def original(**kwargs):
                kwargs["output_path"].write_text(
                    json.dumps({"source_files": kwargs.get("source_files", {}), "generation_mode": "old"}),
                    encoding="utf-8",
                )

            wrapped = gr15._manifest_with_risk(original)
            wrapped(output_path=manifest_path, source_files={"state": {"sha256": "x"}})
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["publication_contract"], gr15.PUBLICATION_CONTRACT)
            self.assertEqual(manifest["generation_mode"], gr15.GENERATION_MODE)
            self.assertEqual(manifest["clinical_risk_evidence_schema"], risk.SCHEMA)
            self.assertFalse(manifest["clinical_risk_network_used"])
            self.assertFalse(manifest["clinical_risk_patient_identifiers_sent_external"])
            self.assertIn("clinical_risk_evidence", manifest["source_files"])
            self.assertIn("clinical_risk_state", manifest["source_files"])

    def test_deployment_manifest_tracks_risk_layer(self) -> None:
        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        for path in (
            "src/clinical_risk_evidence.py",
            "src/clinical_risk_evidence_hotfix1.py",
            "src/clinical_risk_runtime.ps1",
            "src/genome_report_hotfix15.py",
            "tests/test_clinical_risk_evidence.py",
            "tests/test_report_clinical_risk.py",
        ):
            self.assertIn(path, manifest["deployment_paths"])


if __name__ == "__main__":
    unittest.main()
