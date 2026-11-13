from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import clinical_risk_evidence_hotfix1 as risk  # noqa: E402


class ClinicalRiskEvidenceTests(unittest.TestCase):
    def _write_fixture(self, root: Path) -> tuple[Path, Path]:
        patient = root / "patients" / "Demo 6"
        knowledge = patient / "Dante_results" / "knowledge_enrichment"
        knowledge.mkdir(parents=True)
        (knowledge / "Dante_knowledge_enrichment.state.json").write_text(
            json.dumps(
                {
                    "schema": risk.STAGE24_SCHEMA,
                    "status": "complete",
                    "rules_schema": risk.RULES_SCHEMA,
                    "input_profile": {"kind": "GENOTYPING_ARRAY_VCF", "assembly": "GRCh37"},
                }
            ),
            encoding="utf-8",
        )
        columns = [
            "finding_id", "event_key", "reason_code", "ref", "alt", "genotype", "zygosity",
            "gene", "resolved_rsid", "clinvar_significance", "clinvar_review_stars",
            "clinvar_conditions", "disease_name",
        ]
        rows = [
            {
                "finding_id": "F1",
                "event_key": "GRCh37:7:140453136:A:T",
                "reason_code": "K1_EXACT_CLINVAR_PLP",
                "ref": "A",
                "alt": "T",
                "genotype": "0/1",
                "zygosity": "HET_ALT",
                "gene": "BRAF",
                "resolved_rsid": "rs123",
                "clinvar_significance": "Pathogenic",
                "clinvar_review_stars": "3",
                "clinvar_conditions": "Noonan syndrome",
                "disease_name": "Noonan syndrome",
            }
        ]
        with (knowledge / "Dante_knowledge_findings_open.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

        open_root = root / "reference_data" / "open_reference"
        clingen = open_root / "ClinGen"
        clingen.mkdir(parents=True)
        adult = clingen / "actionability_adult.json"
        pediatric = clingen / "actionability_pediatric.json"
        adult.write_text(
            json.dumps(
                {
                    "rows": [
                        {
                            "geneOrVariant": "BRAF",
                            "condition": "Noonan syndrome",
                            "likelihood": "Penetrance approximately 80% by age 40",
                            "docId": "AC-NOONAN",
                        },
                        {
                            "geneOrVariant": "BRAF",
                            "condition": "Melanoma",
                            "likelihood": "Lifetime risk 90%",
                            "docId": "AC-MELANOMA",
                        },
                        {
                            "geneOrVariant": "BRAF",
                            "condition": "Noonan syndrome",
                            "likelihood": "3",
                            "docId": "AC-SCORE-ONLY",
                        },
                    ]
                }
            ),
            encoding="utf-8",
        )
        pediatric.write_text(json.dumps({"rows": []}), encoding="utf-8")
        manifest = {
            "schema": "open-reference-bundle-v2",
            "status": "complete",
            "sources": {
                "ClinGen": {
                    "artifacts": {
                        "actionability_adult": {"path": "ClinGen/actionability_adult.json"},
                        "actionability_pediatric": {"path": "ClinGen/actionability_pediatric.json"},
                    }
                }
            },
        }
        (open_root / "open_reference_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return patient, knowledge

    def test_local_builder_materializes_patient_allele_classification_and_matching_risk(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            patient, _ = self._write_fixture(root)
            self.assertEqual(risk.build(root, patient), 0)
            out = risk.output_paths(patient)
            payload = json.loads(out["json"].read_text(encoding="utf-8"))
            event = payload["events"]["GRCh37:7:140453136:A:T"]
            self.assertEqual(event["patient_alt_allele"], "T")
            self.assertEqual(event["patient_genotype_display"], "A/T")
            self.assertEqual(event["clinvar_classification"], "Pathogenic")
            self.assertEqual(len(event["metrics"]), 1)
            self.assertEqual(event["metrics"][0]["value"], "80")
            self.assertEqual(event["metrics"][0]["condition"], "Noonan syndrome")
            self.assertEqual(event["metrics"][0]["scope"], "gene-condition")
            self.assertNotIn("90", event["risk_summary_ru"])

            state = json.loads(out["state"].read_text(encoding="utf-8"))
            self.assertFalse(state["network_used"])
            self.assertFalse(state["patient_identifiers_sent_external"])
            self.assertEqual(state["counts"]["events_with_quantitative_metric"], 1)
            self.assertEqual(risk.is_ready(root, patient), (True, "READY"))

    def test_stage24_change_invalidates_risk_evidence_without_touching_stage24(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            patient, knowledge = self._write_fixture(root)
            risk.build(root, patient)
            self.assertTrue(risk.is_ready(root, patient)[0])
            open_path = knowledge / "Dante_knowledge_findings_open.tsv"
            open_path.write_text(open_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            ready, reason = risk.is_ready(root, patient)
            self.assertFalse(ready)
            self.assertEqual(reason, "STALE INPUTS")
            stage_state = json.loads((knowledge / "Dante_knowledge_enrichment.state.json").read_text(encoding="utf-8"))
            self.assertEqual(stage_state["status"], "complete")

    def test_bare_clingen_likelihood_score_is_not_quantitative_risk(self) -> None:
        self.assertIsNone(risk.parse_metric("$.likelihood", "3"))
        self.assertIsNone(risk.parse_metric("$.overallScore", "3"))


if __name__ == "__main__":
    unittest.main()
