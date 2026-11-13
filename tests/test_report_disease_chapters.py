from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import genome_report_hotfix11 as gr11  # noqa: E402
import genome_report_base_v2 as base  # noqa: E402
import report_disease_chapters as chapters  # noqa: E402


class DiseaseChapterTaxonomyTests(unittest.TestCase):
    def test_representative_clinvar_conditions_map_to_expected_chapters(self) -> None:
        cases = {
            "Breast_and_colorectal_cancer,_susceptibility_to": "D01_ONCOLOGY",
            "Familial_hypercholesterolemia": "D02_CARDIOVASCULAR",
            "Maturity-onset_diabetes_of_the_young": "D03_ENDOCRINE",
            "Maple_syrup_urine_disease": "D04_INBORN_METABOLIC",
            "Mitochondrial_complex_I_deficiency,_nuclear_type_1": "D04_INBORN_METABOLIC",
            "Amyotrophic_lateral_sclerosis": "D05_NEUROLOGY",
            "Intellectual_developmental_disorder,_X-linked_111": "D06_NEURODEVELOPMENTAL",
            "Hereditary_spherocytosis_type_3": "D07_HEMATOLOGY",
            "Rheumatoid_arthritis": "D08_IMMUNE",
            "Nephrolithiasis,_calcium_oxalate": "D09_RENAL",
            "Hereditary_pancreatitis": "D10_GASTRO",
            "Pulmonary_disease,_chronic_obstructive,_susceptibility_to": "D11_RESPIRATORY",
            "Retinitis_pigmentosa_58": "D12_OPHTHALMOLOGY",
            "Autosomal_recessive_nonsyndromic_hearing_loss_59": "D13_HEARING_ENT",
            "Diastrophic_dysplasia": "D14_MUSCULOSKELETAL",
            "Oculocutaneous_albinism_type_1A": "D15_DERMATOLOGY",
            "Bilateral_cryptorchidism": "D16_REPRODUCTIVE",
            "Noonan_syndrome_and_Noonan-related_syndrome": "D17_CONGENITAL_MULTISYSTEM",
            "COVID-19_susceptibility": "D18_INFECTIOUS",
            "Clopidogrel_response": "D19_PHARMACOGENETICS",
            "C3HEX,_ability_to_smell": "D20_TRAITS",
            "not_provided": "D99_UNCLASSIFIED",
        }
        for condition, expected in cases.items():
            with self.subTest(condition=condition):
                self.assertEqual(chapters.classify_condition(condition), expected)

    def test_mixed_condition_event_has_one_primary_and_secondary_chapter(self) -> None:
        hnf1a = chapters.assign_conditions([
            "Maturity-onset_diabetes_of_the_young|Type_2_diabetes_mellitus|Nonpapillary_renal_cell_carcinoma"
        ])
        self.assertEqual(hnf1a["primary"], "D03_ENDOCRINE")
        self.assertIn("D01_ONCOLOGY", hnf1a["secondary"])

        ugt1a1 = chapters.assign_conditions([
            "Gilbert_syndrome|Crigler-Najjar_syndrome|Irinotecan_response"
        ])
        self.assertEqual(ugt1a1["primary"], "D04_INBORN_METABOLIC")
        self.assertIn("D19_PHARMACOGENETICS", ugt1a1["secondary"])

    def test_generic_noise_does_not_drive_chapter_selection(self) -> None:
        info = chapters.assign_conditions([
            "not_provided|not_specified|ABCC8-related_disorder|Maturity-onset_diabetes_of_the_young"
        ])
        self.assertEqual(info["primary"], "D03_ENDOCRINE")
        self.assertEqual(info["conditions"], ["Maturity-onset diabetes of the young"])


class DiseaseChapterPublicationTests(unittest.TestCase):
    def test_same_event_across_reason_codes_is_one_card_and_keeps_all_rows(self) -> None:
        rows = [
            {
                "finding_id": "OF1", "event_key": "E1", "reason_code": "K1_EXACT_CLINVAR_CONFLICT",
                "resolved_rsid": "rs28385701", "gene": "PCSK9", "priority_rank": "4",
                "clinvar_conditions": "Familial_hypercholesterolemia",
            },
            {
                "finding_id": "OF2", "event_key": "E1", "reason_code": "K4_CLINGEN_ACTIONABLE",
                "resolved_rsid": "rs28385701", "gene": "PCSK9", "priority_rank": "4",
                "disease_name": "Familial_hypercholesterolemia",
            },
        ]
        events = gr11.merge_reason_events(base.group_reason_events(rows), base.OPEN_CATEGORY_ORDER)
        self.assertEqual(len(events), 1)
        self.assertEqual(len(events[0]["rows"]), 2)
        self.assertEqual(set(events[0]["reason_codes"]), {"K1_EXACT_CLINVAR_CONFLICT", "K4_CLINGEN_ACTIONABLE"})

    def test_chapter_grouping_uses_stage24_conditions_and_preserves_display_conditions(self) -> None:
        rows = [
            {
                "finding_id": "OF1", "event_key": "E1", "reason_code": "K1_EXACT_CLINVAR_CONFLICT",
                "resolved_rsid": "rs6413458", "gene": "APOB", "priority_rank": "4",
                "clinvar_conditions": "Cardiovascular_phenotype|Familial_hypercholesterolemia|not_provided",
            },
            {
                "finding_id": "OF2", "event_key": "E2", "reason_code": "K1_EXACT_CLINVAR_CONFLICT",
                "resolved_rsid": "rs17879961", "gene": "CHEK2", "priority_rank": "4",
                "clinvar_conditions": "Hereditary_cancer-predisposing_syndrome|Colorectal_cancer|Prostate_cancer",
            },
        ]
        chaptered, metadata = gr11.group_events_by_chapter(base.group_reason_events(rows), base.OPEN_CATEGORY_ORDER)
        self.assertEqual([e["resolved_rsid"] for e in chaptered["D02_CARDIOVASCULAR"]], ["rs6413458"])
        self.assertEqual([e["resolved_rsid"] for e in chaptered["D01_ONCOLOGY"]], ["rs17879961"])
        self.assertIn("Familial hypercholesterolemia", metadata["E1"]["display_conditions"])
        self.assertNotIn("not provided", [x.lower() for x in metadata["E1"]["display_conditions"]])

    def test_renderer_contract_contains_requested_chapter_index_table_and_card_fields(self) -> None:
        source = (SRC / "genome_report_hotfix11.py").read_text(encoding="utf-8")
        self.assertIn("Краткий перечень вариантов", source)
        self.assertIn('p("rsID / вариант"', source)
        self.assertIn('p("Ген(ы)"', source)
        self.assertIn('p("Связанные состояния"', source)
        self.assertIn('("Клинический раздел", primary)', source)
        self.assertIn('("Дополнительные разделы", secondary)', source)
        self.assertIn("group_events_by_chapter", source)
        self.assertIn("merge_reason_events", source)

    def test_manifest_records_disease_chapter_schema_without_changing_stage24_rules(self) -> None:
        source = (SRC / "genome_report_hotfix11.py").read_text(encoding="utf-8")
        entrypoint = (SRC / "genome_report.py").read_text(encoding="utf-8")
        self.assertIn('"disease_chapter_schema"', source)
        self.assertIn("dante-disease-chapters-v1", (SRC / "report_disease_chapters.py").read_text(encoding="utf-8"))
        self.assertIn("genome_report_hotfix11", entrypoint)
        self.assertIn("_impl.PUBLICATION_CONTRACT = _impl.h10.PUBLICATION_CONTRACT", entrypoint)
        self.assertNotIn("dante-knowledge-rules-v9", source)

    def test_deployment_manifest_contains_chapter_modules_and_contract_doc(self) -> None:
        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        for path in [
            "src/report_disease_chapters.py",
            "src/genome_report_hotfix11.py",
            "tests/test_report_disease_chapters.py",
            "docs/DISEASE_CHAPTER_TAXONOMY_RU.md",
        ]:
            self.assertIn(path, manifest["deployment_paths"])


if __name__ == "__main__":
    unittest.main()
