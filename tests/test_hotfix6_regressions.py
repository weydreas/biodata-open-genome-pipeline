from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path
from unittest import mock
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))
import knowledge_enrichment as ke  # noqa: E402


class Hotfix6RegressionTests(unittest.TestCase):
    def test_rules_schema_v8(self) -> None:
        self.assertEqual(ke.RULES_SCHEMA, "dante-knowledge-rules-v8")

    def test_alphamissense_exact_match_materializes_research_reason(self) -> None:
        event = ke.PatientEvent("SNV_INDEL", "v1", {"CATALOG"}, "1", 100, 100, ref="A", alt="G")
        event.raw_context["alpha_missense_target"] = "yes"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            work = root / "work"
            work.mkdir()
            data = root / "alpha.tsv.gz"
            data.write_text("placeholder", encoding="ascii")

            def fake_tabix(image, source, region_file, output, mounts, include_header=False, attempts=3):
                output.write_text("chr1\t100\tA\tG\tENST\tP\t0.93\tpathogenic\n", encoding="utf-8")

            with mock.patch.object(ke, "source_artifact", return_value=data), \
                 mock.patch.object(ke, "docker_tabix", side_effect=fake_tabix):
                counts = ke.apply_alphamissense({event.event_key: event}, root, {}, "image", work)

        self.assertEqual(counts["matched"], 1)
        self.assertEqual(counts["research_overlays_materialized"], 1)
        self.assertIn("K8_ALPHAMISSENSE_SUPPORT", event.research_reasons)
        rows = ke.research_overlays({event.event_key: event}, [], type("K", (), {"profile": lambda self, gene: None})())
        self.assertEqual(rows[0]["evidence_sources"], "AlphaMissense")

    def test_cross_chrom_bnd_breakpoints_are_preserved_for_giab(self) -> None:
        event = ke.PatientEvent("SV", "x", {"ANNOTSV_SV"}, "1", 100, 101)
        event.raw_context = {
            "sv_type": "TRA",
            "structural_members": [
                {"chrom": "1", "start": 100, "end": 101, "partner_chrom": "14", "partner_position": 200},
                {"chrom": "14", "start": 200, "end": 201, "partner_chrom": "1", "partner_position": 100},
            ],
        }
        points = ke.structural_breakpoint_positions(event)
        self.assertIn(("1", 100), points)
        self.assertIn(("14", 200), points)
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / "mate.bed.gz"
            with gzip.open(bed, "wt", encoding="utf-8") as handle:
                handle.write("chr14\t199\t201\n")
            ke.apply_giab({"x": event}, [("MATE_TEST", bed)])
        self.assertIn("MATE_TEST", event.giab_flags)

    def test_structural_patient_priority_requires_independent_confirmation(self) -> None:
        event = ke.PatientEvent("SV", "x", {"ANNOTSV_SV"}, "7", 100, 200)
        event.raw_context = {"sv_type": "DEL", "technical_confidence": "SINGLE_CALLER_STRONG"}
        rank, relevance = ke.patient_priority("K5_DOSAGE_SENSITIVE_CNV", event, "FLNC", None, "UNKNOWN", 2)
        self.assertEqual(rank, 5)
        self.assertIn("TECHNICAL_SINGLE_CALLER", relevance)
        event.raw_context["technical_confidence"] = "CROSS_CALLER_CONFIRMED"
        rank, _ = ke.patient_priority("K5_DOSAGE_SENSITIVE_CNV", event, "FLNC", None, "UNKNOWN", 2)
        self.assertEqual(rank, 2)

    def test_menu_statuses_use_fixed_column_dot_leaders(self) -> None:
        wrapper = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8-sig")
        self.assertIn("$statusColumn = 78", wrapper)
        self.assertIn('"." * ($padding - 2)', wrapper)
        self.assertIn('Write-Host $leader -NoNewline -ForegroundColor DarkGray', wrapper)
        self.assertIn('Write-Host "[$Status]" -ForegroundColor $color', wrapper)

    def test_interactive_launcher_returns_to_menu_after_success(self) -> None:
        runtime_or_wrapper = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8-sig")
        launcher = (ROOT / "start.cmd").read_text(encoding="utf-8-sig")
        self.assertIn('"exit" { $exitCode = 10 }', runtime_or_wrapper)
        self.assertIn(":PIPELINE_MENU", launcher)
        self.assertIn('-PerformanceProfile max -NoPause', launcher)
        self.assertIn('if "%PIPELINE_EXIT%"=="10"', launcher)
        self.assertIn('if not "%PIPELINE_EXIT%"=="0"', launcher)
        self.assertIn("goto PIPELINE_MENU", launcher)
        self.assertNotIn('powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%DANTE_PS1%" -PerformanceProfile max\nif errorlevel 1 pause', launcher)


if __name__ == "__main__":
    unittest.main()
