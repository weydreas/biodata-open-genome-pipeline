from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import genome_report_hotfix10 as gr10  # noqa: E402


class ReportVersioningTests(unittest.TestCase):
    def test_existing_any_package_member_advances_version(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            v, main, ref, manifest = gr10.allocate_report_package(root, "Dante_genome_report_COMMERCIAL")
            self.assertEqual(v, 1)
            self.assertEqual(main.name, "Dante_genome_report_COMMERCIAL.pdf")
            self.assertEqual(ref.name, "Dante_genome_report_COMMERCIAL_COMMON_MAJOR_REFERENCE.pdf")
            self.assertEqual(manifest.name, "Dante_genome_report_COMMERCIAL.manifest.json")

            main.write_bytes(b"old pdf")
            v, main, ref, manifest = gr10.allocate_report_package(root, "Dante_genome_report_COMMERCIAL")
            self.assertEqual(v, 2)
            self.assertEqual(main.name, "Dante_genome_report_COMMERCIAL_v2.pdf")
            self.assertEqual(ref.name, "Dante_genome_report_COMMERCIAL_COMMON_MAJOR_REFERENCE_v2.pdf")
            self.assertEqual(manifest.name, "Dante_genome_report_COMMERCIAL_v2.manifest.json")

            ref.write_bytes(b"partial v2")
            v, main, ref, manifest = gr10.allocate_report_package(root, "Dante_genome_report_COMMERCIAL")
            self.assertEqual(v, 3)
            self.assertEqual(main.name, "Dante_genome_report_COMMERCIAL_v3.pdf")

    def test_versioned_manifest_records_no_overwrite_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            main = root / "R_v2.pdf"
            ref = root / "R_COMMON_MAJOR_REFERENCE_v2.pdf"
            manifest_path = root / "R_v2.manifest.json"
            main.write_bytes(b"%PDF test main")
            ref.write_bytes(b"%PDF test ref")
            state = {"pipeline_release": "9.0.0", "rules_schema": "dante-knowledge-rules-v8"}
            gr10._write_versioned_manifest(
                output_path=manifest_path,
                version=2,
                mode="COMMERCIAL",
                patient="Test",
                state=state,
                root=root,
                source_files={},
                research_layer_opened=False,
                main_grouped={},
                common_grouped={},
                category_order=[],
                main_pdf=main,
                reference_pdf=ref,
            )
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["publication_contract"], "dante-report-publication-v10")
            self.assertEqual(manifest["report_version"], 2)
            self.assertEqual(manifest["version_suffix"], "_v2")
            self.assertFalse(manifest["versioning_policy"]["overwrite_existing_outputs"])
            self.assertTrue(manifest["versioning_policy"]["partial_package_consumes_version"])
            self.assertEqual(manifest["output"]["name"], "R_v2.pdf")
            self.assertEqual(manifest["reference_appendix"]["name"], "R_COMMON_MAJOR_REFERENCE_v2.pdf")

    def test_powershell_ready_validator_uses_latest_versioned_manifest(self) -> None:
        source = (SRC / "report_runtime.ps1").read_text(encoding="utf-8-sig")
        self.assertIn("function Get-LatestGenomeReportManifest", source)
        self.assertIn("_v(?<version>[0-9]+)", source)
        self.assertIn("dante-report-publication-v10", source)
        self.assertIn("$manifest.output.name", source)
        self.assertIn("$manifest.reference_appendix.name", source)
        self.assertIn("Existing report files are never overwritten", source)

    def test_deployment_contains_versioned_output_module(self) -> None:
        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        self.assertIn("src/genome_report_hotfix10.py", manifest["deployment_paths"])
        self.assertIn("src/report_runtime.ps1", manifest["deployment_paths"])


if __name__ == "__main__":
    unittest.main()
