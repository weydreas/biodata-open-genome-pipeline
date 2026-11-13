from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import patient_folder_discovery as discovery  # noqa: E402
import patient_manager as pm  # noqa: E402


class ManualPatientFolderDiscoveryTests(unittest.TestCase):
    def test_bare_patient_folder_is_adopted_non_destructively(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            patient = root / "patients" / "Demo 2"
            patient.mkdir(parents=True)
            marker = patient / "keep_me.txt"
            marker.write_text("do not move", encoding="utf-8")
            root_source = patient / "sample.raw.vcf.gz"
            root_source.write_bytes(b"source stays where user put it")

            adopted, ignored = discovery.discover_manual_patient_folders(root)

            self.assertEqual(adopted, ["Demo 2"])
            self.assertEqual(ignored, [])
            self.assertTrue((patient / "patient.json").is_file())
            self.assertTrue((patient / "input").is_dir())
            self.assertTrue((patient / "input" / "reference_panel_sources").is_dir())
            self.assertEqual(marker.read_text(encoding="utf-8"), "do not move")
            self.assertTrue(root_source.is_file())
            self.assertFalse((patient / "input" / root_source.name).exists())
            self.assertIn("Demo 2", pm.list_patients(root))

            meta = json.loads((patient / "patient.json").read_text(encoding="utf-8"))
            self.assertEqual(meta["schema"], pm.SCHEMA)
            self.assertEqual(meta["name"], "Demo 2")

    def test_discovery_is_idempotent_and_does_not_rewrite_existing_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            patient = root / "patients" / "Demo 3"
            patient.mkdir(parents=True)
            meta = patient / "patient.json"
            meta.write_text('{"schema":"custom","name":"Demo 3","sentinel":7}\n', encoding="utf-8")

            adopted, ignored = discovery.discover_manual_patient_folders(root)

            self.assertEqual(adopted, [])
            self.assertEqual(ignored, [])
            self.assertEqual(json.loads(meta.read_text(encoding="utf-8"))["sentinel"], 7)

    def test_hidden_and_invalid_manual_folders_are_not_adopted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            patients = root / "patients"
            (patients / ".cache").mkdir(parents=True)
            # Colon is invalid on Windows; on Linux CI it can exist and exercises the guard.
            (patients / "Bad:Name").mkdir()

            adopted, ignored = discovery.discover_manual_patient_folders(root)

            self.assertEqual(adopted, [])
            self.assertEqual(ignored, ["Bad:Name"])
            self.assertFalse((patients / ".cache" / "patient.json").exists())
            self.assertFalse((patients / "Bad:Name" / "patient.json").exists())

    def test_second_discovery_run_is_noop(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            patient = root / "patients" / "Demo 4"
            patient.mkdir(parents=True)

            first, _ = discovery.discover_manual_patient_folders(root)
            metadata_before = (patient / "patient.json").read_text(encoding="utf-8")
            second, _ = discovery.discover_manual_patient_folders(root)

            self.assertEqual(first, ["Demo 4"])
            self.assertEqual(second, [])
            self.assertEqual((patient / "patient.json").read_text(encoding="utf-8"), metadata_before)

    def test_startup_runs_discovery_before_patient_selection(self) -> None:
        start = (ROOT / "start.cmd").read_text(encoding="utf-8")
        discovery_call = '"%DANTE_PYTHON%" "%DANTE_PATIENT_DISCOVERY%" --root "%DANTE_ROOT%"'
        select_call = '"%DANTE_PYTHON%" "%DANTE_PATIENT_MANAGER%" select --root "%DANTE_ROOT%"'
        self.assertIn('set "DANTE_PATIENT_DISCOVERY=%DANTE_SRC%\\patient_folder_discovery.py"', start)
        self.assertIn(discovery_call, start)
        self.assertIn(select_call, start)
        self.assertLess(start.index(discovery_call), start.index(select_call))

    def test_deployment_manifest_updates_patient_discovery_code(self) -> None:
        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        paths = manifest["deployment_paths"]
        self.assertIn("src/patient_manager.py", paths)
        self.assertIn("src/patient_folder_discovery.py", paths)
        self.assertIn("tests/test_patient_folder_discovery.py", paths)


if __name__ == "__main__":
    unittest.main()
