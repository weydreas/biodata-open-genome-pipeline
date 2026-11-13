from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import patient_manager as pm  # noqa: E402


class PatientMigrationPrivacyTests(unittest.TestCase):
    def test_default_patient_is_synthetic(self) -> None:
        self.assertEqual(pm.DEFAULT_PATIENT, "Demo 1")

    def test_completed_migration_keeps_existing_patient_name_and_data(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            existing = pm.ensure_patient(root, "Demo 2")
            results = existing / "Dante_results"
            results.mkdir()
            (results / "sentinel.txt").write_text("existing result", encoding="utf-8")
            marker = existing / "migration_v8.5.0.json"
            original = {
                "schema": pm.MIGRATION_SCHEMA,
                "status": "complete",
                "patient": "Demo 2",
                "moved": [],
            }
            pm.write_json_atomic(marker, original)
            pm.activate(root, "Demo 2")
            before = marker.read_bytes()
            result = pm.bootstrap_legacy_patient(root)
            self.assertEqual(result, original)
            self.assertEqual(marker.read_bytes(), before)
            self.assertEqual((results / "sentinel.txt").read_text(encoding="utf-8"), "existing result")
            self.assertFalse((root / "patients" / "Demo 1").exists())
            self.assertEqual(json.loads((root / ".active_patient.json").read_text(encoding="utf-8"))["name"], "Demo 2")

    def test_initialized_patient_without_migration_marker_is_not_duplicated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            existing = pm.ensure_patient(root, "Demo 2")
            before = (existing / "patient.json").read_bytes()
            result = pm.bootstrap_legacy_patient(root)
            self.assertEqual(result["status"], "already_initialized")
            self.assertEqual(result["moved"], [])
            self.assertEqual((existing / "patient.json").read_bytes(), before)
            self.assertEqual(pm.list_patients(root), ["Demo 2"])
            self.assertFalse((root / "patients" / "Demo 1").exists())

    def test_initial_legacy_migration_targets_demo_1_without_recomputation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = root / "Dante_results"
            work.mkdir()
            (work / "sentinel.txt").write_text("ready", encoding="utf-8")
            result = pm.bootstrap_legacy_patient(root)
            dest = root / "patients" / "Demo 1" / "Dante_results" / "sentinel.txt"
            self.assertEqual(result["patient"], "Demo 1")
            self.assertTrue(dest.is_file())
            self.assertEqual(dest.read_text(encoding="utf-8"), "ready")
            self.assertEqual(pm.bootstrap_legacy_patient(root), result)

    def test_multiple_migration_markers_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("Demo 2", "Demo 3"):
                patient = pm.ensure_patient(root, name)
                pm.write_json_atomic(patient / "migration_v8.5.0.json", {"status": "complete", "moved": []})
            with self.assertRaisesRegex(RuntimeError, "Multiple legacy migration markers"):
                pm.bootstrap_legacy_patient(root)


if __name__ == "__main__":
    unittest.main()
