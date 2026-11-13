from __future__ import annotations

import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import genome_report_hotfix12 as gr12  # noqa: E402


class ReportTableLayoutTests(unittest.TestCase):
    def test_layout_context_enables_split_in_row_and_restores_table(self) -> None:
        calls: list[dict] = []

        def original_table(*args, **kwargs):
            calls.append(dict(kwargs))
            return (args, kwargs)

        fake_platypus = types.ModuleType("reportlab.platypus")
        fake_platypus.Table = original_table
        fake_reportlab = types.ModuleType("reportlab")
        fake_reportlab.platypus = fake_platypus

        with mock.patch.dict(sys.modules, {"reportlab": fake_reportlab, "reportlab.platypus": fake_platypus}):
            with gr12.splittable_report_tables():
                self.assertIsNot(fake_platypus.Table, original_table)
                fake_platypus.Table([["a"]], splitByRow=1)
                self.assertTrue(calls[-1]["splitInRow"])
            self.assertIs(fake_platypus.Table, original_table)

    def test_explicit_split_in_row_setting_is_not_overridden(self) -> None:
        calls: list[dict] = []

        def original_table(*args, **kwargs):
            calls.append(dict(kwargs))
            return None

        fake_platypus = types.ModuleType("reportlab.platypus")
        fake_platypus.Table = original_table
        fake_reportlab = types.ModuleType("reportlab")
        fake_reportlab.platypus = fake_platypus
        with mock.patch.dict(sys.modules, {"reportlab": fake_reportlab, "reportlab.platypus": fake_platypus}):
            with gr12.splittable_report_tables():
                fake_platypus.Table([["a"]], splitInRow=0)
        self.assertEqual(calls[-1]["splitInRow"], 0)

    def test_stable_entrypoint_routes_through_hotfix12_and_keeps_ready_contract(self) -> None:
        source = (SRC / "genome_report.py").read_text(encoding="utf-8")
        self.assertIn("genome_report_hotfix12", source)
        self.assertIn("_impl.PUBLICATION_CONTRACT = _impl.h10.PUBLICATION_CONTRACT", source)
        self.assertIn("_impl.h11.PUBLICATION_CONTRACT = _impl.h10.PUBLICATION_CONTRACT", source)
        hotfix = (SRC / "genome_report_hotfix12.py").read_text(encoding="utf-8")
        self.assertIn('kwargs.setdefault("splitInRow", 1)', hotfix)
        self.assertIn("preserving the complete condition list", hotfix)

    def test_deployment_manifest_contains_layout_hotfix_and_test(self) -> None:
        import json

        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        self.assertIn("src/genome_report_hotfix12.py", manifest["deployment_paths"])
        self.assertIn("tests/test_report_table_layout.py", manifest["deployment_paths"])

    @unittest.skipUnless(importlib.util.find_spec("reportlab") is not None, "reportlab not installed in CI runtime")
    def test_reportlab_can_split_a_single_oversized_condition_row(self) -> None:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.platypus import Paragraph, SimpleDocTemplate, Table

        styles = getSampleStyleSheet()
        p = lambda text: Paragraph(text, styles["BodyText"])
        long_conditions = "; ".join(
            f"Condition {i} with a deliberately long phenotype description" for i in range(260)
        )
        rows = [[p("rsID / variant"), p("Gene(s)"), p("Associated conditions")]]
        for i in range(75):
            rows.append([p(f"rs{i}"), p("GENE"), p(long_conditions if i == 12 else "short condition")])

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "oversized-row.pdf"
            table = Table(
                rows,
                colWidths=[80, 80, 330],
                repeatRows=1,
                splitByRow=1,
                splitInRow=1,
            )
            SimpleDocTemplate(str(output), pagesize=A4).build([table])
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
