from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

import sys
sys.path.insert(0, str(SRC))

import genome_report_hotfix16 as gr16  # noqa: E402


class ClickableTocReportTests(unittest.TestCase):
    def test_toc_hierarchy_contract(self) -> None:
        self.assertEqual(gr16._toc_level("Юридический статус и границы использования"), 0)
        self.assertEqual(gr16._toc_level("Детальный реестр Stage-24 находок по клиническим разделам"), 0)
        self.assertEqual(gr16._toc_level("1. Онкология — 4 событий"), 1)
        self.assertEqual(gr16._toc_level("13. Врожденные пороки — 5 событий"), 1)

    def test_stable_entrypoint_routes_through_hotfix16(self) -> None:
        source = (SRC / "genome_report.py").read_text(encoding="utf-8")
        self.assertIn("genome_report_hotfix16", source)
        self.assertIn("genome_report_hotfix15", source)
        self.assertIn("clinical_risk_evidence_hotfix1", source)

    def test_runtime_requires_clickable_toc_and_backlink_contract(self) -> None:
        source = (SRC / "clinical_risk_runtime.ps1").read_text(encoding="utf-8")
        self.assertIn("dante-report-clickable-toc-backlinks-v17", source)
        self.assertIn("dante-report-clickable-toc-v2", source)
        self.assertIn("clickable_internal_links", source)
        self.assertIn("pdf_outline_bookmarks", source)
        self.assertIn("heading_backlinks_to_contents", source)
        self.assertIn("dante_toc_home", source)

    def test_windows_ready_contract_does_not_compare_localized_backlink_label(self) -> None:
        source = (SRC / "clinical_risk_runtime.ps1").read_text(encoding="utf-8")
        self.assertNotIn('backlink_label -ne "вернуться"', source)
        self.assertIn("non-ASCII source literals unsuitable as machine-readable READY contract values", source)
        self.assertIn('backlink_destination -ne "dante_toc_home"', source)

    @unittest.skipUnless(importlib.util.find_spec("reportlab") is not None, "reportlab not installed")
    def test_heading_backlink_preserves_clean_toc_label(self) -> None:
        from reportlab import platypus
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet

        styles = getSampleStyleSheet()
        styles.add(ParagraphStyle(name="DanteH1", parent=styles["Heading1"], fontName="Helvetica-Bold"))
        heading = platypus.Paragraph("1. Онкология — 4 событий", styles["DanteH1"])
        decorated = gr16._decorate_heading_backlink(heading, platypus)
        self.assertIn(gr16.TOC_BACKLINK_LABEL, decorated.getPlainText())
        self.assertEqual(gr16._plain(decorated), "1. Онкология — 4 событий")
        self.assertIn(f'#{gr16.TOC_HOME_DESTINATION}', str(decorated.text))

    def test_manifest_records_toc_semantics(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "report.manifest.json"

            def original_factory(original):
                def wrapped(**kwargs):
                    original(**kwargs)
                return wrapped

            def original(**kwargs):
                kwargs["output_path"].write_text(
                    json.dumps({"publication_contract": "old", "generation_mode": "old"}),
                    encoding="utf-8",
                )

            writer = gr16._manifest_with_toc(original_factory)(original)
            writer(output_path=output)
            manifest = json.loads(output.read_text(encoding="utf-8"))
            toc = manifest["table_of_contents"]
            self.assertEqual(manifest["publication_contract"], gr16.PUBLICATION_CONTRACT)
            self.assertEqual(toc["schema"], gr16.TOC_SCHEMA)
            self.assertTrue(toc["clickable_internal_links"])
            self.assertTrue(toc["pdf_outline_bookmarks"])
            self.assertTrue(toc["heading_backlinks_to_contents"])
            self.assertEqual(toc["backlink_destination"], gr16.TOC_HOME_DESTINATION)
            self.assertEqual(toc["backlink_label"], gr16.TOC_BACKLINK_LABEL)
            self.assertEqual(toc["page_numbers"], "final_reportlab_multibuild_pages")

    @unittest.skipUnless(importlib.util.find_spec("reportlab") is not None, "reportlab not installed")
    def test_reportlab_smoke_contains_outline_and_bidirectional_link_annotations(self) -> None:
        self.assertEqual(gr16.self_test(), 0)

    def test_deployment_manifest_tracks_toc_layer(self) -> None:
        manifest = json.loads((ROOT / "release" / "pipeline_manifest.json").read_text(encoding="utf-8"))
        self.assertIn("src/genome_report_hotfix16.py", manifest["deployment_paths"])
        self.assertIn("tests/test_report_clickable_toc.py", manifest["deployment_paths"])
        self.assertIn("docs/REPORT_CLICKABLE_TOC_CONTRACT_RU.md", manifest["deployment_paths"])


if __name__ == "__main__":
    unittest.main()
