"""Clickable table-of-contents publication layer over report hotfix15.

This hotfix changes only PDF publication/navigation. Stage 24 scientific findings
and the Clinical Risk Evidence artifact remain unchanged. ReportLab performs a
multi-pass build so the table of contents contains final page numbers, clickable
internal links, matching PDF outline/bookmark entries, and a compact return link
from each TOC-listed report heading back to the contents page.
"""
from __future__ import annotations

from contextlib import contextmanager
import json
import re
import sys
from pathlib import Path

import clinical_risk_evidence_hotfix1 as risk_impl
import genome_report_hotfix15 as h15

# Keep the context-isolated Clinical Risk Evidence extractor as the production
# implementation inherited from hotfix15.
h15.risk = risk_impl

PUBLICATION_CONTRACT = "dante-report-clickable-toc-backlinks-v17"
TOC_SCHEMA = "dante-report-clickable-toc-v2"
TOC_HOME_DESTINATION = "dante_toc_home"
TOC_BACKLINK_LABEL = "вернуться"
GENERATION_MODE = h15.GENERATION_MODE
_ORIGINAL_RENDERER = h15.render_report_pdf
_ORIGINAL_MANIFEST_FACTORY = h15._manifest_with_risk


def _plain(flowable) -> str:
    preserved = getattr(flowable, "_dante_toc_label", None)
    if preserved is not None:
        return str(preserved).strip()
    if hasattr(flowable, "getPlainText"):
        try:
            return str(flowable.getPlainText()).strip()
        except Exception:
            pass
    return str(getattr(flowable, "text", "") or "").strip()


def _toc_level(text: str) -> int:
    """Disease chapters are children of the detailed registry; other H1s are top-level."""
    return 1 if re.match(r"^\d+\.\s+", text or "") else 0


def _is_report_h1(flowable) -> bool:
    style = getattr(flowable, "style", None)
    return bool(style and getattr(style, "name", "") == "DanteH1")


def _decorate_heading_backlink(flowable, platypus):
    """Append a small internal return link without polluting TOC/bookmark labels."""
    if not _is_report_h1(flowable):
        return flowable
    label = _plain(flowable)
    raw_text = str(getattr(flowable, "text", "") or label)
    backlink = (
        f'&nbsp;&nbsp;<font size="6.2">['
        f'<link href="#{TOC_HOME_DESTINATION}">{TOC_BACKLINK_LABEL}</link>]'
        f'</font>'
    )
    decorated = platypus.Paragraph(raw_text + backlink, flowable.style)
    decorated._dante_toc_label = label
    return decorated


def _inject_toc(story: list, platypus, *, mm) -> list:
    """Insert TOC after title/meta block and add return links to TOC-listed H1s."""
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.pdfbase import pdfmetrics
    from reportlab.platypus.tableofcontents import TableOfContents

    # Production renderer registers DanteSans/DanteSansBold before it creates the
    # document. Standalone/self-test callers may not; use built-ins only there.
    registered = set(pdfmetrics.getRegisteredFontNames())
    regular_font = "DanteSans" if "DanteSans" in registered else "Helvetica"
    bold_font = "DanteSansBold" if "DanteSansBold" in registered else "Helvetica-Bold"

    decorated_story = [_decorate_heading_backlink(item, platypus) for item in story]

    toc = TableOfContents()
    toc.levelStyles = [
        ParagraphStyle(
            "DanteTOCLevel0",
            fontName=regular_font,
            fontSize=8.4,
            leading=10.8,
            leftIndent=0,
            firstLineIndent=0,
            spaceBefore=1.2,
            spaceAfter=1.2,
        ),
        ParagraphStyle(
            "DanteTOCLevel1",
            fontName=regular_font,
            fontSize=7.6,
            leading=9.7,
            leftIndent=6 * mm,
            firstLineIndent=0,
            spaceBefore=0.4,
            spaceAfter=0.4,
        ),
    ]
    toc.dotsMinLevel = 0
    title_style = ParagraphStyle(
        "DanteTOCTitle",
        fontName=bold_font,
        fontSize=12.2,
        leading=15,
        spaceBefore=3,
        spaceAfter=5,
    )
    note_style = ParagraphStyle(
        "DanteTOCNote",
        fontName=regular_font,
        fontSize=7.1,
        leading=9.2,
        spaceAfter=4,
    )
    block = [
        platypus.Paragraph(f'<a name="{TOC_HOME_DESTINATION}"/>Содержание', title_style),
        platypus.Paragraph(
            "Названия разделов и номера страниц являются внутренними ссылками PDF. "
            "У заголовков разделов есть компактная ссылка для возврата к содержанию.",
            note_style,
        ),
        toc,
        platypus.PageBreak(),
    ]
    insert_at = next((i for i, item in enumerate(decorated_story) if _is_report_h1(item)), len(decorated_story))
    return list(decorated_story[:insert_at]) + block + list(decorated_story[insert_at:])


@contextmanager
def clickable_toc_documents():
    """Replace only SimpleDocTemplate during rendering with a TOC-aware equivalent."""
    try:
        from reportlab import platypus
        from reportlab.lib.units import mm
        from reportlab.pdfgen import canvas
        from reportlab.platypus import BaseDocTemplate, Frame, PageTemplate, Paragraph
    except ModuleNotFoundError:
        yield
        return

    original_simple = platypus.SimpleDocTemplate

    class DanteTOCSimpleDocTemplate(original_simple):
        """SimpleDocTemplate-compatible multi-pass renderer with clickable TOC/bookmarks."""

        def beforeDocument(self):
            try:
                super().beforeDocument()
            except AttributeError:
                pass
            # multiBuild repeats the document; stable per-pass sequence => stable keys.
            self._dante_heading_seq = 0

        def afterFlowable(self, flowable):
            try:
                super().afterFlowable(flowable)
            except AttributeError:
                pass
            if not isinstance(flowable, Paragraph) or not _is_report_h1(flowable):
                return
            text = _plain(flowable)
            if not text or text == "Содержание":
                return
            self._dante_heading_seq = int(getattr(self, "_dante_heading_seq", 0)) + 1
            key = f"dante_toc_{self._dante_heading_seq:03d}"
            level = _toc_level(text)
            # bookmarkPage provides the destination used by TableOfContents links.
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(text, key, level=level, closed=False)
            self.notify("TOCEntry", (level, text, self.page, key))

        def build(self, flowables, onFirstPage=None, onLaterPages=None, canvasmaker=canvas.Canvas):
            # BaseDocTemplate.multiBuild invokes self.build for each internal pass.
            # The inner call must perform one normal BaseDocTemplate build instead
            # of recursively starting another multiBuild.
            if getattr(self, "_dante_toc_inner_build", False):
                return BaseDocTemplate.build(self, flowables, canvasmaker=canvasmaker)

            from reportlab.platypus.doctemplate import _doNothing

            if onFirstPage is None:
                onFirstPage = _doNothing
            if onLaterPages is None:
                onLaterPages = _doNothing

            self._calc()
            frame = Frame(self.leftMargin, self.bottomMargin, self.width, self.height, id="normal")
            self.pageTemplates = []
            self.addPageTemplates(
                [
                    PageTemplate(id="First", frames=frame, onPage=onFirstPage, pagesize=self.pagesize),
                    PageTemplate(id="Later", frames=frame, onPage=onLaterPages, pagesize=self.pagesize),
                ]
            )
            if onFirstPage is _doNothing and hasattr(self, "onFirstPage"):
                self.pageTemplates[0].beforeDrawPage = self.onFirstPage
            if onLaterPages is _doNothing and hasattr(self, "onLaterPages"):
                self.pageTemplates[1].beforeDrawPage = self.onLaterPages

            story = _inject_toc(list(flowables), platypus, mm=mm)
            self._dante_toc_inner_build = True
            try:
                return BaseDocTemplate.multiBuild(self, story, canvasmaker=canvasmaker)
            finally:
                self._dante_toc_inner_build = False

    platypus.SimpleDocTemplate = DanteTOCSimpleDocTemplate
    try:
        yield
    finally:
        platypus.SimpleDocTemplate = original_simple


def render_report_pdf(output, **kwargs):
    with clickable_toc_documents():
        return _ORIGINAL_RENDERER(output, **kwargs)


def _manifest_with_toc(original_factory):
    """Extend hotfix15 manifest provenance and force old PDFs publication-stale."""
    def factory(original):
        wrapped = original_factory(original)

        def writer(**kwargs):
            wrapped(**kwargs)
            output_path = Path(kwargs["output_path"])
            manifest = json.loads(output_path.read_text(encoding="utf-8"))
            manifest["publication_contract"] = PUBLICATION_CONTRACT
            manifest["generation_mode"] = GENERATION_MODE
            manifest["table_of_contents"] = {
                "schema": TOC_SCHEMA,
                "position": "after_title_block_before_report_sections",
                "page_numbers": "final_reportlab_multibuild_pages",
                "clickable_internal_links": True,
                "pdf_outline_bookmarks": True,
                "included_heading_style": "DanteH1",
                "disease_chapters_as_level_1": True,
                "heading_backlinks_to_contents": True,
                "backlink_destination": TOC_HOME_DESTINATION,
                "backlink_label": TOC_BACKLINK_LABEL,
            }
            output_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

        return writer

    return factory


def self_test() -> int:
    assert _toc_level("1. Онкология — 3 событий") == 1
    assert _toc_level("Knowledge sources и ограничения лицензий") == 0
    try:
        from reportlab import platypus
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
    except ModuleNotFoundError:
        print("SELF-TEST PASS | clickable TOC + backlinks static contract (ReportLab unavailable)")
        return 0

    import tempfile

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="DanteH1", parent=styles["Heading1"], fontName="Helvetica-Bold"))
    original_heading = platypus.Paragraph("1. Онкология — 1 событий", styles["DanteH1"])
    injected = _inject_toc([original_heading], platypus, mm=mm)
    decorated = next(item for item in injected if _is_report_h1(item))
    assert _plain(decorated) == "1. Онкология — 1 событий"
    assert TOC_BACKLINK_LABEL in decorated.getPlainText()
    toc_title = next(item for item in injected if getattr(getattr(item, "style", None), "name", "") == "DanteTOCTitle")
    assert TOC_HOME_DESTINATION in str(getattr(toc_title, "text", ""))

    with tempfile.TemporaryDirectory() as directory:
        pdf = Path(directory) / "toc-smoke.pdf"
        story = [
            platypus.Paragraph("Dante report", styles["Title"]),
            platypus.Paragraph("Юридический статус", styles["DanteH1"]),
            platypus.Paragraph("Body", styles["BodyText"]),
            platypus.PageBreak(),
            platypus.Paragraph("1. Онкология — 1 событий", styles["DanteH1"]),
            platypus.Paragraph("Body", styles["BodyText"]),
        ]
        with clickable_toc_documents():
            doc = platypus.SimpleDocTemplate(str(pdf), pagesize=A4)
            doc.build(story)
        payload = pdf.read_bytes()
        assert len(payload) > 1000
        assert b"/Outlines" in payload
        assert b"/Annots" in payload
    print("SELF-TEST PASS | final-page clickable TOC + heading return links + outline bookmarks")
    return 0


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "self-test":
        return self_test()

    original_renderer = h15.render_report_pdf
    original_manifest_factory = h15._manifest_with_risk
    original_contract = h15.PUBLICATION_CONTRACT
    h15.render_report_pdf = globals().get("render_report_pdf", render_report_pdf)
    h15._manifest_with_risk = _manifest_with_toc(original_manifest_factory)
    h15.PUBLICATION_CONTRACT = PUBLICATION_CONTRACT
    try:
        return h15.main()
    finally:
        h15.render_report_pdf = original_renderer
        h15._manifest_with_risk = original_manifest_factory
        h15.PUBLICATION_CONTRACT = original_contract


for _name in dir(h15):
    if _name.startswith("__") or _name in {"main", "render_report_pdf"}:
        continue
    globals().setdefault(_name, getattr(h15, _name))


if __name__ == "__main__":
    raise SystemExit(main())
