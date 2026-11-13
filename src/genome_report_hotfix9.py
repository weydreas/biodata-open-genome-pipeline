#!/usr/bin/env python3
"""Hotfix9 publication layer for deterministic Stage-24 genome PDFs.

Scientific Stage-24 rules are unchanged. This module changes only publication:
- events with Max AF > 0.70 are removed from the main PDF and placed in a
  separate reference appendix;
- the same rule is applied to research overlays using their open-event AF context;
- event headings include Stage-24 gene(s) and the patient's actual allele(s)
  reconstructed from REF/ALT + genotype.

No AI/API calls are made.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import genome_report_base_v2 as base

REPORT_SCHEMA = "dante-genome-report-v2"
COMMON_MAJOR_AF_THRESHOLD = 0.70
COMMON_MAJOR_RULE = "max_af > 0.70"

# Re-export stable helpers/contracts used by regression tests and callers.
validate_usage = base.validate_usage
group_reason_events = base.group_reason_events
reason_counts = base.reason_counts
OPEN_CATEGORY_ORDER = base.OPEN_CATEGORY_ORDER
RESEARCH_CATEGORY_ORDER = base.RESEARCH_CATEGORY_ORDER
OPEN_CATEGORY_META = base.OPEN_CATEGORY_META
RESEARCH_CATEGORY_META = base.RESEARCH_CATEGORY_META
USAGE_SCHEMA = base.USAGE_SCHEMA
STAGE24_SCHEMA = base.STAGE24_SCHEMA
RULES_SCHEMA = base.RULES_SCHEMA
REPORTLAB_VERSION = base.REPORTLAB_VERSION


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _float_values(value: str) -> list[float]:
    if not value or str(value).strip() in {"", ".", "-"}:
        return []
    text = str(value).strip()
    try:
        return [float(text)]
    except ValueError:
        values: list[float] = []
        for token in re.split(r"[;,|\s]+", text):
            try:
                values.append(float(token))
            except ValueError:
                continue
        return values


def event_max_af(event: dict, open_context: list[dict[str, str]] | None = None) -> float | None:
    values: list[float] = []
    for row in event.get("rows", []):
        values.extend(_float_values(base.first(row, "max_af")))
    for row in open_context or []:
        values.extend(_float_values(base.first(row, "max_af")))
    return max(values) if values else None


def split_common_major_events(
    grouped: dict[str, list[dict]],
    open_by_event: dict[str, list[dict[str, str]]] | None = None,
) -> tuple[dict[str, list[dict]], dict[str, list[dict]]]:
    """Split publication events without changing Stage-24 source data.

    Strict user-facing publication rule: if the available Max AF for an event is
    greater than 0.70, the event is not highlighted in the main report. It is
    retained in a separate reference appendix. Unknown AF stays in the main report.
    """
    main: dict[str, list[dict]] = {}
    common: dict[str, list[dict]] = {}
    for reason, events in grouped.items():
        for event in events:
            context = (open_by_event or {}).get(event.get("event_key", ""), [])
            af = event_max_af(event, context)
            target = common if af is not None and af > COMMON_MAJOR_AF_THRESHOLD else main
            target.setdefault(reason, []).append(event)
    return main, common


def _patient_allele(event: dict) -> str:
    rows = event.get("rows") or []
    if not rows:
        return ""
    row = rows[0]
    genotype = base.first(row, "genotype") or str(event.get("genotype", "") or "").strip()
    ref = base.first(row, "ref") or str(event.get("ref", "") or "").strip()
    alt = base.first(row, "alt") or str(event.get("alt", "") or "").strip()
    if not genotype:
        return ""
    if not ref or not alt:
        return genotype

    alleles = [ref] + [x.strip() for x in alt.split(",")]
    tokens = re.split(r"([/|])", genotype)
    rendered: list[str] = []
    for token in tokens:
        if token in {"/", "|"}:
            rendered.append(token)
            continue
        token = token.strip()
        if token == ".":
            rendered.append("?")
            continue
        try:
            idx = int(token)
        except ValueError:
            rendered.append(token)
            continue
        rendered.append(alleles[idx] if 0 <= idx < len(alleles) else token)
    return "".join(rendered)


def event_heading_label(event: dict) -> str:
    label = base.event_label(event)
    genes = base.unique(base.first(row, "gene") for row in event.get("rows", []))
    if genes:
        shown = genes[:4]
        suffix = f", +{len(genes) - 4}" if len(genes) > 4 else ""
        label += " (" + ", ".join(shown) + suffix + ")"
    allele = _patient_allele(event)
    if allele:
        label += f" | аллель пациента {allele}"
    return label


def _count_rows(grouped: dict[str, list[dict]]) -> int:
    return sum(len(event.get("rows", [])) for events in grouped.values() for event in events)


def _count_events(grouped: dict[str, list[dict]]) -> int:
    return len({event.get("event_key") for events in grouped.values() for event in events})


def _font_paths() -> tuple[Path, Path]:
    regular_candidates = [
        Path(r"C:\Windows\Fonts\arial.ttf"),
        Path(r"C:\Windows\Fonts\segoeui.ttf"),
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
        Path("/usr/share/fonts/dejavu/DejaVuSans.ttf"),
    ]
    bold_candidates = [
        Path(r"C:\Windows\Fonts\arialbd.ttf"),
        Path(r"C:\Windows\Fonts\segoeuib.ttf"),
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        Path("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf"),
    ]
    regular = next((p for p in regular_candidates if p.is_file()), None)
    bold = next((p for p in bold_candidates if p.is_file()), None)
    if not regular or not bold:
        raise RuntimeError("No local Cyrillic-capable report font found")
    return regular, bold


def render_report_pdf(
    output: Path,
    *,
    title: str,
    patient: str,
    mode: str,
    state: dict,
    usage: dict,
    stage_rows: list[dict],
    grouped: dict[str, list[dict]],
    category_meta: dict[str, tuple[str, str]],
    category_order: list[str],
    excluded_grouped: dict[str, list[dict]],
    reference_appendix_name: str,
    is_reference_appendix: bool = False,
    open_by_event: dict[str, list[dict[str, str]]] | None = None,
) -> None:
    try:
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_CENTER
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont
        from reportlab.platypus import KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
    except ImportError as exc:
        raise RuntimeError(f"reportlab is required for PDF generation (recommended {REPORTLAB_VERSION})") from exc

    regular, bold = _font_paths()
    pdfmetrics.registerFont(TTFont("DanteSans", str(regular)))
    pdfmetrics.registerFont(TTFont("DanteSansBold", str(bold)))

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="DanteTitle", parent=styles["Title"], fontName="DanteSansBold", fontSize=17, leading=20, alignment=TA_CENTER, spaceAfter=8))
    styles.add(ParagraphStyle(name="DanteH1", parent=styles["Heading1"], fontName="DanteSansBold", fontSize=11, leading=14, spaceBefore=7, spaceAfter=4))
    styles.add(ParagraphStyle(name="DanteH2", parent=styles["Heading2"], fontName="DanteSansBold", fontSize=8.7, leading=11, spaceBefore=5, spaceAfter=2))
    styles.add(ParagraphStyle(name="DanteBody", parent=styles["BodyText"], fontName="DanteSans", fontSize=8.3, leading=11, spaceAfter=3))
    styles.add(ParagraphStyle(name="DanteSmall", parent=styles["BodyText"], fontName="DanteSans", fontSize=7.1, leading=9.2, spaceAfter=1.5))
    styles.add(ParagraphStyle(name="DanteTiny", parent=styles["BodyText"], fontName="DanteSans", fontSize=6.2, leading=8, spaceAfter=1))
    styles.add(ParagraphStyle(name="DanteWarn", parent=styles["BodyText"], fontName="DanteSans", fontSize=8.1, leading=10.4, backColor=colors.whitesmoke, borderPadding=6, spaceAfter=5))

    def p(text: str, style: str = "DanteBody") -> Paragraph:
        return Paragraph(text, styles[style])

    def esc(value: object) -> str:
        return html.escape(str(value if value is not None else ""))

    def links_markup(links: list[tuple[str, str]]) -> str:
        return " | ".join(f'<link href="{html.escape(url, quote=True)}">{esc(label)}</link>' for label, url in links)

    mode_is_research = mode == "research"
    legal = base.RESEARCH_LEGAL if mode_is_research else base.COMMERCIAL_LEGAL
    if is_reference_appendix:
        banner = "RESEARCH ONLY / COMMON MAJOR REFERENCE APPENDIX / NO AI" if mode_is_research else "OPEN / COMMON MAJOR REFERENCE APPENDIX / NO AI"
    else:
        banner = "RESEARCH ONLY / PRIMARY REPORT / NO AI" if mode_is_research else "OPEN / COMMERCIAL-SAFE / PRIMARY REPORT / NO AI"

    counts = base.reason_counts(grouped, category_order)
    included_events = _count_events(grouped)
    included_rows = _count_rows(grouped)
    excluded_events = _count_events(excluded_grouped)
    excluded_rows = _count_rows(excluded_grouped)

    story = [
        p(esc(title), "DanteTitle"),
        p(f"Пациент: {esc(patient)}", "DanteBody"),
        p(esc(banner), "DanteWarn"),
        p(
            f"Pipeline: {esc(state.get('pipeline_release','?'))} | Rules: {esc(state.get('rules_schema','?'))} | "
            f"Assembly: {esc(state.get('assembly','?'))} | Report schema: {REPORT_SCHEMA}",
            "DanteSmall",
        ),
        Spacer(1, 3 * mm),
    ]

    story.append(p("Юридический статус и границы использования", "DanteH1"))
    for item in legal:
        story.append(p(esc(item), "DanteBody"))

    story.append(p("Правило публикационной фильтрации", "DanteH1"))
    if is_reference_appendix:
        story.append(p(
            f"Это справочное приложение содержит события, у которых доступный Stage-24 Max AF строго больше {COMMON_MAJOR_AF_THRESHOLD:.2f}. "
            "Они не удаляются из Stage 24 и сохраняются для аудита, но не подсвечиваются в основном отчете как индивидуальные находки. "
            "Термин common/major allele здесь означает только высокую наблюдаемую популяционную частоту; он не является самостоятельным доказательством функциональной нейтральности.",
            "DanteWarn",
        ))
    else:
        story.append(p(
            f"Основной отчет не показывает события с Stage-24 Max AF > {COMMON_MAJOR_AF_THRESHOLD:.2f}. "
            f"Они автоматически вынесены в отдельное справочное приложение {esc(reference_appendix_name)}. "
            "Stage-24 TSV остается неизменным и полным; фильтр применяется только при сборке PDF.",
            "DanteWarn",
        ))

    story.append(p("Сводка", "DanteH1"))
    summary_data = [
        [p("Показатель", "DanteSmall"), p("Значение", "DanteSmall")],
        [p("События в этом PDF", "DanteSmall"), p(str(included_events), "DanteSmall")],
        [p("Строки Stage-24 findings/overlays в этом PDF", "DanteSmall"), p(str(included_rows), "DanteSmall")],
        [p("События в парном common-major документе", "DanteSmall"), p(str(excluded_events), "DanteSmall")],
        [p("Строки в парном common-major документе", "DanteSmall"), p(str(excluded_rows), "DanteSmall")],
        [p("Порог common-major", "DanteSmall"), p(COMMON_MAJOR_RULE, "DanteSmall")],
        [p("AI/API при генерации", "DanteSmall"), p("НЕ ИСПОЛЬЗОВАЛИСЬ", "DanteSmall")],
    ]
    table = Table(summary_data, colWidths=[105 * mm, 69 * mm], repeatRows=1)
    table.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), "DanteSans"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(table)

    story.append(p("Категории находок", "DanteH1"))
    count_data = [[p("Reason code", "DanteSmall"), p("События", "DanteSmall"), p("Evidence rows", "DanteSmall")]]
    for item in counts:
        count_data.append([p(esc(item["reason"]), "DanteSmall"), p(str(item["events"]), "DanteSmall"), p(str(item["rows"]), "DanteSmall")])
    if len(count_data) == 1:
        count_data.append([p("-", "DanteSmall"), p("0", "DanteSmall"), p("0", "DanteSmall")])
    count_table = Table(count_data, colWidths=[105 * mm, 33 * mm, 36 * mm], repeatRows=1)
    count_table.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), "DanteSans"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke),
        ("GRID", (0, 0), (-1, -1), 0.2, colors.lightgrey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    story.append(count_table)

    story.append(PageBreak())
    story.append(p("Детальный реестр Stage-24 находок" if not is_reference_appendix else "Справочный реестр common-major событий", "DanteH1"))
    story.append(p(
        "В заголовке каждой карточки после rsID/координаты показаны ген(ы) из Stage 24 и фактический аллель пациента, восстановленный из REF/ALT и genotype. Конкретные значения ниже получены из Stage-24 TSV и подключенных баз.",
        "DanteBody",
    ))

    ordered_reasons = category_order + sorted(set(grouped) - set(category_order))
    global_index = 0
    for reason in ordered_reasons:
        events = grouped.get(reason, [])
        if not events:
            continue
        title_reason, description = category_meta.get(reason, (reason, "Категория присутствует в Stage 24, но для нее не задано специальное stock-описание."))
        story.append(p(f"{esc(title_reason)} - {len(events)} событий / {sum(len(e['rows']) for e in events)} строк", "DanteH1"))
        story.append(p(esc(description), "DanteWarn"))
        for event in events:
            global_index += 1
            heading = f"{global_index}. {event_heading_label(event)} | {reason}"
            context = (open_by_event or {}).get(event["event_key"]) if mode_is_research else None
            fields = base.event_fields(event, context)
            block = [p(esc(heading), "DanteH2")]
            for key, value in fields:
                if value:
                    block.append(p(f"<b>{esc(key)}:</b> {esc(value)}", "DanteSmall"))
            relations = base.relation_lines(event["rows"], reason)
            if relations:
                label = "Research rows" if mode_is_research else "Gene/disease relations"
                block.append(p(f"<b>{label} ({len(relations)}):</b>", "DanteSmall"))
                for line in relations:
                    block.append(p("- " + esc(line), "DanteTiny"))
            links = base.event_database_links(event)
            if mode_is_research:
                seen_links = set(links)
                for open_row in (open_by_event or {}).get(event["event_key"], []):
                    for link in base.row_database_links(open_row, base.first(open_row, "reason_code")):
                        if link not in seen_links:
                            seen_links.add(link)
                            links.append(link)
            if links:
                block.append(p("<b>Базы:</b> " + links_markup(links), "DanteSmall"))
            sources = base.unique(source for row in event["rows"] for source in base.split_multi(base.first(row, "evidence_sources")))
            if sources:
                block.append(p("<b>Evidence sources из Stage 24:</b> " + esc("; ".join(sources)), "DanteTiny"))
            story.append(KeepTogether(block[:3]))
            story.extend(block[3:])
            story.append(Spacer(1, 1.2 * mm))

    story.append(PageBreak())
    story.append(p("Knowledge sources и ограничения лицензий", "DanteH1"))
    for line in base.source_manifest_lines(usage):
        story.append(p(esc(line), "DanteSmall"))

    story.append(p("Воспроизводимость: pipeline items 1-24", "DanteH1"))
    data = [[p("#", "DanteTiny"), p("Этап", "DanteTiny"), p("Артефакт / provenance", "DanteTiny"), p("Статус", "DanteTiny")]]
    for row in stage_rows:
        data.append([p(str(row["item"]), "DanteTiny"), p(esc(row["name"]), "DanteTiny"), p(esc(row["artifact"]), "DanteTiny"), p(esc(row["status"]), "DanteTiny")])
    stages = Table(data, colWidths=[8 * mm, 49 * mm, 88 * mm, 29 * mm], repeatRows=1, splitByRow=1)
    stages.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), "DanteSans"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke),
        ("GRID", (0, 0), (-1, -1), 0.2, colors.lightgrey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    story.append(stages)

    story.append(p("Ограничения текущего аналитического покрытия", "DanteH1"))
    for item in [
        "Порог Max AF > 0.70 является publication filter, а не новой научной классификацией Stage 24.",
        "Полные CYP/PGx haplotypes, HLA typing, repeat expansions, global phasing и mtDNA heteroplasmy требуют отдельных валидированных модулей.",
        "Для структурных вариантов биологическая значимость и техническая достоверность вызова оцениваются раздельно.",
        "Категории K2/K3/K4 характеризуют gene-level evidence context и не являются автоматической variant-level pathogenicity classification.",
    ]:
        story.append(p(esc(item), "DanteBody"))

    footer_text = banner

    def footer(canvas, doc_):
        canvas.saveState()
        canvas.setFont("DanteSans", 6.2)
        canvas.drawString(14 * mm, 8 * mm, footer_text)
        canvas.drawRightString(A4[0] - 14 * mm, 8 * mm, f"стр. {doc_.page}")
        canvas.restoreState()

    output.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(
        str(output), pagesize=A4, rightMargin=14 * mm, leftMargin=14 * mm,
        topMargin=13 * mm, bottomMargin=13 * mm,
        title=title, author="Biodata / Dante WGS pipeline",
    )
    doc.build(story, onFirstPage=footer, onLaterPages=footer)


def _manifest_filter(main_grouped: dict[str, list[dict]], common_grouped: dict[str, list[dict]]) -> dict:
    return {
        "rule": COMMON_MAJOR_RULE,
        "threshold": COMMON_MAJOR_AF_THRESHOLD,
        "source_field": "max_af",
        "unknown_af_policy": "KEEP_IN_PRIMARY",
        "stage24_mutated": False,
        "primary_unique_events": _count_events(main_grouped),
        "primary_rows": _count_rows(main_grouped),
        "reference_unique_events": _count_events(common_grouped),
        "reference_rows": _count_rows(common_grouped),
    }


def _write_manifest(
    *,
    output_path: Path,
    mode: str,
    patient: str,
    state: dict,
    root: Path,
    source_files: dict,
    research_layer_opened: bool,
    main_grouped: dict[str, list[dict]],
    common_grouped: dict[str, list[dict]],
    category_order: list[str],
    main_pdf: Path,
    reference_pdf: Path,
) -> None:
    manifest = {
        "schema": REPORT_SCHEMA,
        "mode": mode,
        "generation_mode": "DETERMINISTIC_STAGE24_ONLY",
        "ai_used": False,
        "network_used_by_report_generator": False,
        # Backward-compatible package-level meaning: all Stage-24 rows are present
        # across the primary PDF + common-major reference appendix.
        "all_stage24_findings_included": True,
        "all_stage24_findings_accounted_for": True,
        "primary_pdf_filtered": True,
        "patient": patient,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "pipeline_release": state.get("pipeline_release"),
        "rules_schema": state.get("rules_schema"),
        "pipeline_commit": base.pipeline_commit(root),
        "source_files": source_files,
        "research_layer_opened": research_layer_opened,
        "unique_events": _count_events(main_grouped),
        "finding_rows": _count_rows(main_grouped),
        "category_counts": base.reason_counts(main_grouped, category_order),
        "common_major_filter": _manifest_filter(main_grouped, common_grouped),
        "output": {"name": main_pdf.name, "sha256": sha256(main_pdf)},
        "reference_appendix": {"name": reference_pdf.name, "sha256": sha256(reference_pdf)},
    }
    output_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate deterministic split WGS PDFs with common-major publication filter")
    parser.add_argument("mode", choices=("commercial", "research", "self-test"))
    parser.add_argument("--root", type=Path)
    parser.add_argument("--patient-root", type=Path)
    parser.add_argument("--accept-research-restrictions", action="store_true")
    return parser.parse_args()


def self_test() -> int:
    rows = [
        {
            "finding_id": "OF_COMMON", "priority_rank": "4", "event_type": "SNV_INDEL",
            "event_key": "GRCh37:11:17435840:A:G", "chrom": "11", "position": "17435840",
            "ref": "A", "alt": "G", "genotype": "1/1", "zygosity": "HOM_ALT",
            "gene": "ABCC8", "max_af": "0.9834", "resolved_rsid": "rs4148632",
            "reason_code": "K1_EXACT_CLINVAR_CONFLICT",
        },
        {
            "finding_id": "OF_COMMON_2", "priority_rank": "4", "event_type": "SNV_INDEL",
            "event_key": "GRCh37:11:17435840:A:G", "chrom": "11", "position": "17435840",
            "ref": "A", "alt": "G", "genotype": "1/1", "zygosity": "HOM_ALT",
            "gene": "LOC110121471", "max_af": "0.9834", "resolved_rsid": "rs4148632",
            "reason_code": "K1_EXACT_CLINVAR_CONFLICT",
        },
        {
            "finding_id": "OF_RARE", "priority_rank": "2", "event_type": "SNV_INDEL",
            "event_key": "GRCh37:1:100:C:T", "chrom": "1", "position": "100",
            "ref": "C", "alt": "T", "genotype": "0/1", "zygosity": "HET_ALT",
            "gene": "GENE1", "max_af": "0.001", "resolved_rsid": "rs1",
            "reason_code": "K3_RARE_MISSENSE_VALIDATED_GENE",
        },
    ]
    grouped = base.group_reason_events(rows)
    primary, common = split_common_major_events(grouped)
    assert _count_events(primary) == 1
    assert _count_events(common) == 1
    common_event = common["K1_EXACT_CLINVAR_CONFLICT"][0]
    heading = event_heading_label(common_event)
    assert "rs4148632" in heading
    assert "(ABCC8, LOC110121471)" in heading
    assert "аллель пациента G/G" in heading
    assert event_max_af(common_event) == 0.9834
    print("SELF-TEST PASS | Max AF > 0.70 reference split + gene/patient-allele headings + no AI/API")
    return 0


def main() -> int:
    args = parse_args()
    if args.mode == "self-test":
        return self_test()
    if not args.root or not args.patient_root:
        raise RuntimeError("--root and --patient-root are required")

    root = args.root.resolve()
    patient_root = args.patient_root.resolve()
    patient = patient_root.name
    knowledge = patient_root / "Dante_results" / "knowledge_enrichment"
    state_path = knowledge / "Dante_knowledge_enrichment.state.json"
    open_path = knowledge / "Dante_knowledge_findings_open.tsv"
    research_path = knowledge / "Dante_knowledge_findings_research.tsv"
    open_usage_path = knowledge / "Dante_knowledge_usage_manifest.json"
    research_usage_path = knowledge / "Dante_knowledge_usage_research_manifest.json"

    for required in (state_path, open_path, open_usage_path):
        if not required.is_file():
            raise FileNotFoundError(required)
    state = base.read_json(state_path)
    if state.get("status") != "complete" or state.get("schema") != STAGE24_SCHEMA or state.get("rules_schema") != RULES_SCHEMA:
        raise RuntimeError("Stage 24 is not READY under rules-v8")

    out_root = patient_root / "Dante_results" / "reports"
    out_root.mkdir(parents=True, exist_ok=True)
    stage_rows = base.stage_presence(root, patient_root)
    open_rows = base.read_tsv(open_path)
    open_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in open_rows:
        open_by_event[base.first(row, "event_key") or base.first(row, "finding_id")].append(row)

    if args.mode == "commercial":
        usage = base.read_json(open_usage_path)
        base.validate_usage("commercial", usage)
        all_grouped = base.group_reason_events(open_rows)
        main_grouped, common_grouped = split_common_major_events(all_grouped)
        main_pdf = out_root / "Dante_genome_report_COMMERCIAL.pdf"
        reference_pdf = out_root / "Dante_genome_report_COMMERCIAL_COMMON_MAJOR_REFERENCE.pdf"
        render_report_pdf(
            main_pdf,
            title="Dante WGS - итоговый геномный отчет",
            patient=patient, mode="commercial", state=state, usage=usage, stage_rows=stage_rows,
            grouped=main_grouped, category_meta=OPEN_CATEGORY_META, category_order=OPEN_CATEGORY_ORDER,
            excluded_grouped=common_grouped, reference_appendix_name=reference_pdf.name,
        )
        render_report_pdf(
            reference_pdf,
            title="Dante WGS - справочное приложение: common major alleles",
            patient=patient, mode="commercial", state=state, usage=usage, stage_rows=stage_rows,
            grouped=common_grouped, category_meta=OPEN_CATEGORY_META, category_order=OPEN_CATEGORY_ORDER,
            excluded_grouped=main_grouped, reference_appendix_name=main_pdf.name, is_reference_appendix=True,
        )
        source_files = {
            "state": {"name": state_path.name, "sha256": sha256(state_path)},
            "open_findings": {"name": open_path.name, "sha256": sha256(open_path)},
            "usage": {"name": open_usage_path.name, "sha256": sha256(open_usage_path)},
        }
        manifest_path = out_root / "Dante_genome_report_COMMERCIAL.manifest.json"
        _write_manifest(
            output_path=manifest_path, mode="COMMERCIAL", patient=patient, state=state, root=root,
            source_files=source_files, research_layer_opened=False, main_grouped=main_grouped,
            common_grouped=common_grouped, category_order=OPEN_CATEGORY_ORDER,
            main_pdf=main_pdf, reference_pdf=reference_pdf,
        )
        print(f"Created COMMERCIAL primary report: {main_pdf}")
        print(f"Created common-major reference appendix: {reference_pdf}")
        print(
            f"Publication filter {COMMON_MAJOR_RULE}: main={_count_events(main_grouped)} events/{_count_rows(main_grouped)} rows; "
            f"reference={_count_events(common_grouped)} events/{_count_rows(common_grouped)} rows. Stage 24 unchanged; OpenAI/API not used."
        )
        return 0

    if not args.accept_research_restrictions:
        if not sys.stdin.isatty():
            raise RuntimeError("Research mode requires --accept-research-restrictions")
        print("RESEARCH ONLY: non-commercial; not for diagnosis, medical decisions or healthcare-service use.")
        if input("Type RESEARCH to continue: ").strip() != "RESEARCH":
            raise RuntimeError("Research report cancelled")
    if not research_path.is_file() or not research_usage_path.is_file():
        raise FileNotFoundError("Research Stage-24 outputs are missing")

    usage = base.read_json(research_usage_path)
    base.validate_usage("research", usage)
    overlays = base.read_tsv(research_path)
    all_grouped = base.group_reason_events(overlays)
    main_grouped, common_grouped = split_common_major_events(all_grouped, open_by_event)
    main_pdf = out_root / "Dante_genome_report_RESEARCH_ONLY.pdf"
    reference_pdf = out_root / "Dante_genome_report_RESEARCH_ONLY_COMMON_MAJOR_REFERENCE.pdf"
    render_report_pdf(
        main_pdf,
        title="Dante WGS - исследовательский геномный отчет",
        patient=patient, mode="research", state=state, usage=usage, stage_rows=stage_rows,
        grouped=main_grouped, category_meta=RESEARCH_CATEGORY_META, category_order=RESEARCH_CATEGORY_ORDER,
        excluded_grouped=common_grouped, reference_appendix_name=reference_pdf.name, open_by_event=open_by_event,
    )
    render_report_pdf(
        reference_pdf,
        title="Dante WGS - RESEARCH ONLY: common major reference appendix",
        patient=patient, mode="research", state=state, usage=usage, stage_rows=stage_rows,
        grouped=common_grouped, category_meta=RESEARCH_CATEGORY_META, category_order=RESEARCH_CATEGORY_ORDER,
        excluded_grouped=main_grouped, reference_appendix_name=main_pdf.name,
        is_reference_appendix=True, open_by_event=open_by_event,
    )
    source_files = {
        "state": {"name": state_path.name, "sha256": sha256(state_path)},
        "open_findings": {"name": open_path.name, "sha256": sha256(open_path)},
        "research_overlays": {"name": research_path.name, "sha256": sha256(research_path)},
        "usage": {"name": research_usage_path.name, "sha256": sha256(research_usage_path)},
    }
    manifest_path = out_root / "Dante_genome_report_RESEARCH_ONLY.manifest.json"
    _write_manifest(
        output_path=manifest_path, mode="RESEARCH_ONLY", patient=patient, state=state, root=root,
        source_files=source_files, research_layer_opened=True, main_grouped=main_grouped,
        common_grouped=common_grouped, category_order=RESEARCH_CATEGORY_ORDER,
        main_pdf=main_pdf, reference_pdf=reference_pdf,
    )
    print(f"Created RESEARCH ONLY primary report: {main_pdf}")
    print(f"Created research common-major reference appendix: {reference_pdf}")
    print(
        f"Publication filter {COMMON_MAJOR_RULE}: main={_count_events(main_grouped)} events/{_count_rows(main_grouped)} rows; "
        f"reference={_count_events(common_grouped)} events/{_count_rows(common_grouped)} rows. Stage 24 unchanged; OpenAI/API not used."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
