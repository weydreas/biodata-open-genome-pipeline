#!/usr/bin/env python3
"""Disease-chaptered deterministic publication layer for Stage-24 genome reports.

Builds on hotfix10 output versioning and hotfix9 common-major filtering. Stage-24
scientific rules and source files are unchanged. Publication is reorganized into
clinically readable disease chapters using a deterministic no-AI taxonomy over the
existing ClinVar conditions / disease relations.
"""
from __future__ import annotations

import html
import json
import sys
from collections import defaultdict
from pathlib import Path

from genome_report_hotfix10 import *  # noqa: F401,F403
import genome_report_hotfix10 as h10
import genome_report_hotfix9 as h9
import genome_report_base_v2 as base
import report_disease_chapters as taxonomy

PUBLICATION_CONTRACT = "dante-report-publication-v11"
CHAPTER_ASSIGNMENT_RULE = (
    "condition-level deterministic keyword taxonomy; one primary chapter per event by "
    "classified-condition vote, stable chapter order as tie-breaker; secondary chapters preserved"
)


def _reason_rank(reason: str, order: list[str]) -> int:
    try:
        return order.index(reason)
    except ValueError:
        return len(order) + 100


def merge_reason_events(grouped: dict[str, list[dict]], category_order: list[str]) -> list[dict]:
    """Collapse repeated publication cards for the same Stage-24 event without losing rows."""
    merged: dict[str, dict] = {}
    seen_rows: dict[str, set[tuple[str, str, str]]] = defaultdict(set)
    for reason, events in grouped.items():
        for source in events:
            key = str(source.get("event_key") or base.event_label(source))
            if key not in merged:
                event = dict(source)
                event["rows"] = []
                event["reason_codes"] = []
                merged[key] = event
            event = merged[key]
            if reason not in event["reason_codes"]:
                event["reason_codes"].append(reason)
            for row in source.get("rows", []):
                row_key = (
                    base.first(row, "finding_id", "overlay_id") or "",
                    base.first(row, "reason_code") or reason,
                    json.dumps(row, ensure_ascii=False, sort_keys=True),
                )
                if row_key in seen_rows[key]:
                    continue
                seen_rows[key].add(row_key)
                event["rows"].append(row)
            try:
                event["priority_rank"] = min(int(event.get("priority_rank", 99)), int(source.get("priority_rank", 99)))
            except (TypeError, ValueError):
                pass

    result: list[dict] = []
    for event in merged.values():
        event["reason_codes"] = sorted(event["reason_codes"], key=lambda x: (_reason_rank(x, category_order), x))
        event["reason_code"] = event["reason_codes"][0] if event["reason_codes"] else "UNCLASSIFIED"
        event["rows"].sort(key=lambda row: (_reason_rank(base.first(row, "reason_code") or event["reason_code"], category_order), base.first(row, "finding_id", "overlay_id")))
        result.append(event)
    result.sort(key=lambda e: (
        int(e.get("priority_rank", 99)) if str(e.get("priority_rank", "")).isdigit() else 99,
        _reason_rank(e.get("reason_code", ""), category_order),
        str(e.get("chrom", "")), str(e.get("position", "")), str(e.get("event_key", "")),
    ))
    return result


def _raw_condition_values(event: dict, open_context: list[dict[str, str]] | None = None) -> list[str]:
    values: list[str] = []
    for row in list(event.get("rows", [])) + list(open_context or []):
        for key in ("clinvar_conditions", "disease_name", "condition_name", "condition"):
            value = base.first(row, key)
            if value:
                values.append(value)
    return values


def _fallback_display_conditions(values: list[str]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        for raw in taxonomy.split_conditions(value):
            normalized = taxonomy.normalize_condition(raw)
            if normalized in {"not provided", "not specified", "unspecified", "see cases"}:
                continue
            shown = raw.replace("_", " ")
            key = shown.lower()
            if key not in seen:
                seen.add(key)
                result.append(shown)
    return result


def event_chapter_info(event: dict, open_context: list[dict[str, str]] | None = None) -> dict[str, object]:
    values = _raw_condition_values(event, open_context)
    info = taxonomy.assign_conditions(values)
    display_conditions = list(info["conditions"])
    if not display_conditions:
        display_conditions = _fallback_display_conditions(values)
    info["display_conditions"] = display_conditions
    return info


def group_events_by_chapter(grouped: dict[str, list[dict]], category_order: list[str], open_by_event: dict[str, list[dict[str, str]]] | None = None) -> tuple[dict[str, list[dict]], dict[str, dict[str, object]]]:
    chaptered: dict[str, list[dict]] = defaultdict(list)
    metadata: dict[str, dict[str, object]] = {}
    for event in merge_reason_events(grouped, category_order):
        context = (open_by_event or {}).get(event.get("event_key", ""), [])
        info = event_chapter_info(event, context)
        metadata[event.get("event_key", base.event_label(event))] = info
        chaptered[str(info["primary"])].append(event)
    return dict(chaptered), metadata


def chapter_counts(grouped: dict[str, list[dict]], category_order: list[str], open_by_event: dict[str, list[dict[str, str]]] | None = None) -> list[dict[str, object]]:
    chaptered, _ = group_events_by_chapter(grouped, category_order, open_by_event)
    result: list[dict[str, object]] = []
    for code in taxonomy.CHAPTER_ORDER:
        events = chaptered.get(code, [])
        if events:
            result.append({"code": code, "title": taxonomy.chapter_title(code), "events": len(events), "evidence_rows": sum(len(e.get("rows", [])) for e in events)})
    return result


def _augment_fields(fields: list[tuple[str, str]], info: dict[str, object]) -> list[tuple[str, str]]:
    primary = taxonomy.chapter_title(str(info["primary"]))
    secondary = "; ".join(taxonomy.chapter_title(str(code)) for code in info.get("secondary", []))
    result: list[tuple[str, str]] = []
    inserted = False
    for key, value in fields:
        result.append((key, value))
        if key == "Гены":
            result.append(("Клинический раздел", primary))
            if secondary:
                result.append(("Дополнительные разделы", secondary))
            inserted = True
    if not inserted:
        result.insert(0, ("Клинический раздел", primary))
        if secondary:
            result.insert(1, ("Дополнительные разделы", secondary))
    return result


def _event_relations(event: dict) -> list[str]:
    lines: list[str] = []
    seen: set[str] = set()
    reason_codes = event.get("reason_codes") or [event.get("reason_code", "UNCLASSIFIED")]
    for reason in reason_codes:
        rows = [row for row in event.get("rows", []) if (base.first(row, "reason_code") or reason) == reason]
        if not rows:
            rows = event.get("rows", [])
        for line in base.relation_lines(rows, reason):
            rendered = f"[{reason}] {line}" if len(reason_codes) > 1 else line
            if rendered not in seen:
                seen.add(rendered)
                lines.append(rendered)
    return lines


def render_report_pdf(output: Path, *, title: str, patient: str, mode: str, state: dict, usage: dict, stage_rows: list[dict], grouped: dict[str, list[dict]], category_meta: dict[str, tuple[str, str]], category_order: list[str], excluded_grouped: dict[str, list[dict]], reference_appendix_name: str, is_reference_appendix: bool = False, open_by_event: dict[str, list[dict[str, str]]] | None = None) -> None:
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
        raise RuntimeError(f"reportlab is required for PDF generation (recommended {h9.REPORTLAB_VERSION})") from exc

    regular, bold = h9._font_paths()
    pdfmetrics.registerFont(TTFont("DanteSans", str(regular)))
    pdfmetrics.registerFont(TTFont("DanteSansBold", str(bold)))
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="DanteTitle", parent=styles["Title"], fontName="DanteSansBold", fontSize=17, leading=20, alignment=TA_CENTER, spaceAfter=8))
    styles.add(ParagraphStyle(name="DanteH1", parent=styles["Heading1"], fontName="DanteSansBold", fontSize=11.4, leading=14.2, spaceBefore=7, spaceAfter=4))
    styles.add(ParagraphStyle(name="DanteH2", parent=styles["Heading2"], fontName="DanteSansBold", fontSize=8.8, leading=11, spaceBefore=5, spaceAfter=2))
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
    banner = ("RESEARCH ONLY / COMMON MAJOR REFERENCE APPENDIX / NO AI" if mode_is_research else "OPEN / COMMON MAJOR REFERENCE APPENDIX / NO AI") if is_reference_appendix else ("RESEARCH ONLY / PRIMARY REPORT / NO AI" if mode_is_research else "OPEN / COMMERCIAL-SAFE / PRIMARY REPORT / NO AI")
    included_events, included_rows = h9._count_events(grouped), h9._count_rows(grouped)
    excluded_events, excluded_rows = h9._count_events(excluded_grouped), h9._count_rows(excluded_grouped)
    reason_counts = base.reason_counts(grouped, category_order)
    chaptered, chapter_meta = group_events_by_chapter(grouped, category_order, open_by_event)

    story = [p(esc(title), "DanteTitle"), p(f"Пациент: {esc(patient)}", "DanteBody"), p(esc(banner), "DanteWarn"), p(f"Pipeline: {esc(state.get('pipeline_release','?'))} | Rules: {esc(state.get('rules_schema','?'))} | Assembly: {esc(state.get('assembly','?'))} | Report schema: {h9.REPORT_SCHEMA}", "DanteSmall"), Spacer(1, 3 * mm), p("Юридический статус и границы использования", "DanteH1")]
    for item in legal:
        story.append(p(esc(item), "DanteBody"))
    story.append(p("Правило публикационной фильтрации", "DanteH1"))
    if is_reference_appendix:
        story.append(p(f"Это справочное приложение содержит события с доступным Stage-24 Max AF > {h9.COMMON_MAJOR_AF_THRESHOLD:.2f}. Они сохранены для аудита, но не подсвечиваются в основном отчете как индивидуальные находки.", "DanteWarn"))
    else:
        story.append(p(f"Основной отчет не показывает события с Stage-24 Max AF > {h9.COMMON_MAJOR_AF_THRESHOLD:.2f}; они вынесены в {esc(reference_appendix_name)}. Stage-24 TSV остается неизменным; фильтр применяется только при публикации PDF.", "DanteWarn"))

    story.append(p("Клиническая структура отчета", "DanteH1"))
    story.append(p("События сгруппированы без ИИ по ClinVar conditions и Disease relations, уже присутствующим в Stage 24. Каждый вариант печатается один раз в основном клиническом разделе; если его состояния относятся к нескольким системам, дополнительные разделы сохраняются внутри карточки. Разные reason codes Stage 24 для одного event объединяются в одной карточке без потери evidence rows.", "DanteWarn"))
    summary_data = [[p("Показатель", "DanteSmall"), p("Значение", "DanteSmall")], [p("События в этом PDF", "DanteSmall"), p(str(included_events), "DanteSmall")], [p("Строки Stage-24 findings/overlays", "DanteSmall"), p(str(included_rows), "DanteSmall")], [p("События в парном common-major документе", "DanteSmall"), p(str(excluded_events), "DanteSmall")], [p("Строки в парном common-major документе", "DanteSmall"), p(str(excluded_rows), "DanteSmall")], [p("Disease chapter schema", "DanteSmall"), p(taxonomy.TAXONOMY_SCHEMA, "DanteSmall")], [p("AI/API при генерации", "DanteSmall"), p("НЕ ИСПОЛЬЗОВАЛИСЬ", "DanteSmall")]]
    table = Table(summary_data, colWidths=[105 * mm, 69 * mm], repeatRows=1)
    table.setStyle(TableStyle([("FONTNAME", (0,0), (-1,-1), "DanteSans"), ("BACKGROUND", (0,0), (-1,0), colors.whitesmoke), ("GRID", (0,0), (-1,-1), 0.25, colors.lightgrey), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 3), ("RIGHTPADDING", (0,0), (-1,-1), 3)]))
    story.append(table)

    story.append(p("Разделы заболеваний", "DanteH1"))
    chdata = [[p("Раздел", "DanteSmall"), p("События", "DanteSmall"), p("Evidence rows", "DanteSmall")]]
    for code in taxonomy.CHAPTER_ORDER:
        events = chaptered.get(code, [])
        if events:
            chdata.append([p(esc(taxonomy.chapter_title(code)), "DanteSmall"), p(str(len(events)), "DanteSmall"), p(str(sum(len(e.get('rows', [])) for e in events)), "DanteSmall")])
    cht = Table(chdata, colWidths=[105 * mm, 33 * mm, 36 * mm], repeatRows=1, splitByRow=1)
    cht.setStyle(TableStyle([("FONTNAME", (0,0), (-1,-1), "DanteSans"), ("BACKGROUND", (0,0), (-1,0), colors.whitesmoke), ("GRID", (0,0), (-1,-1), 0.25, colors.lightgrey), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 3), ("RIGHTPADDING", (0,0), (-1,-1), 3)]))
    story.append(cht)

    story.append(p("Категории evidence Stage 24", "DanteH1"))
    rcdata = [[p("Reason code", "DanteSmall"), p("События", "DanteSmall"), p("Evidence rows", "DanteSmall")]]
    for item in reason_counts:
        rcdata.append([p(esc(item["reason"]), "DanteSmall"), p(str(item["events"]), "DanteSmall"), p(str(item["rows"]), "DanteSmall")])
    rct = Table(rcdata, colWidths=[105 * mm, 33 * mm, 36 * mm], repeatRows=1, splitByRow=1)
    rct.setStyle(TableStyle([("FONTNAME", (0,0), (-1,-1), "DanteSans"), ("BACKGROUND", (0,0), (-1,0), colors.whitesmoke), ("GRID", (0,0), (-1,-1), 0.25, colors.lightgrey), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 3), ("RIGHTPADDING", (0,0), (-1,-1), 3)]))
    story.append(rct)

    story.append(PageBreak())
    story.append(p("Детальный реестр Stage-24 находок по клиническим разделам", "DanteH1"))
    global_index = 0
    chapter_number = 0
    for code in taxonomy.CHAPTER_ORDER:
        events = chaptered.get(code, [])
        if not events:
            continue
        chapter_number += 1
        if chapter_number > 1:
            story.append(PageBreak())
        story.append(p(f"{chapter_number}. {esc(taxonomy.chapter_title(code))} — {len(events)} событий", "DanteH1"))
        story.append(p(esc(taxonomy.chapter_description(code)), "DanteWarn"))
        story.append(p("Краткий перечень вариантов", "DanteH2"))
        index_data = [[p("rsID / вариант", "DanteTiny"), p("Ген(ы)", "DanteTiny"), p("Связанные состояния", "DanteTiny")]]
        for event in events:
            info = chapter_meta[event.get("event_key", base.event_label(event))]
            genes = base.unique(base.first(row, "gene") for row in event.get("rows", []))
            conditions = list(info.get("display_conditions", []))
            index_data.append([p(esc(base.event_label(event)), "DanteTiny"), p(esc(", ".join(genes) or "-"), "DanteTiny"), p(esc("; ".join(conditions) or "-"), "DanteTiny")])
        idx_table = Table(index_data, colWidths=[38 * mm, 42 * mm, 94 * mm], repeatRows=1, splitByRow=1)
        idx_table.setStyle(TableStyle([("FONTNAME", (0,0), (-1,-1), "DanteSans"), ("BACKGROUND", (0,0), (-1,0), colors.whitesmoke), ("GRID", (0,0), (-1,-1), 0.2, colors.lightgrey), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 2), ("RIGHTPADDING", (0,0), (-1,-1), 2), ("TOPPADDING", (0,0), (-1,-1), 2), ("BOTTOMPADDING", (0,0), (-1,-1), 2)]))
        story.append(idx_table)
        story.append(Spacer(1, 2 * mm))

        by_reason: dict[str, list[dict]] = defaultdict(list)
        for event in events:
            by_reason[event.get("reason_code", "UNCLASSIFIED")].append(event)
        reason_order = category_order + sorted(set(by_reason) - set(category_order))
        for reason in reason_order:
            reason_events = by_reason.get(reason, [])
            if not reason_events:
                continue
            reason_title, reason_description = category_meta.get(reason, (reason, "Категория Stage 24 без отдельного stock-описания."))
            story.append(p(f"{esc(reason_title)} — {len(reason_events)} событий", "DanteH2"))
            story.append(p(esc(reason_description), "DanteSmall"))
            for event in reason_events:
                global_index += 1
                info = chapter_meta[event.get("event_key", base.event_label(event))]
                reason_codes = event.get("reason_codes") or [reason]
                heading = f"{global_index}. {h9.event_heading_label(event)} | {'; '.join(reason_codes)}"
                block = [p(esc(heading), "DanteH2")]
                context = (open_by_event or {}).get(event.get("event_key", ""), [])
                fields = _augment_fields(base.event_fields(event, context if mode_is_research else None), info)
                for key, value in fields:
                    if value:
                        block.append(p(f"<b>{esc(key)}:</b> {esc(value)}", "DanteSmall"))
                relations = _event_relations(event)
                if relations:
                    label = "Research rows / relations" if mode_is_research else "Gene/disease relations"
                    block.append(p(f"<b>{label} ({len(relations)}):</b>", "DanteSmall"))
                    for line in relations:
                        block.append(p("- " + esc(line), "DanteTiny"))
                links = base.event_database_links(event)
                if mode_is_research:
                    seen_links = set(links)
                    for open_row in context:
                        for link in base.row_database_links(open_row, base.first(open_row, "reason_code")):
                            if link not in seen_links:
                                seen_links.add(link)
                                links.append(link)
                if links:
                    block.append(p("<b>Базы:</b> " + links_markup(links), "DanteSmall"))
                sources = base.unique(source for row in event.get("rows", []) for source in base.split_multi(base.first(row, "evidence_sources")))
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
    stages.setStyle(TableStyle([("FONTNAME", (0,0), (-1,-1), "DanteSans"), ("BACKGROUND", (0,0), (-1,0), colors.whitesmoke), ("GRID", (0,0), (-1,-1), 0.2, colors.lightgrey), ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 2), ("RIGHTPADDING", (0,0), (-1,-1), 2)]))
    story.append(stages)
    story.append(p("Ограничения интерпретации", "DanteH1"))
    for item in ["Клинические разделы являются навигационной структурой отчета, а не диагнозами и не новой классификацией патогенности.", "Один вариант может иметь состояния из нескольких разделов; он печатается один раз в primary-разделе, а дополнительные разделы явно перечисляются в карточке.", "Автоматическая группировка использует только текст Stage-24 conditions/disease relations и поэтому может оставлять неспецифичные записи в разделе 'Недостаточно данных для клинической группировки'.", "Текущий pipeline не заменяет отдельные валидированные анализы PRS/GWAS, PGx star-alleles, HLA, repeat expansions, phasing и mtDNA heteroplasmy, если эти модули отдельно не реализованы."]:
        story.append(p(esc(item), "DanteBody"))

    footer_text = banner
    def footer(canvas, doc_):
        canvas.saveState(); canvas.setFont("DanteSans", 6.2); canvas.drawString(14 * mm, 8 * mm, footer_text); canvas.drawRightString(A4[0] - 14 * mm, 8 * mm, f"стр. {doc_.page}"); canvas.restoreState()
    output.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(str(output), pagesize=A4, rightMargin=14 * mm, leftMargin=14 * mm, topMargin=13 * mm, bottomMargin=13 * mm, title=title, author="Biodata / deterministic Stage-24 renderer")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)


def _write_manifest_v11(*, output_path: Path, version: int, mode: str, patient: str, state: dict, root: Path, source_files: dict, research_layer_opened: bool, main_grouped: dict[str, list[dict]], common_grouped: dict[str, list[dict]], category_order: list[str], main_pdf: Path, reference_pdf: Path, open_by_event: dict[str, list[dict[str, str]]] | None = None) -> None:
    h10._write_versioned_manifest(output_path=output_path, version=version, mode=mode, patient=patient, state=state, root=root, source_files=source_files, research_layer_opened=research_layer_opened, main_grouped=main_grouped, common_grouped=common_grouped, category_order=category_order, main_pdf=main_pdf, reference_pdf=reference_pdf)
    manifest = json.loads(output_path.read_text(encoding="utf-8"))
    manifest["publication_contract"] = PUBLICATION_CONTRACT
    manifest["disease_chapter_schema"] = taxonomy.TAXONOMY_SCHEMA
    manifest["disease_chapter_assignment"] = CHAPTER_ASSIGNMENT_RULE
    manifest["one_primary_chapter_per_event"] = True
    manifest["secondary_chapters_preserved"] = True
    manifest["duplicate_event_cards_across_reason_codes"] = False
    manifest["chapter_counts"] = chapter_counts(main_grouped, category_order, open_by_event)
    manifest["reference_appendix_chapter_counts"] = chapter_counts(common_grouped, category_order, open_by_event)
    output_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def self_test() -> int:
    taxonomy.self_test(); h10.self_test()
    rows = [
        {"event_key":"E1", "reason_code":"K1_EXACT_CLINVAR_CONFLICT", "finding_id":"F1", "resolved_rsid":"rs1", "gene":"PCSK9", "clinvar_conditions":"Familial_hypercholesterolemia", "priority_rank":"4"},
        {"event_key":"E1", "reason_code":"K4_CLINGEN_ACTIONABLE", "finding_id":"F2", "resolved_rsid":"rs1", "gene":"PCSK9", "disease_name":"Familial_hypercholesterolemia", "priority_rank":"4"},
        {"event_key":"E2", "reason_code":"K1_EXACT_CLINVAR_CONFLICT", "finding_id":"F3", "resolved_rsid":"rs2", "gene":"CHEK2", "clinvar_conditions":"Breast_cancer|Colorectal_cancer", "priority_rank":"4"},
    ]
    grouped = base.group_reason_events(rows)
    events = merge_reason_events(grouped, base.OPEN_CATEGORY_ORDER)
    assert len(events) == 2
    pcsk9 = next(e for e in events if e.get("resolved_rsid") == "rs1")
    assert set(pcsk9["reason_codes"]) == {"K1_EXACT_CLINVAR_CONFLICT", "K4_CLINGEN_ACTIONABLE"}
    chaptered, meta = group_events_by_chapter(grouped, base.OPEN_CATEGORY_ORDER)
    assert any(e.get("resolved_rsid") == "rs1" for e in chaptered["D02_CARDIOVASCULAR"])
    assert any(e.get("resolved_rsid") == "rs2" for e in chaptered["D01_ONCOLOGY"])
    assert meta["E1"]["primary"] == "D02_CARDIOVASCULAR"
    print("SELF-TEST PASS | disease-chaptered publication + reason-code event dedup")
    return 0


def main() -> int:
    args = h9.parse_args()
    if args.mode == "self-test": return self_test()
    if not args.root or not args.patient_root: raise RuntimeError("--root and --patient-root are required")
    root, patient_root = args.root.resolve(), args.patient_root.resolve(); patient = patient_root.name
    knowledge = patient_root / "Dante_results" / "knowledge_enrichment"
    state_path = knowledge / "Dante_knowledge_enrichment.state.json"; open_path = knowledge / "Dante_knowledge_findings_open.tsv"; research_path = knowledge / "Dante_knowledge_findings_research.tsv"; open_usage_path = knowledge / "Dante_knowledge_usage_manifest.json"; research_usage_path = knowledge / "Dante_knowledge_usage_research_manifest.json"
    for required in (state_path, open_path, open_usage_path):
        if not required.is_file(): raise FileNotFoundError(required)
    state = base.read_json(state_path)
    if state.get("status") != "complete" or state.get("schema") != h9.STAGE24_SCHEMA or state.get("rules_schema") != h9.RULES_SCHEMA: raise RuntimeError("Stage 24 is not READY under rules-v8")
    out_root = patient_root / "Dante_results" / "reports"; out_root.mkdir(parents=True, exist_ok=True)
    stage_rows = base.stage_presence(root, patient_root); open_rows = base.read_tsv(open_path)
    open_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in open_rows: open_by_event[base.first(row, "event_key") or base.first(row, "finding_id")].append(row)

    if args.mode == "commercial":
        usage = base.read_json(open_usage_path); base.validate_usage("commercial", usage)
        all_grouped = base.group_reason_events(open_rows); main_grouped, common_grouped = h9.split_common_major_events(all_grouped)
        version, main_pdf, reference_pdf, manifest_path = h10.allocate_report_package(out_root, "Dante_genome_report_COMMERCIAL")
        render_report_pdf(main_pdf, title="Dante WGS - итоговый геномный отчет", patient=patient, mode="commercial", state=state, usage=usage, stage_rows=stage_rows, grouped=main_grouped, category_meta=h9.OPEN_CATEGORY_META, category_order=h9.OPEN_CATEGORY_ORDER, excluded_grouped=common_grouped, reference_appendix_name=reference_pdf.name)
        render_report_pdf(reference_pdf, title="Dante WGS - справочное приложение: common major alleles", patient=patient, mode="commercial", state=state, usage=usage, stage_rows=stage_rows, grouped=common_grouped, category_meta=h9.OPEN_CATEGORY_META, category_order=h9.OPEN_CATEGORY_ORDER, excluded_grouped=main_grouped, reference_appendix_name=main_pdf.name, is_reference_appendix=True)
        source_files = {"state":{"name":state_path.name,"sha256":h9.sha256(state_path)}, "open_findings":{"name":open_path.name,"sha256":h9.sha256(open_path)}, "usage":{"name":open_usage_path.name,"sha256":h9.sha256(open_usage_path)}}
        _write_manifest_v11(output_path=manifest_path, version=version, mode="COMMERCIAL", patient=patient, state=state, root=root, source_files=source_files, research_layer_opened=False, main_grouped=main_grouped, common_grouped=common_grouped, category_order=h9.OPEN_CATEGORY_ORDER, main_pdf=main_pdf, reference_pdf=reference_pdf)
        print(f"Created COMMERCIAL disease-chaptered report package v{version}: {main_pdf.name}"); print(f"Created common-major reference appendix: {reference_pdf.name}"); print(f"Created manifest: {manifest_path.name}"); print(f"Clinical structure: {taxonomy.TAXONOMY_SCHEMA}; no AI/API; Stage 24 unchanged."); return 0

    if not args.accept_research_restrictions:
        if not sys.stdin.isatty(): raise RuntimeError("Research mode requires --accept-research-restrictions")
        print("RESEARCH ONLY: non-commercial; not for diagnosis, medical decisions or healthcare-service use.")
        if input("Type RESEARCH to continue: ").strip() != "RESEARCH": raise RuntimeError("Research report cancelled")
    if not research_path.is_file() or not research_usage_path.is_file(): raise FileNotFoundError("Research Stage-24 outputs are missing")
    usage = base.read_json(research_usage_path); base.validate_usage("research", usage); overlays = base.read_tsv(research_path)
    all_grouped = base.group_reason_events(overlays); main_grouped, common_grouped = h9.split_common_major_events(all_grouped, open_by_event)
    version, main_pdf, reference_pdf, manifest_path = h10.allocate_report_package(out_root, "Dante_genome_report_RESEARCH_ONLY")
    render_report_pdf(main_pdf, title="Dante WGS - исследовательский геномный отчет", patient=patient, mode="research", state=state, usage=usage, stage_rows=stage_rows, grouped=main_grouped, category_meta=h9.RESEARCH_CATEGORY_META, category_order=h9.RESEARCH_CATEGORY_ORDER, excluded_grouped=common_grouped, reference_appendix_name=reference_pdf.name, open_by_event=open_by_event)
    render_report_pdf(reference_pdf, title="Dante WGS - RESEARCH ONLY: common major reference appendix", patient=patient, mode="research", state=state, usage=usage, stage_rows=stage_rows, grouped=common_grouped, category_meta=h9.RESEARCH_CATEGORY_META, category_order=h9.RESEARCH_CATEGORY_ORDER, excluded_grouped=main_grouped, reference_appendix_name=main_pdf.name, is_reference_appendix=True, open_by_event=open_by_event)
    source_files = {"state":{"name":state_path.name,"sha256":h9.sha256(state_path)}, "open_findings":{"name":open_path.name,"sha256":h9.sha256(open_path)}, "research_overlays":{"name":research_path.name,"sha256":h9.sha256(research_path)}, "usage":{"name":research_usage_path.name,"sha256":h9.sha256(research_usage_path)}}
    _write_manifest_v11(output_path=manifest_path, version=version, mode="RESEARCH_ONLY", patient=patient, state=state, root=root, source_files=source_files, research_layer_opened=True, main_grouped=main_grouped, common_grouped=common_grouped, category_order=h9.RESEARCH_CATEGORY_ORDER, main_pdf=main_pdf, reference_pdf=reference_pdf, open_by_event=open_by_event)
    print(f"Created RESEARCH ONLY disease-chaptered report package v{version}: {main_pdf.name}"); print(f"Created research common-major reference appendix: {reference_pdf.name}"); print(f"Created manifest: {manifest_path.name}"); print(f"Clinical structure: {taxonomy.TAXONOMY_SCHEMA}; no AI/API; Stage 24 unchanged."); return 0


if __name__ == "__main__": raise SystemExit(main())
