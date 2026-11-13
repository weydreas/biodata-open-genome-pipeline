#!/usr/bin/env python3
"""Deterministic final PDF reports from Stage 24.

Publication modes are physically separated:
- commercial: reads only Stage-24 open/commercial-safe findings + open usage manifest.
- research: reads open context + Stage-24 research overlays and requires explicit acknowledgement.

No OpenAI/API research is performed. The PDFs are a deterministic presentation of
Stage-24 structured evidence with stock section explanations, database-derived
annotations, source links, provenance identifiers and pipeline-stage references.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from urllib.parse import quote

REPORT_SCHEMA = "dante-genome-report-v2"
STAGE24_SCHEMA = "dante-knowledge-enrichment-v1"
RULES_SCHEMA = "dante-knowledge-rules-v8"
USAGE_SCHEMA = "knowledge-usage-v1"
REPORTLAB_VERSION = "4.2.5"

COMMERCIAL_LEGAL = [
    "Информационно-аналитический отчет по данным полногеномного секвенирования (WGS). Отчет не является медицинской услугой, медицинским заключением, диагнозом, назначением лечения или заменой очной консультации врача.",
    "Отчет автоматически формируется из детерминированных результатов Stage 24. При его создании не используются генеративный ИИ, внешняя интерпретация или передача геномных данных в OpenAI/API.",
    "Классификация вариантов и связи ген-заболевание отражают состояние подключенных баз знаний на дату их установки. Они могут изменяться при обновлении ClinVar, ClinGen, GenCC и других источников.",
    "Клинически значимые решения должны приниматься врачом с учетом анамнеза, фенотипа и иных обследований. Находки, способные повлиять на медицинские решения, целесообразно подтверждать валидированным клинико-диагностическим методом, особенно при ограниченной технической уверенности вызова.",
    "Отсутствие значимой находки в этом отчете не исключает наследственное заболевание или генетический риск: текущий pipeline имеет ограничения по фазированию, повторным экспансиям, HLA, фармакогенетическим гаплотипам, митохондриальной гетероплазмии и другим классам анализа, если они не были отдельно реализованы и валидированы.",
]

RESEARCH_LEGAL = [
    "RESEARCH ONLY / ТОЛЬКО ДЛЯ ИССЛЕДОВАТЕЛЬСКОГО ИСПОЛЬЗОВАНИЯ.",
    "Этот отчет не предназначен для коммерческого использования, диагностики, медицинского решения или оказания медицинской услуги.",
    "Отчет автоматически формируется из Stage-24 research-only overlays. При его создании не используются генеративный ИИ, внешняя интерпретация или передача геномных данных в OpenAI/API.",
    "Он может содержать исследовательские источники и вычислительные предикторы с отдельными лицензионными и/или use-case ограничениями, включая PanelApp, AlphaMissense и CADD. Ограничения конкретных источников перечислены в разделе лицензирования.",
    "Исследовательские предикторы и panel membership не являются клинической классификацией варианта. Они используются только как дополнительный слой гипотез поверх детерминированного open-контекста Stage 24.",
]

PIPELINE_ITEMS = [
    (1, "VEP 116 GRCh37 cache", "annotation_data/vep_rsid/"),
    (2, "Dante 23andMe-style variant export", "Dante_results/Dante_23andMe_style.txt.gz"),
    (3, "AnnotSV public annotations", "annotation_data/annotsv/"),
    (4, "CNV/SV AnnotSV annotation", "Dante_results/Dante_CNV_AnnotSV.tsv; Dante_results/Dante_SV_AnnotSV.tsv"),
    (5, "Stage 1 orchestration shortcut", "items 1-4"),
    (6, "GRCh37 reference FASTA", "annotation_data/vep_rsid/fasta/"),
    (7, "ClinVar GRCh37", "annotation_data/clinvar_grch37/clinvar.vcf.gz"),
    (8, "Normalized SNP/INDEL/residual layers", "Dante_work/Dante_*_normalized.vcf.gz"),
    (9, "Full VEP annotation", "Dante_results/full_vep/"),
    (10, "Unified annotated small-variant catalog", "Dante_results/Dante_small_variants_annotated.tsv.gz"),
    (11, "Deterministic shortlist + genome-wide uncertain coding sweep", "Dante_results/Dante_priority_candidates.tsv"),
    (12, "Stage 2 orchestration shortcut", "items 6-11"),
    (13, "Reference-genotype target panel", "Dante_qc/Dante_reference_targets.tsv"),
    (14, "Merge known array genotypes", "Dante_results/Dante_23andMe_reference_complete.txt.gz"),
    (15, "FASTQ reference-aware target alignment", "Dante_work/Dante_reference_targets.bam"),
    (16, "Strict reference/alternate target calls", "Dante_work/Dante_reference_calls.vcf.gz"),
    (17, "Tiered reference-complete export / FASTQ evidence", "Dante_results/Dante_reference_FASTQ_evidence.tsv.gz"),
    (18, "Stage 3 orchestration shortcut", "items 13-17"),
    (19, "Telomere estimate", "Dante_results/Dante_telomere_report.txt"),
    (20, "Reserved in current v9 menu", "no analytical output"),
    (21, "Automation shortcut: stages 1+2", "items 1-12"),
    (22, "Open/commercial-safe knowledge bundle", "reference_data/open_reference/"),
    (23, "Research-only knowledge bundle", "reference_data/research_reference/"),
    (24, "Full-genome knowledge enrichment", "Dante_results/knowledge_enrichment/"),
]

SOURCE_URLS = {
    "clinvar": ("ClinVar", "https://www.ncbi.nlm.nih.gov/clinvar/"),
    "vep/catalog": ("Ensembl VEP", "https://www.ensembl.org/info/docs/tools/vep/index.html"),
    "vep": ("Ensembl VEP", "https://www.ensembl.org/info/docs/tools/vep/index.html"),
    "clingen": ("ClinGen", "https://clinicalgenome.org/"),
    "clingen actionability": ("ClinGen Actionability", "https://actionability.clinicalgenome.org/"),
    "clingen dosage": ("ClinGen Dosage", "https://clinicalgenome.org/working-groups/dosage-sensitivity-curation/"),
    "gencc": ("GenCC", "https://thegencc.org/"),
    "hgnc": ("HGNC", "https://www.genenames.org/"),
    "mondo": ("Mondo", "https://mondo.monarchinitiative.org/"),
    "giab": ("Genome in a Bottle", "https://www.nist.gov/programs-projects/genome-bottle"),
    "dbsnp": ("dbSNP", "https://www.ncbi.nlm.nih.gov/snp/"),
    "annotsv": ("AnnotSV", "https://lbgi.fr/AnnotSV/"),
    "panelapp": ("PanelApp", "https://panelapp.genomicsengland.co.uk/"),
    "alphamissense": ("AlphaMissense", "https://github.com/google-deepmind/alphamissense"),
    "cadd": ("CADD", "https://cadd.gs.washington.edu/"),
}

OPEN_CATEGORY_META = {
    "K1_EXACT_CLINVAR_PLP": (
        "K1 - ClinVar: точное Pathogenic / Likely pathogenic совпадение",
        "Точный аллель chr:pos:ref:alt найден в ClinVar и имеет агрегированную клиническую значимость Pathogenic/Likely pathogenic без смешения с benign или trait-only категорией. Это сильная аннотация конкретного варианта, но сама по себе не устанавливает диагноз и требует учета зиготности, наследования, фенотипа и качества исходного вызова.",
    ),
    "K1_EXACT_CLINVAR_PLP_COMMON": (
        "K1 - ClinVar P/LP, но распространенный аллель",
        "Точный аллель имеет ClinVar P/LP-аннотацию, однако наблюдаемая популяционная частота не соответствует автоматическому приоритету редкого менделирующего варианта. Категория сохранена для прозрачности, но намеренно понижена в приоритете.",
    ),
    "K1_EXACT_CLINVAR_CONFLICT": (
        "K1 - ClinVar: конфликтующие интерпретации",
        "Точный аллель найден в ClinVar, но агрегированная запись содержит конфликтующие клинические интерпретации. Pipeline фиксирует сам конфликт и параметры записи, но не пытается автоматически разрешить спор между submitter-ами.",
    ),
    "K1_EXACT_CLINVAR_MIXED": (
        "K1 - ClinVar: смешанная клиническая категория",
        "В агрегированной записи ClinVar одновременно присутствуют P/LP и trait/association-type интерпретации. Такая запись не повышается автоматически до чистой P/LP-находки.",
    ),
    "K1_EXACT_CLINVAR_ASSOCIATION": (
        "K1 - ClinVar: association / risk factor / protective / drug response",
        "Точный аллель найден в ClinVar, но запись относится к ассоциации, фактору риска, защитному эффекту, drug response или другой trait-категории, а не к чистой менделирующей патогенности. Эти находки описывают контекст и не равны диагнозу.",
    ),
    "K2_RARE_LOF_VALIDATED_GENE": (
        "K2 - редкий predicted loss-of-function в валидированном disease gene",
        "VEP определяет loss-of-function последствие, популяционная частота поддерживает редкость, а ген имеет сильную/definitive gene-disease связь в open knowledge layer (ClinGen/GenCC). Значимость конкретного варианта зависит от механизма заболевания, транскрипта, зиготности и наследования.",
    ),
    "K2_LOF_VALIDATED_GENE_AF_UNKNOWN": (
        "K2 - predicted loss-of-function в валидированном disease gene, AF неизвестна",
        "Loss-of-function последствие находится в гене с сильной gene-disease валидностью, но популяционная частота варианта не установлена. Неизвестная AF является неопределенностью, а не доказательством редкости.",
    ),
    "K3_RARE_MISSENSE_VALIDATED_GENE": (
        "K3 - редкий missense в валидированном disease gene",
        "VEP определяет missense-последствие, вариант редкий по доступной AF, а ген имеет сильную/definitive gene-disease связь. Нахождение missense в disease gene является основанием для дальнейшей оценки, но не доказывает патогенность конкретной замены.",
    ),
    "K3_MISSENSE_VALIDATED_GENE_AF_UNKNOWN": (
        "K3 - missense в валидированном disease gene, AF неизвестна",
        "Missense-последствие находится в гене с сильной gene-disease валидностью, но популяционная частота варианта не установлена. Категория сохраняется как кандидат с повышенной неопределенностью.",
    ),
    "K4_CLINGEN_ACTIONABLE": (
        "K4 - ClinGen Actionability context",
        "Функциональный вариант расположен в гене/состоянии, для которого ClinGen содержит actionability-контекст. Actionability относится к gene-condition pair и не является доказательством патогенности именно этого варианта.",
    ),
    "K4_CLINGEN_ACTIONABLE_AF_UNKNOWN": (
        "K4 - ClinGen Actionability context, AF неизвестна",
        "Как и K4, но популяционная частота конкретного функционального варианта неизвестна. Категория не повышает вариант до клинически значимого без variant-level доказательств.",
    ),
    "K5_DOSAGE_SENSITIVE_CNV": (
        "K5 - CNV в dosage-sensitive гене/регионе",
        "CNV пересекает ген или регион с evidence по дозовой чувствительности. В отчете отдельно показываются биологический контекст и техническая уверенность CNV-call; single-caller события не приравниваются к независимо подтвержденным.",
    ),
    "K6_HIGH_EVIDENCE_SV_GENE": (
        "K6 - структурный вариант с сильным gene-disease контекстом",
        "Структурный вариант пересекает gene-disease relation с сильной evidence-категорией. Для SV особенно важно различать биологическую значимость затронутого гена и техническую достоверность самого структурного вызова.",
    ),
    "K6_STRUCTURAL_GENE_CANDIDATE": (
        "K6 - структурный gene candidate",
        "Структурный вариант пересекает потенциально релевантный ген/локус, но evidence или техническая уверенность недостаточны для более сильной категории. Такие события предназначены прежде всего для прозрачного реестра и ручной проверки при необходимости.",
    ),
}

RESEARCH_CATEGORY_META = {
    "K7_PANELAPP_RESEARCH": (
        "K7 - PanelApp research overlay",
        "Ген/вариант попадает в PanelApp panel context. Panel membership показывает экспертный disease-panel контекст, но не является variant-level классификацией и не доказывает причинность. Этот слой остается research-only.",
    ),
    "K8_ALPHAMISSENSE_SUPPORT": (
        "K8 - AlphaMissense research support",
        "Для missense-варианта доступен вычислительный AlphaMissense score/class. Это predictor evidence: он может поддерживать гипотезу, но не заменяет клиническую классификацию, segregation, functional evidence или phenotype concordance.",
    ),
    "K8_CADD_SUPPORT_TARGET": (
        "K8 - CADD research support",
        "Для варианта доступен CADD score, отражающий относительный deleteriousness/ranking по интегрированным признакам. Высокий CADD не является самостоятельным доказательством патогенности или медицинской actionability.",
    ),
}

OPEN_CATEGORY_ORDER = list(OPEN_CATEGORY_META)
RESEARCH_CATEGORY_ORDER = list(RESEARCH_CATEGORY_META)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def first(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip() not in {"", "."}:
            return str(value).strip()
    return ""


def int_or(value: str, default: int = 999) -> int:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return default


def split_multi(value: str) -> list[str]:
    if not value or value.strip() in {"", ".", "-"}:
        return []
    return [x.strip() for x in re.split(r"[;,|]", value) if x.strip() and x.strip() not in {".", "-"}]


def unique(values: Iterable[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        value = (value or "").strip()
        if not value or value in {".", "-"} or value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def validate_usage(mode: str, usage: dict) -> list[str]:
    if usage.get("schema") != USAGE_SCHEMA:
        raise RuntimeError(f"Unsupported usage manifest schema: {usage.get('schema')}")
    notices = []
    for item in usage.get("sources_used", []):
        source = str(item.get("source", "UNKNOWN"))
        scope = str(item.get("publication_scope", "UNKNOWN")).upper()
        prohibited = [str(x).lower() for x in item.get("prohibited_uses", [])]
        if mode == "commercial" and scope != "COMMERCIAL_ALLOWED":
            raise RuntimeError(f"COMMERCIAL publication blocked by source {source}: scope={scope}")
        if mode == "research" and scope == "UNKNOWN":
            raise RuntimeError(f"RESEARCH publication blocked by unknown source scope: {source}")
        if prohibited:
            notices.append(f"{source}: " + "; ".join(prohibited))
    if mode == "commercial" and usage.get("scope") != "OPEN_FINDINGS_ONLY":
        raise RuntimeError(f"Commercial report requires OPEN_FINDINGS_ONLY usage manifest, got {usage.get('scope')}")
    return notices


def group_reason_events(rows: list[dict[str, str]]) -> dict[str, list[dict]]:
    grouped: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        reason = first(row, "reason_code") or "UNCLASSIFIED"
        event_key = first(row, "event_key") or first(row, "finding_id", "overlay_id") or "UNKNOWN"
        grouped[reason][event_key].append(row)
    result: dict[str, list[dict]] = {}
    for reason, events in grouped.items():
        items = []
        for event_key, event_rows in events.items():
            ordered = sorted(event_rows, key=lambda r: (
                int_or(first(r, "priority_rank")),
                first(r, "gene"),
                first(r, "disease_name"),
                first(r, "finding_id", "overlay_id"),
            ))
            rep = ordered[0]
            items.append({
                "event_key": event_key,
                "reason_code": reason,
                "priority_rank": min(int_or(first(r, "priority_rank")) for r in ordered),
                "event_type": first(rep, "event_type"),
                "chrom": first(rep, "chrom"),
                "position": first(rep, "position"),
                "end": first(rep, "end"),
                "ref": first(rep, "ref"),
                "alt": first(rep, "alt"),
                "genotype": first(rep, "genotype"),
                "zygosity": first(rep, "zygosity"),
                "resolved_rsid": first(rep, "resolved_rsid", "rsid"),
                "source_layers": ";".join(unique(first(r, "source_layers") for r in ordered)),
                "technical_confidence": first(rep, "technical_confidence"),
                "giab_flags": first(rep, "giab_flags"),
                "rows": ordered,
            })
        result[reason] = sorted(items, key=lambda e: (
            e["priority_rank"], e["chrom"], int_or(e["position"]), e["event_key"]
        ))
    return result


def event_label(event: dict) -> str:
    if event.get("resolved_rsid"):
        return str(event["resolved_rsid"])
    chrom, pos = event.get("chrom", "?"), event.get("position", "?")
    ref, alt = event.get("ref", "?"), event.get("alt", "?")
    if ref not in {"", "."} or alt not in {"", "."}:
        return f"{chrom}:{pos}:{ref or '?'}:{alt or '?'}"
    end = event.get("end")
    return f"{chrom}:{pos}-{end or pos}"


def pipeline_commit(root: Path) -> str:
    repo = root / ".pipeline_repo"
    if not (repo / ".git").exists():
        return "unknown"
    try:
        value = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            text=True, stderr=subprocess.DEVNULL, timeout=5,
        ).strip()
        return value if re.fullmatch(r"[0-9a-f]{40}", value) else "unknown"
    except Exception:
        return "unknown"


def stage_presence(root: Path, patient_root: Path) -> list[dict]:
    rows = []
    for number, name, artifact in PIPELINE_ITEMS:
        if number in {5, 12, 18, 20, 21}:
            status = "META/RESERVED"
        else:
            probes = []
            for token in [x.strip() for x in artifact.split(";")]:
                if token.startswith("Dante_"):
                    token = "Dante_results/" + token
                base = patient_root if token.startswith(("Dante_results/", "Dante_work/", "Dante_qc/")) else root
                path = base / token
                if "*" in token:
                    probes.append(any(base.glob(token)))
                else:
                    probes.append(path.exists())
            status = "PRESENT" if probes and all(probes) else "PARTIAL/NOT CHECKED"
        rows.append({"item": number, "name": name, "artifact": artifact, "status": status})
    return rows


def stage_links_for_event(event: dict, mode: str) -> str:
    items = {24}
    event_type = str(event.get("event_type", ""))
    if event_type in {"SNV_INDEL", "SMALL"}:
        items.update({7, 8, 9, 10})
    if event_type in {"STRUCTURAL", "SV", "CNV"}:
        items.add(4)
    layers = str(event.get("source_layers", ""))
    if "STAGE17" in layers or "REFERENCE" in layers:
        items.add(17)
    items.add(22)
    if mode == "research":
        items.add(23)
    return ", ".join(str(x) for x in sorted(items))


def row_database_links(row: dict[str, str], reason: str) -> list[tuple[str, str]]:
    links: list[tuple[str, str]] = []
    seen = set()

    def add(label: str, url: str) -> None:
        key = (label, url)
        if key not in seen:
            seen.add(key)
            links.append(key)

    rsid = first(row, "resolved_rsid", "rsid")
    if rsid and re.fullmatch(r"rs\d+", rsid, flags=re.I):
        add(f"dbSNP {rsid}", f"https://www.ncbi.nlm.nih.gov/snp/{quote(rsid)}")
        if reason.startswith("K1_") or first(row, "clinvar_significance"):
            add(f"ClinVar {rsid}", f"https://www.ncbi.nlm.nih.gov/clinvar/?term={quote(rsid)}")

    sources = []
    for raw in unique(first(r, "evidence_sources") for r in [row]):
        sources.extend(split_multi(raw))
    if reason.startswith("K1_"):
        sources.extend(["ClinVar", "VEP/CATALOG"])
    if reason.startswith(("K2_", "K3_")):
        sources.extend(["VEP/CATALOG", "ClinGen", "GenCC"])
    if reason.startswith("K4_"):
        sources.extend(["VEP/CATALOG", "ClinGen Actionability"])
    if reason.startswith("K5_"):
        sources.extend(["AnnotSV", "ClinGen Dosage"])
    if reason.startswith("K6_"):
        sources.extend(["AnnotSV", "ClinGen", "GenCC"])
    if reason.startswith("K7_"):
        sources.append("PanelApp")
    if reason == "K8_ALPHAMISSENSE_SUPPORT":
        sources.append("AlphaMissense")
    if reason == "K8_CADD_SUPPORT_TARGET":
        sources.append("CADD")
    if first(row, "giab_flags"):
        sources.append("GIAB")

    for source in unique(sources):
        key = source.strip().lower()
        item = SOURCE_URLS.get(key)
        if item:
            add(*item)
    return links


def event_database_links(event: dict) -> list[tuple[str, str]]:
    all_links: list[tuple[str, str]] = []
    seen = set()
    for row in event["rows"]:
        for link in row_database_links(row, event["reason_code"]):
            if link not in seen:
                seen.add(link)
                all_links.append(link)
    return all_links


def relation_lines(rows: list[dict[str, str]], reason: str) -> list[str]:
    lines = []
    seen = set()
    research = reason.startswith(("K7_", "K8_"))
    for row in rows:
        if research:
            parts = []
            gene = first(row, "gene")
            if gene:
                parts.append(f"gene={gene}")
            if reason.startswith("K7_"):
                panels = first(row, "panelapp_panels")
                conf = first(row, "panelapp_confidence")
                moi = first(row, "panelapp_inheritance")
                if panels:
                    parts.append(f"PanelApp panels={panels}")
                if conf:
                    parts.append(f"confidence={conf}")
                if moi:
                    parts.append(f"inheritance={moi}")
            if reason == "K8_ALPHAMISSENSE_SUPPORT":
                score = first(row, "alpha_missense_pathogenicity")
                clazz = first(row, "alpha_missense_class")
                if score:
                    parts.append(f"AlphaMissense score={score}")
                if clazz:
                    parts.append(f"class={clazz}")
            if reason == "K8_CADD_SUPPORT_TARGET":
                raw = first(row, "cadd_raw")
                phred = first(row, "cadd_phred")
                if raw:
                    parts.append(f"CADD raw={raw}")
                if phred:
                    parts.append(f"PHRED={phred}")
            text = "; ".join(parts) or "research overlay"
        else:
            parts = []
            gene = first(row, "gene")
            disease = first(row, "disease_name")
            disease_id = first(row, "disease_id")
            classification = first(row, "gene_disease_classification")
            inheritance = first(row, "inheritance")
            if gene:
                parts.append(gene)
            if disease:
                parts.append(disease)
            if disease_id:
                parts.append(f"[{disease_id}]")
            if classification:
                parts.append(f"gene-disease={classification}")
            if inheritance:
                parts.append(f"MOI={inheritance}")
            text = " | ".join(parts)
            if not text:
                text = "variant-level finding"
        if text not in seen:
            seen.add(text)
            lines.append(text)
    return lines


def event_fields(event: dict, open_context: list[dict[str, str]] | None = None) -> list[tuple[str, str]]:
    rows = event["rows"]
    rep = rows[0]
    genes = unique(first(row, "gene") for row in rows)
    diseases = unique(first(row, "disease_name") for row in rows)
    fields = [
        ("Event", event["event_key"]),
        ("Тип", event.get("event_type", "") or "-"),
        ("Приоритет Stage 24", str(event.get("priority_rank", "?"))),
        ("Гены", ", ".join(genes) or "-"),
    ]
    if first(rep, "genotype"):
        fields.append(("Генотип / зиготность", f"{first(rep, 'genotype')} / {first(rep, 'zygosity') or '-'}"))
    consequence = first(rep, "consequence")
    impact = first(rep, "impact")
    if consequence or impact:
        fields.append(("VEP", " | ".join(x for x in [consequence, impact] if x)))
    af = first(rep, "max_af")
    if af:
        fields.append(("Max AF", af))
    cv = first(rep, "clinvar_significance")
    if cv:
        stars = first(rep, "clinvar_review_stars")
        fields.append(("ClinVar", cv + (f" | review stars={stars}" if stars else "")))
    conditions = first(rep, "clinvar_conditions")
    if conditions:
        fields.append(("ClinVar conditions", conditions))
    if diseases:
        fields.append(("Disease relations", "; ".join(diseases)))
    actionability = unique(first(row, "actionability") for row in rows)
    if actionability:
        fields.append(("Actionability", "; ".join(actionability)))
    dosage = unique(first(row, "dosage_evidence") for row in rows)
    if dosage:
        fields.append(("Dosage evidence", "; ".join(dosage)))
    tech = first(rep, "technical_confidence")
    if tech:
        fields.append(("Техническая уверенность", tech))
    giab = first(rep, "giab_flags")
    if giab:
        fields.append(("GIAB flags", giab))
    if open_context:
        open_rep = open_context[0]
        open_genes = unique(first(row, "gene") for row in open_context)
        open_reasons = unique(first(row, "reason_code") for row in open_context)
        open_diseases = unique(first(row, "disease_name") for row in open_context)
        open_gt = first(open_rep, "genotype")
        open_zyg = first(open_rep, "zygosity")
        if open_genes:
            fields.append(("Open context genes", ", ".join(open_genes)))
        if open_reasons:
            fields.append(("Open context reasons", ", ".join(open_reasons)))
        if open_gt or open_zyg:
            fields.append(("Open genotype / zygosity", f"{open_gt or '-'} / {open_zyg or '-'}"))
        open_consequence = first(open_rep, "consequence")
        if open_consequence:
            fields.append(("Open VEP consequence", open_consequence))
        open_af = first(open_rep, "max_af")
        if open_af:
            fields.append(("Open max AF", open_af))
        open_clinvar = first(open_rep, "clinvar_significance")
        if open_clinvar:
            fields.append(("Open ClinVar", open_clinvar))
        if open_diseases:
            fields.append(("Open disease relations", "; ".join(open_diseases)))
    fields.append(("Pipeline items", stage_links_for_event(event, "research" if event["reason_code"].startswith(("K7_", "K8_")) else "commercial")))
    ids = unique(first(row, "finding_id", "overlay_id") for row in rows)
    if ids:
        fields.append(("Stage-24 IDs", "; ".join(ids)))
    return fields


def reason_counts(grouped: dict[str, list[dict]], order: list[str]) -> list[dict]:
    rows = []
    for reason in order + sorted(set(grouped) - set(order)):
        events = grouped.get(reason, [])
        if not events:
            continue
        rows.append({
            "reason": reason,
            "events": len(events),
            "rows": sum(len(event["rows"]) for event in events),
        })
    return rows


def source_manifest_lines(usage: dict) -> list[str]:
    lines = []
    for item in usage.get("sources_used", []):
        parts = [str(item.get("source", "UNKNOWN")), str(item.get("publication_scope", "UNKNOWN"))]
        if item.get("license"):
            parts.append(f"license={item.get('license')}")
        if item.get("prohibited_uses"):
            parts.append("prohibited=" + ", ".join(map(str, item.get("prohibited_uses"))))
        lines.append(" | ".join(parts))
    return lines


def render_pdf(
    output: Path,
    title: str,
    patient: str,
    mode: str,
    legal: list[str],
    state: dict,
    usage: dict,
    stage_rows: list[dict],
    grouped: dict[str, list[dict]],
    category_meta: dict[str, tuple[str, str]],
    category_order: list[str],
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
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle, KeepTogether
    except ImportError as exc:
        raise RuntimeError(f"reportlab is required for PDF generation (recommended {REPORTLAB_VERSION})") from exc

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
    bold = next((p for p in bold_candidates if p.is_file()), regular)
    if regular is None:
        raise RuntimeError("No Cyrillic-capable local TTF font found (Arial/Segoe UI/DejaVu Sans)")
    pdfmetrics.registerFont(TTFont("DanteSans", str(regular)))
    pdfmetrics.registerFont(TTFont("DanteSansBold", str(bold)))

    doc = SimpleDocTemplate(
        str(output), pagesize=A4,
        rightMargin=14 * mm, leftMargin=14 * mm, topMargin=15 * mm, bottomMargin=15 * mm,
        title=title, author="Dante WGS pipeline",
    )
    styles = getSampleStyleSheet()
    for style in styles.byName.values():
        style.fontName = "DanteSans"
    styles.add(ParagraphStyle(name="DanteTitle", parent=styles["Title"], fontName="DanteSansBold", fontSize=18, leading=22, alignment=TA_CENTER, spaceAfter=10))
    styles.add(ParagraphStyle(name="DanteH1", parent=styles["Heading1"], fontName="DanteSansBold", fontSize=13, leading=16, spaceBefore=10, spaceAfter=5))
    styles.add(ParagraphStyle(name="DanteH2", parent=styles["Heading2"], fontName="DanteSansBold", fontSize=10, leading=12.5, spaceBefore=7, spaceAfter=3))
    styles.add(ParagraphStyle(name="DanteBody", parent=styles["BodyText"], fontName="DanteSans", fontSize=8.2, leading=10.7, spaceAfter=4))
    styles.add(ParagraphStyle(name="DanteSmall", parent=styles["BodyText"], fontName="DanteSans", fontSize=6.8, leading=8.5, spaceAfter=2))
    styles.add(ParagraphStyle(name="DanteTiny", parent=styles["BodyText"], fontName="DanteSans", fontSize=6.2, leading=7.5, spaceAfter=1.5))
    styles.add(ParagraphStyle(name="DanteWarn", parent=styles["BodyText"], fontName="DanteSans", fontSize=8.2, leading=10.5, backColor=colors.whitesmoke, borderPadding=6, spaceAfter=5))

    def p(text: str, style="DanteBody"):
        return Paragraph(text, styles[style])

    def esc(value: object) -> str:
        return html.escape(str(value if value is not None else ""))

    def links_markup(links: list[tuple[str, str]]) -> str:
        return " | ".join(f'<link href="{html.escape(url, quote=True)}">{esc(label)}</link>' for label, url in links)

    banner = "OPEN / COMMERCIAL-SAFE / DETERMINISTIC / NO AI" if mode == "commercial" else "RESEARCH ONLY / NON-COMMERCIAL / DETERMINISTIC / NO AI"
    all_events = {event["event_key"] for events in grouped.values() for event in events}
    counts = reason_counts(grouped, category_order)
    finding_rows = sum(item["rows"] for item in counts)

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

    story.append(p("Что именно содержит этот PDF", "DanteH1"))
    if mode == "commercial":
        story.append(p(
            "Это детерминированное представление полного open/commercial-safe слоя Stage 24. PDF не выбирает первые N вариантов и не выполняет дополнительную AI-интерпретацию: в него включены все уникальные события, имеющие open reason codes, с относящимися к ним Stage-24 gene/disease relations.",
            "DanteBody",
        ))
    else:
        story.append(p(
            "Это детерминированное представление полного research-only overlay слоя Stage 24. В него включены все уникальные события K7/K8 и все соответствующие research rows. Для ориентации показывается только локальный open context того же event_key; полный open-реестр остается в отдельном коммерческом PDF.",
            "DanteBody",
        ))

    scan = state.get("scan_stats", {})
    story.append(p("Сводка", "DanteH1"))
    summary_data = [
        [p("Показатель", "DanteSmall"), p("Значение", "DanteSmall")],
        [p("Уникальные события в этом отчете", "DanteSmall"), p(str(len(all_events)), "DanteSmall")],
        [p("Строки Stage-24 findings/overlays в этом отчете", "DanteSmall"), p(str(finding_rows), "DanteSmall")],
        [p("Stage 24 candidate events", "DanteSmall"), p(str(scan.get("candidate_events", "?")), "DanteSmall")],
        [p("Stage 24 open findings", "DanteSmall"), p(str(scan.get("open_findings", "?")), "DanteSmall")],
        [p("Stage 24 research overlays", "DanteSmall"), p(str(scan.get("research_overlays", "?")), "DanteSmall")],
        [p("AI/API при генерации", "DanteSmall"), p("НЕ ИСПОЛЬЗОВАЛИСЬ", "DanteSmall")],
    ]
    table = Table(summary_data, colWidths=[92 * mm, 82 * mm], repeatRows=1)
    table.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), "DanteSans"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]))
    story.append(table)

    story.append(p("Категории находок", "DanteH1"))
    count_data = [[p("Reason code", "DanteSmall"), p("Уникальные события", "DanteSmall"), p("Строки evidence", "DanteSmall")]]
    for item in counts:
        count_data.append([
            p(esc(item["reason"]), "DanteSmall"),
            p(str(item["events"]), "DanteSmall"),
            p(str(item["rows"]), "DanteSmall"),
        ])
    count_table = Table(count_data, colWidths=[105 * mm, 33 * mm, 36 * mm], repeatRows=1)
    count_table.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), "DanteSans"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke),
        ("GRID", (0, 0), (-1, -1), 0.2, colors.lightgrey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("FONTSIZE", (0, 0), (-1, -1), 6.8),
    ]))
    story.append(count_table)

    story.append(PageBreak())
    story.append(p("Детальный реестр Stage-24 находок", "DanteH1"))
    story.append(p(
        "Описание каждого раздела ниже является фиксированным правилом интерпретации категории. Конкретные значения внутри карточек получены из Stage-24 TSV и подключенных баз; рядом приведены ссылки на соответствующие базы данных.",
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
            label = event_label(event)
            heading = f"{global_index}. {label} | {reason}"
            fields = event_fields(event, (open_by_event or {}).get(event["event_key"]) if mode == "research" else None)
            block = [p(esc(heading), "DanteH2")]
            for key, value in fields:
                if value:
                    block.append(p(f"<b>{esc(key)}:</b> {esc(value)}", "DanteSmall"))
            relations = relation_lines(event["rows"], reason)
            if relations:
                block.append(p(f"<b>{'Research rows' if mode == 'research' else 'Gene/disease relations'} ({len(relations)}):</b>", "DanteSmall"))
                for line in relations:
                    block.append(p("- " + esc(line), "DanteTiny"))
            links = event_database_links(event)
            if mode == "research":
                seen_links = set(links)
                for open_row in (open_by_event or {}).get(event["event_key"], []):
                    for link in row_database_links(open_row, first(open_row, "reason_code")):
                        if link not in seen_links:
                            seen_links.add(link)
                            links.append(link)
            if links:
                block.append(p("<b>Базы:</b> " + links_markup(links), "DanteSmall"))
            sources = unique(source for row in event["rows"] for source in split_multi(first(row, "evidence_sources")))
            if sources:
                block.append(p("<b>Evidence sources из Stage 24:</b> " + esc("; ".join(sources)), "DanteTiny"))
            story.append(KeepTogether(block[:3]))
            story.extend(block[3:])
            story.append(Spacer(1, 1.2 * mm))

    story.append(PageBreak())
    story.append(p("Knowledge sources и ограничения лицензий", "DanteH1"))
    for line in source_manifest_lines(usage):
        story.append(p(esc(line), "DanteSmall"))
    source_links = []
    for item in usage.get("sources_used", []):
        key = str(item.get("source", "")).strip().lower()
        if key in SOURCE_URLS:
            source_links.append(SOURCE_URLS[key])
    if source_links:
        dedup_links = list(dict.fromkeys(source_links))
        story.append(p("Ссылки: " + links_markup(dedup_links), "DanteSmall"))

    story.append(p("Воспроизводимость: pipeline items 1-24", "DanteH1"))
    data = [[p("#", "DanteTiny"), p("Этап", "DanteTiny"), p("Артефакт / provenance", "DanteTiny"), p("Статус", "DanteTiny")]]
    for row in stage_rows:
        data.append([
            p(str(row["item"]), "DanteTiny"),
            p(esc(row["name"]), "DanteTiny"),
            p(esc(row["artifact"]), "DanteTiny"),
            p(esc(row["status"]), "DanteTiny"),
        ])
    stages = Table(data, colWidths=[8 * mm, 49 * mm, 88 * mm, 29 * mm], repeatRows=1, splitByRow=1)
    stages.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, -1), "DanteSans"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke),
        ("GRID", (0, 0), (-1, -1), 0.2, colors.lightgrey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 2),
        ("RIGHTPADDING", (0, 0), (-1, -1), 2),
        ("TOPPADDING", (0, 0), (-1, -1), 1.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5),
    ]))
    story.append(stages)

    story.append(p("Ограничения текущего аналитического покрытия", "DanteH1"))
    limitations = [
        "PRS/GWAS не входят в текущий детерминированный Stage-24 clinical/open contract.",
        "Полные CYP/PGx haplotypes, HLA typing, repeat expansions, global phasing и mtDNA heteroplasmy требуют отдельных валидированных модулей.",
        "GIAB difficult-region overlap снижает техническую уверенность, но не удаляет событие автоматически.",
        "Для структурных вариантов биологическая значимость и техническая достоверность вызова оцениваются раздельно.",
        "Категории K2/K3/K4 характеризуют gene-level evidence context и не являются автоматической variant-level pathogenicity classification.",
    ]
    for item in limitations:
        story.append(p(esc(item), "DanteBody"))

    footer_text = "OPEN / COMMERCIAL-SAFE / NO AI" if mode == "commercial" else "RESEARCH ONLY / NON-COMMERCIAL / NO AI"

    def footer(canvas, doc_):
        canvas.saveState()
        canvas.setFont("DanteSans", 6.2)
        canvas.drawString(14 * mm, 8 * mm, footer_text)
        canvas.drawRightString(A4[0] - 14 * mm, 8 * mm, f"стр. {doc_.page}")
        canvas.restoreState()

    output.parent.mkdir(parents=True, exist_ok=True)
    doc.build(story, onFirstPage=footer, onLaterPages=footer)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate deterministic split WGS PDF reports from Stage 24")
    parser.add_argument("mode", choices=("commercial", "research", "self-test"))
    parser.add_argument("--root", type=Path)
    parser.add_argument("--patient-root", type=Path)
    parser.add_argument("--accept-research-restrictions", action="store_true")
    return parser.parse_args()


def self_test() -> int:
    open_usage = {
        "schema": USAGE_SCHEMA,
        "scope": "OPEN_FINDINGS_ONLY",
        "sources_used": [{"source": "ClinVar", "publication_scope": "COMMERCIAL_ALLOWED", "prohibited_uses": []}],
    }
    validate_usage("commercial", open_usage)
    bad = {
        "schema": USAGE_SCHEMA,
        "scope": "OPEN_PLUS_RESEARCH_OVERLAYS",
        "sources_used": [{"source": "PanelApp", "publication_scope": "RESEARCH_ONLY", "prohibited_uses": ["commercial"]}],
    }
    try:
        validate_usage("commercial", bad)
    except RuntimeError:
        pass
    else:
        raise AssertionError("commercial mode accepted research-only source")
    validate_usage("research", bad)
    rows = [
        {"finding_id": "OF1", "priority_rank": "2", "event_type": "SNV_INDEL", "event_key": "GRCh37:1:1:A:G", "gene": "GENE1", "reason_code": "K1_EXACT_CLINVAR_PLP"},
        {"finding_id": "OF2", "priority_rank": "2", "event_type": "SNV_INDEL", "event_key": "GRCh37:1:1:A:G", "gene": "GENE1", "disease_name": "Disease B", "reason_code": "K1_EXACT_CLINVAR_PLP"},
        {"finding_id": "OF3", "priority_rank": "4", "event_type": "SNV_INDEL", "event_key": "GRCh37:2:2:C:T", "gene": "GENE2", "reason_code": "K3_RARE_MISSENSE_VALIDATED_GENE"},
    ]
    grouped = group_reason_events(rows)
    assert len(grouped["K1_EXACT_CLINVAR_PLP"]) == 1
    assert len(grouped["K1_EXACT_CLINVAR_PLP"][0]["rows"]) == 2
    assert len(grouped["K3_RARE_MISSENSE_VALIDATED_GENE"]) == 1
    print("SELF-TEST PASS | deterministic split PDFs + all Stage-24 rows + no AI/API")
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
    state = read_json(state_path)
    if state.get("status") != "complete" or state.get("schema") != STAGE24_SCHEMA or state.get("rules_schema") != RULES_SCHEMA:
        raise RuntimeError("Stage 24 is not READY under rules-v8")

    out_root = patient_root / "Dante_results" / "reports"
    out_root.mkdir(parents=True, exist_ok=True)
    stage_rows = stage_presence(root, patient_root)
    open_rows = read_tsv(open_path)
    open_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in open_rows:
        open_by_event[first(row, "event_key") or first(row, "finding_id")].append(row)

    if args.mode == "commercial":
        usage = read_json(open_usage_path)
        validate_usage("commercial", usage)
        grouped = group_reason_events(open_rows)
        output = out_root / "Dante_genome_report_COMMERCIAL.pdf"
        render_pdf(
            output,
            "Dante WGS - итоговый геномный отчет",
            patient,
            "commercial",
            COMMERCIAL_LEGAL,
            state,
            usage,
            stage_rows,
            grouped,
            OPEN_CATEGORY_META,
            OPEN_CATEGORY_ORDER,
        )
        counts = reason_counts(grouped, OPEN_CATEGORY_ORDER)
        manifest = {
            "schema": REPORT_SCHEMA,
            "mode": "COMMERCIAL",
            "generation_mode": "DETERMINISTIC_STAGE24_ONLY",
            "ai_used": False,
            "network_used_by_report_generator": False,
            "all_stage24_findings_included": True,
            "patient": patient,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "pipeline_release": state.get("pipeline_release"),
            "rules_schema": state.get("rules_schema"),
            "pipeline_commit": pipeline_commit(root),
            "source_files": {
                "state": {"name": state_path.name, "sha256": sha256(state_path)},
                "open_findings": {"name": open_path.name, "sha256": sha256(open_path)},
                "usage": {"name": open_usage_path.name, "sha256": sha256(open_usage_path)},
            },
            "research_layer_opened": False,
            "unique_events": len({e["event_key"] for events in grouped.values() for e in events}),
            "finding_rows": len(open_rows),
            "category_counts": counts,
            "output": {"name": output.name, "sha256": sha256(output)},
        }
        (out_root / "Dante_genome_report_COMMERCIAL.manifest.json").write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        print(f"Created COMMERCIAL deterministic report: {output}")
        print(f"Included all {len(open_rows)} Stage-24 open finding rows; research-only TSV was not opened; OpenAI/API was not used.")
        return 0

    if not args.accept_research_restrictions:
        if not sys.stdin.isatty():
            raise RuntimeError("Research mode requires --accept-research-restrictions")
        print("RESEARCH ONLY: non-commercial; not for diagnosis, medical decisions or healthcare-service use.")
        if input("Type RESEARCH to continue: ").strip() != "RESEARCH":
            raise RuntimeError("Research report cancelled")
    if not research_path.is_file() or not research_usage_path.is_file():
        raise FileNotFoundError("Research Stage-24 outputs are missing")

    usage = read_json(research_usage_path)
    validate_usage("research", usage)
    overlays = read_tsv(research_path)
    grouped = group_reason_events(overlays)
    output = out_root / "Dante_genome_report_RESEARCH_ONLY.pdf"
    render_pdf(
        output,
        "Dante WGS - исследовательский геномный отчет",
        patient,
        "research",
        RESEARCH_LEGAL,
        state,
        usage,
        stage_rows,
        grouped,
        RESEARCH_CATEGORY_META,
        RESEARCH_CATEGORY_ORDER,
        open_by_event=open_by_event,
    )
    counts = reason_counts(grouped, RESEARCH_CATEGORY_ORDER)
    manifest = {
        "schema": REPORT_SCHEMA,
        "mode": "RESEARCH_ONLY",
        "generation_mode": "DETERMINISTIC_STAGE24_ONLY",
        "ai_used": False,
        "network_used_by_report_generator": False,
        "all_stage24_findings_included": True,
        "patient": patient,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "pipeline_release": state.get("pipeline_release"),
        "rules_schema": state.get("rules_schema"),
        "pipeline_commit": pipeline_commit(root),
        "source_files": {
            "state": {"name": state_path.name, "sha256": sha256(state_path)},
            "open_findings": {"name": open_path.name, "sha256": sha256(open_path)},
            "research_overlays": {"name": research_path.name, "sha256": sha256(research_path)},
            "usage": {"name": research_usage_path.name, "sha256": sha256(research_usage_path)},
        },
        "research_layer_opened": True,
        "unique_events": len({e["event_key"] for events in grouped.values() for e in events}),
        "finding_rows": len(overlays),
        "category_counts": counts,
        "output": {"name": output.name, "sha256": sha256(output)},
    }
    (out_root / "Dante_genome_report_RESEARCH_ONLY.manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Created RESEARCH ONLY deterministic report: {output}")
    print(f"Included all {len(overlays)} Stage-24 research overlay rows; OpenAI/API was not used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
