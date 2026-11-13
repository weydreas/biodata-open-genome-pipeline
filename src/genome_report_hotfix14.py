"""Genotype-table TXT publication scope over report hotfix13.

The previous WGS and GENOTYPING_ARRAY_VCF publication paths are delegated
unchanged.  This layer adds a third independent report contract for measured
23andMe-style GRCh37 genotype tables.
"""
from __future__ import annotations

import genome_report_hotfix13 as h13
import genome_report_base_v2 as base

PUBLICATION_SCOPE_HOTFIX = "dante-report-genotype-txt-scope-v14"
TXT_KIND = "GENOTYPE_TABLE_TXT"
TXT_STAGE_ITEMS = {1, 2, 6, 7, 8, 9, 10, 11, 12, 21, 22, 23, 24}
_ORIGINAL_RENDERER = h13.render_report_pdf

TXT_COMMERCIAL_LEGAL = [
    "Информационно-аналитический отчет по данным генотипирования в 23andMe-style TXT (GRCh37/build 37). Это не отчет по полногеномному секвенированию (WGS) и не интерпретация полного генома. Отчет не является медицинской услугой, медицинским заключением, диагнозом, назначением лечения или заменой очной консультации врача.",
    "Область отрицательных выводов ограничена локусами, реально присутствующими в исходной таблице. Отсутствующий локус и строка с no-call (`--`) означают UNKNOWN и не могут интерпретироваться как референсный генотип или отсутствие генетического риска.",
    "Для sequence-level аннотации используются только A/C/G/T генотипы, разрешенные относительно установленного GRCh37 FASTA. Символические I/D, II/DD/DI варианты не содержат последовательности вставки/делеции, поэтому сохраняются в прямом genotype export/QC, но не включаются в VEP/Stage 24 как последовательностно определенные indel-варианты.",
    "Исходная genotype-table не содержит read-level метрик DP/GQ/QUAL/FILTER. Значимые находки следует подтверждать валидированным клинико-диагностическим методом до медицинских решений.",
    "Отчет автоматически формируется из детерминированных результатов Stage 24. При его создании не используются генеративный ИИ, внешняя интерпретация или передача геномных данных в OpenAI/API.",
    "Классификация вариантов и связи ген-заболевание отражают состояние подключенных баз знаний на дату их установки и могут изменяться при обновлении ClinVar, ClinGen, GenCC и других источников.",
]

TXT_RESEARCH_LEGAL = [
    "RESEARCH ONLY / ТОЛЬКО ДЛЯ ИССЛЕДОВАТЕЛЬСКОГО ИСПОЛЬЗОВАНИЯ.",
    "Исходные данные — 23andMe-style genotype TXT на GRCh37/build 37, а не WGS. Отсутствующие локусы и no-call строки остаются UNKNOWN.",
    "Sequence-level аннотация ограничена A/C/G/T генотипами, однозначно разрешенными относительно GRCh37 FASTA; символические I/D indel-коды не интерпретируются как известные последовательности.",
    "Этот отчет не предназначен для коммерческого использования, диагностики, медицинского решения или оказания медицинской услуги.",
    "Отчет автоматически формируется из Stage-24 research-only overlays. При его создании не используются генеративный ИИ, внешняя интерпретация или передача геномных данных в OpenAI/API.",
    "Исследовательские предикторы и panel membership не являются клинической классификацией варианта и используются только как дополнительный слой гипотез.",
]


def is_txt_state(state: dict) -> bool:
    profile = state.get("input_profile") or {}
    return str(profile.get("kind") or "") == TXT_KIND


def render_report_pdf(output, **kwargs):
    state = kwargs.get("state") or {}
    if not is_txt_state(state):
        return _ORIGINAL_RENDERER(output, **kwargs)

    old_commercial = base.COMMERCIAL_LEGAL
    old_research = base.RESEARCH_LEGAL
    base.COMMERCIAL_LEGAL = TXT_COMMERCIAL_LEGAL
    base.RESEARCH_LEGAL = TXT_RESEARCH_LEGAL
    try:
        mode = str(kwargs.get("mode") or "commercial")
        appendix = bool(kwargs.get("is_reference_appendix"))
        if mode == "research":
            kwargs["title"] = (
                "Dante Genotype Table - RESEARCH ONLY: common major reference appendix"
                if appendix else
                "Dante Genotype Table - исследовательский отчет по измеренным локусам"
            )
        else:
            kwargs["title"] = (
                "Dante Genotype Table - справочное приложение: common major alleles"
                if appendix else
                "Dante Genotype Table - отчет по измеренным локусам"
            )
        kwargs["stage_rows"] = [
            row for row in list(kwargs.get("stage_rows") or [])
            if int(row.get("item", -1)) in TXT_STAGE_ITEMS
        ]
        return _ORIGINAL_RENDERER(output, **kwargs)
    finally:
        base.COMMERCIAL_LEGAL = old_commercial
        base.RESEARCH_LEGAL = old_research


def main() -> int:
    old = h13.render_report_pdf
    h13.render_report_pdf = render_report_pdf
    try:
        return h13.main()
    finally:
        h13.render_report_pdf = old


for _name in dir(h13):
    if _name.startswith("__") or _name in {"main", "render_report_pdf"}:
        continue
    globals().setdefault(_name, getattr(h13, _name))


if __name__ == "__main__":
    raise SystemExit(main())
