"""Input-scope aware publication layer over report hotfix12.

For genotyping-array patients the same Stage-24 findings are rendered, but the
report must never describe the source as WGS or imply negative evidence outside
interrogated loci. WGS patients retain the exact validated hotfix12 behavior.
"""
from __future__ import annotations

import genome_report_hotfix12 as h12
import genome_report_base_v2 as base

PUBLICATION_SCOPE_HOTFIX = "dante-report-input-scope-v13"
ARRAY_KIND = "GENOTYPING_ARRAY_VCF"
ARRAY_STAGE_ITEMS = {1, 2, 6, 7, 8, 9, 10, 11, 12, 21, 22, 23, 24}
_ORIGINAL_RENDERER = h12.render_report_pdf

ARRAY_COMMERCIAL_LEGAL = [
    "Информационно-аналитический отчет по данным генотипирования на SNP/variant array (VCF, hg19/GRCh37). Это не отчет по полногеномному секвенированию (WGS). Отчет не является медицинской услугой, медицинским заключением, диагнозом, назначением лечения или заменой очной консультации врача.",
    "Отчет автоматически формируется из детерминированных результатов Stage 24. При его создании не используются генеративный ИИ, внешняя интерпретация или передача геномных данных в OpenAI/API.",
    "Область отрицательных выводов ограничена локусами, реально присутствующими в исходном array VCF. Отсутствие варианта или локуса вне измеренной панели означает UNKNOWN и не может интерпретироваться как референсный генотип или отсутствие генетического риска.",
    "Исходный array VCF может не содержать read-level метрики DP/GQ/QUAL/FILTER. Значимые находки следует подтверждать валидированным клинико-диагностическим методом до медицинских решений.",
    "Классификация вариантов и связи ген-заболевание отражают состояние подключенных баз знаний на дату их установки и могут изменяться при обновлении ClinVar, ClinGen, GenCC и других источников.",
]

ARRAY_RESEARCH_LEGAL = [
    "RESEARCH ONLY / ТОЛЬКО ДЛЯ ИССЛЕДОВАТЕЛЬСКОГО ИСПОЛЬЗОВАНИЯ.",
    "Исходные данные — генотипирование на SNP/variant array (VCF, hg19/GRCh37), а не WGS. Покрытие ограничено локусами исходной панели; отсутствие данных за пределами панели остается UNKNOWN.",
    "Этот отчет не предназначен для коммерческого использования, диагностики, медицинского решения или оказания медицинской услуги.",
    "Отчет автоматически формируется из Stage-24 research-only overlays. При его создании не используются генеративный ИИ, внешняя интерпретация или передача геномных данных в OpenAI/API.",
    "Исследовательские предикторы и panel membership не являются клинической классификацией варианта и используются только как дополнительный слой гипотез.",
]


def is_array_state(state: dict) -> bool:
    profile = state.get("input_profile") or {}
    return str(profile.get("kind") or "") == ARRAY_KIND


def render_report_pdf(output, **kwargs):
    state = kwargs.get("state") or {}
    if not is_array_state(state):
        return _ORIGINAL_RENDERER(output, **kwargs)

    old_commercial = base.COMMERCIAL_LEGAL
    old_research = base.RESEARCH_LEGAL
    base.COMMERCIAL_LEGAL = ARRAY_COMMERCIAL_LEGAL
    base.RESEARCH_LEGAL = ARRAY_RESEARCH_LEGAL
    try:
        mode = str(kwargs.get("mode") or "commercial")
        appendix = bool(kwargs.get("is_reference_appendix"))
        if mode == "research":
            kwargs["title"] = (
                "Dante Genotyping Array - RESEARCH ONLY: common major reference appendix"
                if appendix else
                "Dante Genotyping Array - исследовательский отчет по измеренным локусам"
            )
        else:
            kwargs["title"] = (
                "Dante Genotyping Array - справочное приложение: common major alleles"
                if appendix else
                "Dante Genotyping Array - отчет по измеренным локусам"
            )
        kwargs["stage_rows"] = [
            row for row in list(kwargs.get("stage_rows") or [])
            if int(row.get("item", -1)) in ARRAY_STAGE_ITEMS
        ]
        return _ORIGINAL_RENDERER(output, **kwargs)
    finally:
        base.COMMERCIAL_LEGAL = old_commercial
        base.RESEARCH_LEGAL = old_research


def main() -> int:
    old = h12.render_report_pdf
    h12.render_report_pdf = render_report_pdf
    try:
        return h12.main()
    finally:
        h12.render_report_pdf = old


for _name in dir(h12):
    if _name.startswith("__") or _name in {"main", "render_report_pdf"}:
        continue
    globals().setdefault(_name, getattr(h12, _name))


if __name__ == "__main__":
    raise SystemExit(main())
