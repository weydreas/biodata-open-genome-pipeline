"""Clinical-risk publication layer over report hotfix14.

Stage 24 remains scientifically unchanged. This layer consumes the deterministic
clinical-risk evidence artifact and augments the disease-chapter summary tables
with patient allele/genotype, normalized ClinVar classification and a conservative
quantitative risk/penetrance field. PDF generation itself performs no network or
AI calls.
"""
from __future__ import annotations

from contextlib import contextmanager
import html
import json
from pathlib import Path

import clinical_risk_evidence as risk
import genome_report_hotfix14 as h14
import genome_report_hotfix9 as h9

h13 = h14.h13
h12 = h13.h12
h11 = h12.h11
base = h11.base

PUBLICATION_CONTRACT = "dante-report-clinical-risk-v15"
GENERATION_MODE = "DETERMINISTIC_STAGE24_PLUS_LOCAL_RISK_EVIDENCE"
_ORIGINAL_RENDERER = h14.render_report_pdf
_ACTIVE_EVENTS: dict[str, dict] = {}
_BY_LABEL: dict[str, list[dict]] = {}


def _plain(cell) -> str:
    if hasattr(cell, "getPlainText"):
        try:
            return str(cell.getPlainText())
        except Exception:
            pass
    return str(getattr(cell, "text", cell))


def _register_events(events: dict[str, dict]) -> None:
    global _ACTIVE_EVENTS, _BY_LABEL
    _ACTIVE_EVENTS = dict(events)
    by_label: dict[str, list[dict]] = {}
    for event_key, event in _ACTIVE_EVENTS.items():
        labels = []
        rsid = str(event.get("rsid") or "").strip()
        if rsid:
            labels.append(rsid)
        if event_key.startswith("GRCh37:"):
            labels.append(event_key[len("GRCh37:"):])
        labels.append(event_key)
        for label in labels:
            by_label.setdefault(label, []).append(event)
    _BY_LABEL = by_label


def _risk_event_for_label(label: str, genes_text: str = "") -> dict | None:
    candidates = _BY_LABEL.get((label or "").strip(), [])
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    genes = {g.strip().upper() for g in genes_text.replace(";", ",").split(",") if g.strip()}
    for candidate in candidates:
        if genes & {str(g).upper() for g in candidate.get("genes", [])}:
            return candidate
    return candidates[0]


def _risk_event_for_event(event: dict) -> dict | None:
    key = str(event.get("event_key") or "")
    if key and key in _ACTIVE_EVENTS:
        return _ACTIVE_EVENTS[key]
    return _risk_event_for_label(base.event_label(event))


def _allele_cell(event: dict | None) -> str:
    if not event:
        return "-"
    alt = str(event.get("patient_alt_allele") or "").strip()
    genotype = str(event.get("patient_genotype_display") or "").strip()
    if alt and genotype and genotype != "-":
        return f"{alt}; {genotype}"
    return alt or genotype or "-"


def _risk_cell(event: dict | None) -> str:
    if not event:
        return "-"
    metrics = list(event.get("metrics") or [])
    if metrics:
        return "; ".join(str(item.get("risk_summary_ru") or "") for item in metrics[:2] if item.get("risk_summary_ru")) or "-"
    return "Нет надёжной количественной оценки"


def _classification_cell(event: dict | None) -> str:
    if not event:
        return "-"
    return str(event.get("clinvar_classification") or "Not classified")


def _is_variant_index(data) -> bool:
    if not isinstance(data, list) or not data or not isinstance(data[0], list) or len(data[0]) != 3:
        return False
    return [_plain(cell) for cell in data[0]] == ["rsID / вариант", "Ген(ы)", "Связанные состояния"]


@contextmanager
def enriched_summary_tables():
    """Intercept only the disease-chapter short variant index Table constructor."""
    try:
        from reportlab import platypus
        from reportlab.lib.units import mm
    except ModuleNotFoundError:
        yield
        return

    original_table = platypus.Table
    paragraph = platypus.Paragraph

    def make_para(text: str, template):
        style = getattr(template, "style", None)
        if style is None:
            return text
        return paragraph(html.escape(str(text)), style)

    def table_with_clinical_risk(*args, **kwargs):
        data = args[0] if args else kwargs.get("data")
        if _is_variant_index(data):
            header_template = data[0][0]
            enriched = [[
                make_para("rsID / вариант", header_template),
                make_para("Аллель пациента", header_template),
                make_para("Ген(ы)", header_template),
                make_para("Классификация ClinVar", header_template),
                make_para("Риск / пенетрантность", header_template),
                make_para("Связанные состояния", header_template),
            ]]
            for row in data[1:]:
                label = _plain(row[0])
                genes = _plain(row[1])
                evidence = _risk_event_for_label(label, genes)
                enriched.append([
                    row[0],
                    make_para(_allele_cell(evidence), row[0]),
                    row[1],
                    make_para(_classification_cell(evidence), row[0]),
                    make_para(_risk_cell(evidence), row[0]),
                    row[2],
                ])
            if args:
                args = (enriched, *args[1:])
            else:
                kwargs["data"] = enriched
            kwargs["colWidths"] = [22 * mm, 24 * mm, 19 * mm, 27 * mm, 40 * mm, 42 * mm]
        return original_table(*args, **kwargs)

    platypus.Table = table_with_clinical_risk
    try:
        yield
    finally:
        platypus.Table = original_table


def render_report_pdf(output, **kwargs):
    with enriched_summary_tables():
        return _ORIGINAL_RENDERER(output, **kwargs)


def _event_fields_with_risk(original):
    def wrapped(event: dict, open_context=None):
        fields = list(original(event, open_context))
        evidence = _risk_event_for_event(event)
        if not evidence:
            return fields
        additions = [
            ("Аллель пациента", _allele_cell(evidence)),
            ("Классификация ClinVar (норм.)", _classification_cell(evidence)),
            ("Риск / пенетрантность", _risk_cell(evidence)),
        ]
        insert_at = next((i + 1 for i, (key, _) in enumerate(fields) if key == "Генотип / зиготность"), 4)
        for item in reversed(additions):
            fields.insert(insert_at, item)
        return fields
    return wrapped


def _manifest_with_risk(original):
    def wrapped(**kwargs):
        output_path = Path(kwargs["output_path"])
        risk_dir = output_path.parent.parent / "clinical_risk_evidence"
        evidence_json = risk_dir / "Dante_clinical_risk_evidence.json"
        evidence_tsv = risk_dir / "Dante_clinical_risk_evidence.tsv"
        evidence_state = risk_dir / "Dante_clinical_risk_evidence.state.json"
        for required in (evidence_json, evidence_tsv, evidence_state):
            if not required.is_file():
                raise FileNotFoundError(required)
        source_files = dict(kwargs.get("source_files") or {})
        source_files["clinical_risk_evidence"] = {"name": evidence_json.name, "sha256": risk.sha256_file(evidence_json)}
        source_files["clinical_risk_evidence_tsv"] = {"name": evidence_tsv.name, "sha256": risk.sha256_file(evidence_tsv)}
        source_files["clinical_risk_state"] = {"name": evidence_state.name, "sha256": risk.sha256_file(evidence_state)}
        kwargs["source_files"] = source_files
        original(**kwargs)
        manifest = json.loads(output_path.read_text(encoding="utf-8"))
        state = json.loads(evidence_state.read_text(encoding="utf-8"))
        manifest["publication_contract"] = PUBLICATION_CONTRACT
        manifest["generation_mode"] = GENERATION_MODE
        manifest["clinical_risk_evidence_schema"] = risk.SCHEMA
        manifest["clinical_risk_source_policy"] = state.get("source_policy", risk.SOURCE_POLICY)
        manifest["clinical_risk_network_used"] = bool(state.get("network_used", False))
        manifest["clinical_risk_patient_identifiers_sent_external"] = bool(state.get("patient_identifiers_sent_external", False))
        manifest["clinical_risk_counts"] = state.get("counts", {})
        manifest["clinical_risk_interpretation_boundary"] = (
            "Quantitative risk is shown only when supported by explicit numeric local ClinGen evidence. "
            "Gene-condition estimates are labelled as such and are not silently presented as variant-specific penetrance."
        )
        output_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return wrapped


def _load_or_build(root: Path, patient_root: Path) -> dict[str, dict]:
    ready, _ = risk.is_ready(root, patient_root)
    if not ready:
        risk.build(root, patient_root)
    out = risk.output_paths(patient_root)["json"]
    payload = json.loads(out.read_text(encoding="utf-8"))
    if payload.get("schema") != risk.SCHEMA or payload.get("status") != "complete":
        raise RuntimeError("Clinical-risk evidence artifact is not READY")
    events = payload.get("events")
    if not isinstance(events, dict):
        raise RuntimeError("Clinical-risk evidence events map is invalid")
    return events


def self_test() -> int:
    risk.self_test()
    _register_events({
        "GRCh37:7:140453136:A:T": {
            "event_key": "GRCh37:7:140453136:A:T",
            "rsid": "rs123",
            "genes": ["BRAF"],
            "patient_alt_allele": "T",
            "patient_genotype_display": "A/T",
            "clinvar_classification": "Pathogenic",
            "metrics": [{"risk_summary_ru": "Пенетрантность 80% (ClinGen Actionability, Adult; gene-condition)"}],
        }
    })
    event = _risk_event_for_label("rs123", "BRAF")
    assert event and _allele_cell(event) == "T; A/T"
    assert _classification_cell(event) == "Pathogenic"
    assert "80%" in _risk_cell(event)
    assert _risk_event_for_label("7:140453136:A:T", "BRAF") is event
    print("SELF-TEST PASS | clinical-risk report summary augmentation + provenance contract")
    return 0


def main() -> int:
    args = h9.parse_args()
    if args.mode == "self-test":
        return self_test()
    if not args.root or not args.patient_root:
        raise RuntimeError("--root and --patient-root are required")
    root = args.root.resolve()
    patient_root = args.patient_root.resolve()
    _register_events(_load_or_build(root, patient_root))

    original_renderer = h14.render_report_pdf
    original_fields = base.event_fields
    original_manifest = h11._write_manifest_v11
    h14.render_report_pdf = render_report_pdf
    base.event_fields = _event_fields_with_risk(original_fields)
    h11._write_manifest_v11 = _manifest_with_risk(original_manifest)
    try:
        return h14.main()
    finally:
        h14.render_report_pdf = original_renderer
        base.event_fields = original_fields
        h11._write_manifest_v11 = original_manifest


for _name in dir(h14):
    if _name.startswith("__") or _name in {"main", "render_report_pdf"}:
        continue
    globals().setdefault(_name, getattr(h14, _name))


if __name__ == "__main__":
    raise SystemExit(main())
