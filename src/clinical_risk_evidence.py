#!/usr/bin/env python3
"""Deterministic clinical-risk evidence layer for Dante reports.

This stage sits between Stage 24 and PDF publication. It does not change Stage-24
finding rules or pathogenicity. It materializes report-facing fields from the
already-local patient findings and the locally installed ClinGen Actionability
snapshots:

- patient ALT allele and readable genotype;
- normalized ClinVar classification;
- conservative quantitative risk/penetrance evidence when an explicit number is
  present in a matching ClinGen gene-condition record;
- explicit unavailable state when a quantitative estimate is not supported.

No network requests are made. No patient rsIDs, genes or genotypes leave the local
workspace. Bare ClinGen semi-quantitative scores are never converted into disease
risk. A quantitative metric is accepted only when its source field contains an
explicit percentage/fraction or explicitly-labelled RR/OR/HR value.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Iterator

SCHEMA = "dante-clinical-risk-evidence-v1"
STATE_SCHEMA = "dante-clinical-risk-evidence-state-v1"
STAGE24_SCHEMA = "dante-knowledge-enrichment-v1"
RULES_SCHEMA = "dante-knowledge-rules-v8"
SOURCE_POLICY = "LOCAL_OPEN_SNAPSHOT_ONLY"
MAX_METRICS_PER_EVENT = 3

TSV_COLUMNS = [
    "event_key", "rsid", "gene", "ref", "alt", "genotype", "patient_alt_allele",
    "patient_genotype_display", "zygosity", "clinvar_classification",
    "clinvar_significance_raw", "clinvar_review_stars", "conditions",
    "risk_metric_type", "risk_value", "risk_unit", "risk_age", "risk_condition",
    "risk_scope", "risk_source", "risk_context", "source_doc_id", "source_field",
    "risk_summary_ru",
]

RISK_KEYWORDS = ("penetrance", "risk", "likelihood", "incidence", "prevalence", "probability", "frequency")
BLOCKED_CONTEXT = ("effectiveness", "intervention", "severity", "natureofintervention", "overallscore", "overall_score")
METRIC_RANK = {
    "penetrance": 80,
    "lifetime_risk": 75,
    "cumulative_risk": 70,
    "absolute_risk": 65,
    "probability": 60,
    "incidence": 55,
    "prevalence": 50,
    "relative_risk": 40,
    "hazard_ratio": 35,
    "odds_ratio": 30,
    "frequency": 20,
}
STOPWORDS = {
    "and", "or", "of", "the", "with", "without", "related", "syndrome", "disease",
    "disorder", "type", "types", "familial", "hereditary", "autosomal", "dominant",
    "recessive", "adult", "pediatric", "clinical", "risk", "cancer",
}


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    return {"path": str(path), "length": path.stat().st_size, "sha256": sha256_file(path)}


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def write_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(part, path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def write_tsv_atomic(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    with part.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in TSV_COLUMNS})
    os.replace(part, path)


def split_multi(value: str) -> list[str]:
    value = (value or "").strip().strip('"')
    if not value or value in {".", "-"}:
        return []
    return [v.strip().strip('"') for v in re.split(r"\s*[|;,]\s*", value) if v.strip().strip('"')]


def unique(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        v = (value or "").strip()
        if not v or v in seen:
            continue
        seen.add(v)
        result.append(v)
    return result


def normalize_clinvar_classification(values: Iterable[str]) -> tuple[str, str]:
    raw_values = unique(v for value in values for v in split_multi(value.replace("/", "|")))
    normalized = {re.sub(r"[_\s]+", " ", v).strip().lower() for v in raw_values}
    has_p = "pathogenic" in normalized
    has_lp = "likely pathogenic" in normalized
    has_b = "benign" in normalized
    has_lb = "likely benign" in normalized
    if any("conflicting" in x for x in normalized) or ((has_p or has_lp) and (has_b or has_lb)):
        label = "Conflicting interpretations"
    elif has_p and has_lp:
        label = "Pathogenic / Likely pathogenic"
    elif has_p:
        label = "Pathogenic"
    elif has_lp:
        label = "Likely pathogenic"
    elif has_b and has_lb:
        label = "Benign / Likely benign"
    elif has_b:
        label = "Benign"
    elif has_lb:
        label = "Likely benign"
    elif any("uncertain" in x or x == "vus" for x in normalized):
        label = "Uncertain significance"
    elif any("risk factor" in x for x in normalized):
        label = "Risk factor"
    elif any("protective" in x for x in normalized):
        label = "Protective"
    elif any("drug response" in x for x in normalized):
        label = "Drug response"
    elif raw_values:
        label = " / ".join(raw_values)
    else:
        label = "Not classified"
    return label, " | ".join(raw_values)


def readable_genotype(ref: str, alt: str, genotype: str) -> str:
    gt = (genotype or "").strip().replace("|", "/")
    if not gt or gt in {".", "./."}:
        return "-"
    if re.fullmatch(r"\d+(?:/\d+)?", gt):
        alleles = []
        for token in gt.split("/"):
            if token == "0":
                alleles.append(ref)
            elif token == "1":
                alleles.append(alt)
            else:
                alleles.append("?")
        return "/".join(alleles)
    compact = re.sub(r"[/\s]", "", gt.upper())
    if compact and set(compact) <= {"A", "C", "G", "T"} and len(compact) <= 2:
        return "/".join(compact)
    return gt


def condition_tokens(values: Iterable[str]) -> set[str]:
    tokens: set[str] = set()
    for value in values:
        text = re.sub(r"[_-]+", " ", (value or "").lower())
        for token in re.findall(r"[a-z0-9]+", text):
            if len(token) >= 3 and token not in STOPWORDS:
                tokens.add(token)
    return tokens


def flatten_scalars(value: Any, path: str = "$") -> Iterator[tuple[str, str]]:
    if isinstance(value, dict):
        for key, child in value.items():
            yield from flatten_scalars(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from flatten_scalars(child, f"{path}[{index}]")
    elif value is not None:
        yield path, str(value)


def walk_dicts(value: Any, path: str = "$") -> Iterator[tuple[str, dict[str, Any]]]:
    if isinstance(value, dict):
        yield path, value
        for key, child in value.items():
            yield from walk_dicts(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from walk_dicts(child, f"{path}[{index}]")


def field_values(node: dict[str, Any], key_terms: Iterable[str]) -> list[str]:
    terms = tuple(t.lower() for t in key_terms)
    result: list[str] = []
    for path, value in flatten_scalars(node):
        leaf = re.split(r"\.|\[", path)[-1].rstrip("]").lower()
        if any(term in leaf for term in terms):
            if value and value not in result:
                result.append(value)
    return result


def source_doc_id(node: dict[str, Any]) -> str:
    for path, value in flatten_scalars(node):
        leaf = re.split(r"\.|\[", path)[-1].rstrip("]").lower()
        if leaf in {"docid", "documentid", "reportid", "curationid", "id"} and value.strip():
            return value.strip()
    return ""


def age_from_text(text: str) -> str:
    patterns = [
        r"(?:by|at|before|after)\s+age\s+(\d{1,3})",
        r"age\s+(\d{1,3})",
        r"(?:by|at|before|after)\s+(\d{1,3})\s+years?",
    ]
    lower = text.lower()
    for pattern in patterns:
        match = re.search(pattern, lower)
        if match:
            return match.group(1)
    return ""


def metric_type(path: str, text: str, ratio_hint: str = "") -> str:
    joined = f"{path} {text}".lower()
    if ratio_hint:
        return ratio_hint
    if "penetrance" in joined:
        return "penetrance"
    if "lifetime" in joined and "risk" in joined:
        return "lifetime_risk"
    if "cumulative" in joined and "risk" in joined:
        return "cumulative_risk"
    if "incidence" in joined:
        return "incidence"
    if "prevalence" in joined:
        return "prevalence"
    if "probability" in joined:
        return "probability"
    if "frequency" in joined:
        return "frequency"
    return "absolute_risk"


def parse_metric(path: str, value: str) -> dict[str, str] | None:
    path_lower = path.lower().replace("_", "")
    text = (value or "").strip()
    lower = text.lower()
    if not text or any(token in path_lower for token in BLOCKED_CONTEXT):
        return None

    ratio_patterns = [
        ("relative_risk", r"(?:\brr\b|relative\s+risk)\s*(?:of|=|:)?\s*([0-9]+(?:\.[0-9]+)?)"),
        ("hazard_ratio", r"(?:\bhr\b|hazard\s+ratio)\s*(?:of|=|:)?\s*([0-9]+(?:\.[0-9]+)?)"),
        ("odds_ratio", r"(?:\bor\b|odds\s+ratio)\s*(?:of|=|:)?\s*([0-9]+(?:\.[0-9]+)?)"),
    ]
    for kind, pattern in ratio_patterns:
        match = re.search(pattern, lower, flags=re.IGNORECASE)
        if match:
            return {
                "metric_type": kind,
                "value": match.group(1),
                "unit": "ratio",
                "age": age_from_text(text),
                "evidence_text": text[:1200],
                "source_field": path,
            }

    allowed = any(keyword in path_lower or keyword in lower for keyword in RISK_KEYWORDS)
    if not allowed:
        return None

    range_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(?:-|–|—|to)\s*([0-9]+(?:\.[0-9]+)?)\s*%", text, flags=re.IGNORECASE)
    if range_match:
        value_text = f"{range_match.group(1)}–{range_match.group(2)}"
        return {
            "metric_type": metric_type(path, text),
            "value": value_text,
            "unit": "%",
            "age": age_from_text(text),
            "evidence_text": text[:1200],
            "source_field": path,
        }

    percent_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*%", text)
    if percent_match:
        return {
            "metric_type": metric_type(path, text),
            "value": percent_match.group(1),
            "unit": "%",
            "age": age_from_text(text),
            "evidence_text": text[:1200],
            "source_field": path,
        }

    fraction_match = re.search(r"\b1\s+(?:in|/|per)\s+([0-9]+(?:\.[0-9]+)?)\b", lower)
    if fraction_match:
        return {
            "metric_type": metric_type(path, text),
            "value": f"1 in {fraction_match.group(1)}",
            "unit": "fraction",
            "age": age_from_text(text),
            "evidence_text": text[:1200],
            "source_field": path,
        }

    # Intentionally reject bare numbers such as ClinGen likelihood score 0-3.
    return None


def context_label(context: str) -> str:
    return "Adult" if context.lower().startswith("adult") else "Pediatric"


def risk_label_ru(kind: str) -> str:
    return {
        "penetrance": "Пенетрантность",
        "lifetime_risk": "Пожизненный риск",
        "cumulative_risk": "Кумулятивный риск",
        "absolute_risk": "Абсолютный риск",
        "probability": "Вероятность",
        "incidence": "Заболеваемость/частота события",
        "prevalence": "Распространенность",
        "relative_risk": "RR",
        "hazard_ratio": "HR",
        "odds_ratio": "OR",
        "frequency": "Частота",
    }.get(kind, kind)


def render_metric_ru(metric: dict[str, Any]) -> str:
    unit = metric.get("unit", "")
    value = str(metric.get("value", ""))
    amount = f"{value}{unit}" if unit == "%" else value
    age = f" к возрасту {metric['age']} лет" if metric.get("age") else ""
    scope = "variant-specific" if metric.get("scope") == "variant-specific" else "gene-condition"
    context = metric.get("context", "")
    return f"{risk_label_ru(str(metric.get('metric_type','')))} {amount}{age} (ClinGen Actionability, {context}; {scope})"


def extract_snapshot_metrics(snapshot: dict[str, Any], known_genes: set[str], context: str) -> dict[str, list[dict[str, Any]]]:
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    seen: set[tuple[str, str, str, str, str, str, str]] = set()
    for node_path, node in walk_dicts(snapshot):
        serialized = json.dumps(node, ensure_ascii=False, separators=(",", ":")).upper()
        matched = [gene for gene in known_genes if re.search(rf"(?<![A-Z0-9]){re.escape(gene.upper())}(?![A-Z0-9])", serialized)]
        if not matched:
            continue
        conditions = unique(field_values(node, ("condition", "disease", "outcome", "syndrome", "phenotype")))
        doc_id = source_doc_id(node)
        for scalar_path, scalar_value in flatten_scalars(node, node_path):
            metric = parse_metric(scalar_path, scalar_value)
            if not metric:
                continue
            for gene in matched:
                for condition in (conditions or [""]):
                    key = (
                        gene,
                        metric["metric_type"], metric["value"], metric["unit"], metric.get("age", ""),
                        condition, metric["source_field"],
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    by_gene[gene].append({
                        **metric,
                        "gene": gene,
                        "condition": condition,
                        "context": context_label(context),
                        "source": "ClinGen Actionability",
                        "source_doc_id": doc_id,
                        "node_path": node_path,
                        "node_text_upper": serialized,
                    })
    return by_gene


def metric_score(metric: dict[str, Any], event_tokens: set[str], rsid: str) -> tuple[int, int, int, str, str]:
    condition = str(metric.get("condition", ""))
    candidate_tokens = condition_tokens([condition])
    overlap = len(event_tokens & candidate_tokens)
    scope = "variant-specific" if rsid and rsid.upper() in str(metric.get("node_text_upper", "")) else "gene-condition"
    metric["scope"] = scope
    return (
        100 if scope == "variant-specific" else 0,
        overlap * 10,
        METRIC_RANK.get(str(metric.get("metric_type", "")), 0),
        str(metric.get("context", "")),
        str(metric.get("source_field", "")),
    )


def choose_metrics(candidates: list[dict[str, Any]], event_conditions: list[str], rsid: str) -> list[dict[str, Any]]:
    event_tokens = condition_tokens(event_conditions)
    accepted: list[dict[str, Any]] = []
    for metric in candidates:
        candidate_tokens = condition_tokens([str(metric.get("condition", ""))])
        # A disease-specific quantitative number must have at least one semantic token
        # in common with the patient's Stage-24 condition. If either side lacks a
        # condition label, keep only an exact variant-specific metric.
        exact_variant = bool(rsid and rsid.upper() in str(metric.get("node_text_upper", "")))
        if event_tokens and candidate_tokens and not (event_tokens & candidate_tokens) and not exact_variant:
            continue
        if (not event_tokens or not candidate_tokens) and not exact_variant:
            continue
        metric_score(metric, event_tokens, rsid)
        accepted.append(metric)

    accepted.sort(key=lambda m: metric_score(m, event_tokens, rsid), reverse=True)
    result: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str, str, str]] = set()
    for metric in accepted:
        key = (
            str(metric.get("metric_type", "")), str(metric.get("value", "")), str(metric.get("unit", "")),
            str(metric.get("age", "")), str(metric.get("condition", "")), str(metric.get("context", "")),
        )
        if key in seen:
            continue
        seen.add(key)
        clean = {k: v for k, v in metric.items() if k != "node_text_upper"}
        clean["risk_summary_ru"] = render_metric_ru(clean)
        result.append(clean)
        if len(result) >= MAX_METRICS_PER_EVENT:
            break
    return result


def actionability_paths(root: Path, manifest: dict[str, Any]) -> dict[str, Path]:
    clin = manifest.get("sources", {}).get("ClinGen", {})
    artifacts = clin.get("artifacts", {}) if isinstance(clin, dict) else {}
    result: dict[str, Path] = {}
    for key, context in (("actionability_adult", "Adult"), ("actionability_pediatric", "Pediatric")):
        record = artifacts.get(key, {}) if isinstance(artifacts, dict) else {}
        rel = record.get("path") if isinstance(record, dict) else None
        if not rel:
            raise RuntimeError(f"Open-reference manifest lacks ClinGen/{key}")
        path = root / "reference_data" / "open_reference" / str(rel)
        if not path.is_file():
            raise FileNotFoundError(path)
        result[context] = path
    return result


def input_paths(root: Path, patient_root: Path) -> dict[str, Path]:
    knowledge = patient_root / "Dante_results" / "knowledge_enrichment"
    open_root = root / "reference_data" / "open_reference"
    return {
        "stage24_state": knowledge / "Dante_knowledge_enrichment.state.json",
        "open_findings": knowledge / "Dante_knowledge_findings_open.tsv",
        "open_reference_manifest": open_root / "open_reference_manifest.json",
    }


def output_paths(patient_root: Path) -> dict[str, Path]:
    out = patient_root / "Dante_results" / "clinical_risk_evidence"
    return {
        "dir": out,
        "tsv": out / "Dante_clinical_risk_evidence.tsv",
        "json": out / "Dante_clinical_risk_evidence.json",
        "state": out / "Dante_clinical_risk_evidence.state.json",
    }


def current_input_identities(root: Path, patient_root: Path) -> tuple[dict[str, dict[str, Any]], dict[str, Path], dict[str, Any]]:
    paths = input_paths(root, patient_root)
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    stage_state = read_json(paths["stage24_state"])
    if stage_state.get("status") != "complete" or stage_state.get("schema") != STAGE24_SCHEMA or stage_state.get("rules_schema") != RULES_SCHEMA:
        raise RuntimeError("Stage 24 is not READY under rules-v8")
    manifest = read_json(paths["open_reference_manifest"])
    if manifest.get("status") != "complete" or manifest.get("schema") != "open-reference-bundle-v2":
        raise RuntimeError("Open reference bundle is not READY")
    action_paths = actionability_paths(root, manifest)
    identities = {name: identity(path) for name, path in paths.items()}
    for context, path in action_paths.items():
        identities[f"clingen_actionability_{context.lower()}"] = identity(path)
    return identities, action_paths, stage_state


def is_ready(root: Path, patient_root: Path) -> tuple[bool, str]:
    out = output_paths(patient_root)
    if not out["state"].is_file() or not out["json"].is_file() or not out["tsv"].is_file():
        return False, "MISSING"
    try:
        state = read_json(out["state"])
        if state.get("schema") != STATE_SCHEMA or state.get("status") != "complete":
            return False, "INVALID STATE"
        inputs, _, _ = current_input_identities(root, patient_root)
        if state.get("inputs") != inputs:
            return False, "STALE INPUTS"
        outputs = state.get("outputs", {})
        for key in ("json", "tsv"):
            if outputs.get(key, {}).get("sha256") != sha256_file(out[key]):
                return False, f"STALE {key.upper()}"
        return True, "READY"
    except Exception as exc:
        return False, f"ERROR: {type(exc).__name__}: {exc}"


def event_groups(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row.get("event_key") or row.get("finding_id") or "").strip()
        if key:
            grouped[key].append(row)
    return dict(grouped)


def build(root: Path, patient_root: Path, force: bool = False) -> int:
    ready, status_text = is_ready(root, patient_root)
    if ready and not force:
        print("READY: clinical risk evidence is current")
        return 0

    inputs, action_paths, stage_state = current_input_identities(root, patient_root)
    in_paths = input_paths(root, patient_root)
    rows = read_tsv(in_paths["open_findings"])
    grouped = event_groups(rows)
    known_genes = {row.get("gene", "").strip().upper() for row in rows if row.get("gene", "").strip() not in {"", ".", "-"}}

    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for context, path in action_paths.items():
        snapshot = json.loads(path.read_text(encoding="utf-8"))
        extracted = extract_snapshot_metrics(snapshot, known_genes, context)
        for gene, metrics in extracted.items():
            by_gene[gene].extend(metrics)

    events: dict[str, dict[str, Any]] = {}
    tsv_rows: list[dict[str, Any]] = []
    metric_counts: Counter[str] = Counter()
    events_with_metrics = 0

    for event_key in sorted(grouped):
        event_rows = grouped[event_key]
        ref = next((r.get("ref", "") for r in event_rows if r.get("ref", "") not in {"", ".", "-"}), "")
        alt = next((r.get("alt", "") for r in event_rows if r.get("alt", "") not in {"", ".", "-"}), "")
        genotype = next((r.get("genotype", "") for r in event_rows if r.get("genotype", "") not in {"", ".", "-"}), "")
        zyg = next((r.get("zygosity", "") for r in event_rows if r.get("zygosity", "") not in {"", ".", "-"}), "")
        rsid = next((r.get("resolved_rsid", "") for r in event_rows if r.get("resolved_rsid", "").startswith("rs")), "")
        if not rsid:
            rsid = next((r.get("rsid", "") for r in event_rows if r.get("rsid", "").startswith("rs")), "")
        genes = unique(r.get("gene", "") for r in event_rows if r.get("gene", "") not in {"", ".", "-"})
        conditions = unique(
            item
            for r in event_rows
            for value in (r.get("clinvar_conditions", ""), r.get("disease_name", ""))
            for item in split_multi(value)
        )
        classification, raw_sig = normalize_clinvar_classification(r.get("clinvar_significance", "") for r in event_rows)
        stars = max([int(r.get("clinvar_review_stars", "0") or 0) for r in event_rows] or [0])
        candidates = [metric for gene in genes for metric in by_gene.get(gene.upper(), [])]
        metrics = choose_metrics(candidates, conditions, rsid)
        if metrics:
            events_with_metrics += 1
            for metric in metrics:
                metric_counts[str(metric.get("metric_type", "unknown"))] += 1
            risk_summary = "; ".join(str(m.get("risk_summary_ru", "")) for m in metrics)
        else:
            risk_summary = "Количественная оценка риска не установлена по локальному ClinGen Actionability snapshot"

        event = {
            "event_key": event_key,
            "rsid": rsid,
            "genes": genes,
            "ref": ref,
            "alt": alt,
            "genotype": genotype,
            "patient_alt_allele": alt,
            "patient_genotype_display": readable_genotype(ref, alt, genotype),
            "zygosity": zyg,
            "clinvar_classification": classification,
            "clinvar_significance_raw": raw_sig,
            "clinvar_review_stars": stars,
            "conditions": conditions,
            "metrics": metrics,
            "risk_summary_ru": risk_summary,
        }
        events[event_key] = event

        if metrics:
            for metric in metrics:
                tsv_rows.append({
                    "event_key": event_key, "rsid": rsid, "gene": metric.get("gene", ""), "ref": ref, "alt": alt,
                    "genotype": genotype, "patient_alt_allele": alt, "patient_genotype_display": event["patient_genotype_display"],
                    "zygosity": zyg, "clinvar_classification": classification, "clinvar_significance_raw": raw_sig,
                    "clinvar_review_stars": stars, "conditions": "; ".join(conditions),
                    "risk_metric_type": metric.get("metric_type", ""), "risk_value": metric.get("value", ""),
                    "risk_unit": metric.get("unit", ""), "risk_age": metric.get("age", ""),
                    "risk_condition": metric.get("condition", ""), "risk_scope": metric.get("scope", ""),
                    "risk_source": metric.get("source", ""), "risk_context": metric.get("context", ""),
                    "source_doc_id": metric.get("source_doc_id", ""), "source_field": metric.get("source_field", ""),
                    "risk_summary_ru": metric.get("risk_summary_ru", ""),
                })
        else:
            tsv_rows.append({
                "event_key": event_key, "rsid": rsid, "gene": ",".join(genes), "ref": ref, "alt": alt,
                "genotype": genotype, "patient_alt_allele": alt, "patient_genotype_display": event["patient_genotype_display"],
                "zygosity": zyg, "clinvar_classification": classification, "clinvar_significance_raw": raw_sig,
                "clinvar_review_stars": stars, "conditions": "; ".join(conditions),
                "risk_metric_type": "unavailable", "risk_source": "ClinGen Actionability",
                "risk_scope": "none", "risk_summary_ru": risk_summary,
            })

    out = output_paths(patient_root)
    payload = {
        "schema": SCHEMA,
        "status": "complete",
        "rules_schema": RULES_SCHEMA,
        "source_policy": SOURCE_POLICY,
        "network_used": False,
        "patient_identifiers_sent_external": False,
        "interpretation_boundary": (
            "Quantitative evidence is publication context, not a diagnosis. Gene-condition estimates are not variant-specific "
            "unless the source record explicitly contains the patient's rsID. Bare ClinGen likelihood scores are not risk estimates."
        ),
        "events": events,
    }
    write_json_atomic(out["json"], payload)
    write_tsv_atomic(out["tsv"], tsv_rows)
    state = {
        "schema": STATE_SCHEMA,
        "status": "complete",
        "generated_utc": utcnow(),
        "rules_schema": RULES_SCHEMA,
        "source_policy": SOURCE_POLICY,
        "network_used": False,
        "patient_identifiers_sent_external": False,
        "input_profile": stage_state.get("input_profile", {}),
        "inputs": inputs,
        "outputs": {"json": identity(out["json"]), "tsv": identity(out["tsv"])},
        "counts": {
            "stage24_open_rows": len(rows),
            "events": len(events),
            "events_with_quantitative_metric": events_with_metrics,
            "events_without_quantitative_metric": len(events) - events_with_metrics,
            "metrics_by_type": dict(sorted(metric_counts.items())),
        },
    }
    write_json_atomic(out["state"], state)
    print("CLINICAL RISK EVIDENCE")
    print(f"  Stage-24 events:                {len(events)}")
    print(f"  Events with quantitative risk: {events_with_metrics}")
    print(f"  Unavailable quantitative risk: {len(events) - events_with_metrics}")
    print("  Source:                         local ClinGen Actionability snapshots")
    print("  Network/API:                    NOT USED")
    print(f"  JSON:                           {out['json']}")
    print(f"  TSV:                            {out['tsv']}")
    return 0


def self_test() -> int:
    assert normalize_clinvar_classification(["Pathogenic"])[0] == "Pathogenic"
    assert normalize_clinvar_classification(["Likely_pathogenic"])[0] == "Likely pathogenic"
    assert normalize_clinvar_classification(["Pathogenic|Likely_pathogenic"])[0] == "Pathogenic / Likely pathogenic"
    assert readable_genotype("A", "G", "0/1") == "A/G"
    assert readable_genotype("A", "G", "1/1") == "G/G"
    assert parse_metric("$.likelihood", "3") is None
    m = parse_metric("$.penetrance", "Penetrance is approximately 70% by age 70")
    assert m and m["metric_type"] == "penetrance" and m["value"] == "70" and m["age"] == "70"
    m = parse_metric("$.likelihood", "Lifetime risk 28-76%")
    assert m and m["metric_type"] == "lifetime_risk" and m["value"] == "28–76"
    m = parse_metric("$.risk", "RR 2.4 (95% CI 1.2-4.8)")
    assert m and m["metric_type"] == "relative_risk" and m["value"] == "2.4"
    assert parse_metric("$.effectiveness", "Risk reduced by 60%") is None

    snapshot = {
        "rows": [
            {"geneOrVariant": "BRAF", "condition": "Noonan syndrome", "likelihood": "Penetrance 80% by age 40", "docId": "AC1"},
            {"geneOrVariant": "BRAF", "condition": "Unrelated condition", "likelihood": "Risk 90%", "docId": "AC2"},
        ]
    }
    by_gene = extract_snapshot_metrics(snapshot, {"BRAF"}, "Adult")
    chosen = choose_metrics(by_gene["BRAF"], ["Noonan syndrome"], "rs123")
    assert chosen and chosen[0]["value"] == "80" and chosen[0]["condition"] == "Noonan syndrome"
    assert all(item["condition"] != "Unrelated condition" for item in chosen)
    print("SELF-TEST PASS | clinical-risk evidence classification, genotype and conservative quantitative extraction")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("build", "status", "self-test"))
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--patient-root", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.mode == "self-test":
        return self_test()
    if not args.patient_root:
        raise RuntimeError("--patient-root is required")
    root = args.root.resolve()
    patient_root = args.patient_root.resolve()
    if args.mode == "status":
        ready, text = is_ready(root, patient_root)
        print(text)
        return 0 if ready else 1
    return build(root, patient_root, args.force)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}")
        raise SystemExit(1)
