"""Context-scoped extraction fix for clinical_risk_evidence.

The base v1 recursive scanner is intentionally schema-tolerant, but a parent JSON
object containing several sibling gene-condition records can pair a number from one
sibling with the condition label from another. This overlay replaces only snapshot
metric extraction with an inherited-context traversal: parent metadata flows down,
but sibling metadata never flows sideways.
"""
from __future__ import annotations

import json
import re
from collections import defaultdict
from typing import Any, Iterable

import clinical_risk_evidence as base


def _scalar_strings(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (str, int, float, bool)):
        return [str(value)]
    if isinstance(value, list) and all(isinstance(x, (str, int, float, bool)) or x is None for x in value):
        return [str(x) for x in value if x is not None]
    return []


def _leaf_name(key: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (key or "").lower())


def _local_context(
    node: dict[str, Any],
    known_genes: set[str],
    inherited_genes: set[str],
    inherited_conditions: list[str],
    inherited_doc_id: str,
    inherited_text: str,
) -> tuple[set[str], list[str], str, str]:
    direct_pairs: list[tuple[str, str]] = []
    for key, value in node.items():
        for scalar in _scalar_strings(value):
            direct_pairs.append((str(key), scalar))

    direct_text = " ".join(value for _, value in direct_pairs)
    genes = set(inherited_genes)
    for gene in known_genes:
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])", direct_text, flags=re.IGNORECASE):
            genes.add(gene)

    conditions = list(inherited_conditions)
    for key, value in direct_pairs:
        leaf = _leaf_name(key)
        if any(token in leaf for token in ("condition", "disease", "outcome", "syndrome", "phenotype")):
            for item in base.split_multi(value):
                if item not in conditions:
                    conditions.append(item)

    doc_id = inherited_doc_id
    for key, value in direct_pairs:
        leaf = _leaf_name(key)
        if leaf in {"docid", "documentid", "reportid", "curationid", "id"} and value.strip():
            doc_id = value.strip()
            break

    context_text = (inherited_text + " " + direct_text).strip().upper()
    return genes, conditions, doc_id, context_text


def extract_snapshot_metrics(snapshot: dict[str, Any], known_genes: set[str], context: str) -> dict[str, list[dict[str, Any]]]:
    known = {str(g).strip().upper() for g in known_genes if str(g).strip()}
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    seen: set[tuple[str, str, str, str, str, str, str, str]] = set()

    def walk(
        value: Any,
        path: str,
        genes: set[str],
        conditions: list[str],
        doc_id: str,
        context_text: str,
    ) -> None:
        if isinstance(value, dict):
            local_genes, local_conditions, local_doc, local_text = _local_context(
                value, known, genes, conditions, doc_id, context_text
            )

            # Parse only scalar values owned by this dictionary. Nested child values
            # are handled by recursive calls with inherited context; this is the key
            # property that prevents cross-sibling condition leakage.
            for key, child in value.items():
                scalars = _scalar_strings(child)
                if not scalars:
                    continue
                for index, scalar in enumerate(scalars):
                    scalar_path = f"{path}.{key}" + (f"[{index}]" if len(scalars) > 1 else "")
                    metric = base.parse_metric(scalar_path, scalar)
                    if not metric or not local_genes:
                        continue
                    for gene in sorted(local_genes):
                        for condition in (local_conditions or [""]):
                            dedup = (
                                gene,
                                metric["metric_type"],
                                metric["value"],
                                metric["unit"],
                                metric.get("age", ""),
                                condition,
                                metric["source_field"],
                                local_doc,
                            )
                            if dedup in seen:
                                continue
                            seen.add(dedup)
                            by_gene[gene].append(
                                {
                                    **metric,
                                    "gene": gene,
                                    "condition": condition,
                                    "context": base.context_label(context),
                                    "source": "ClinGen Actionability",
                                    "source_doc_id": local_doc,
                                    "node_path": path,
                                    "node_text_upper": local_text,
                                }
                            )

            for key, child in value.items():
                if isinstance(child, (dict, list)) and not _scalar_strings(child):
                    walk(child, f"{path}.{key}", local_genes, local_conditions, local_doc, local_text)
            return

        if isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, f"{path}[{index}]", genes, conditions, doc_id, context_text)

    walk(snapshot, "$", set(), [], "", "")
    return dict(by_gene)


# Patch the implementation module itself so build()/self_test()/status consumers all
# use the corrected function even when they call through base.main().
base.extract_snapshot_metrics = extract_snapshot_metrics


# Re-export the stable public API.
for _name in dir(base):
    if _name.startswith("__"):
        continue
    globals().setdefault(_name, getattr(base, _name))

# Ensure our replacement wins over the re-exported original name.
globals()["extract_snapshot_metrics"] = extract_snapshot_metrics


def self_test() -> int:
    # Retain all base contracts first, but under the corrected extractor.
    base.self_test()
    snapshot = {
        "rows": [
            {"gene": "BRAF", "condition": "Noonan syndrome", "risk": "Penetrance 80% by age 40", "docId": "A"},
            {"gene": "BRAF", "condition": "Melanoma", "risk": "Lifetime risk 90%", "docId": "B"},
        ]
    }
    by_gene = extract_snapshot_metrics(snapshot, {"BRAF"}, "Adult")
    selected = base.choose_metrics(by_gene["BRAF"], ["Noonan syndrome"], "rs-does-not-appear")
    assert len(selected) == 1
    assert selected[0]["condition"] == "Noonan syndrome"
    assert selected[0]["value"] == "80"
    print("SELF-TEST PASS | sibling disease contexts remain isolated")
    return 0


def main() -> int:
    # base.main resolves base.extract_snapshot_metrics dynamically and therefore uses
    # the patched implementation above.
    if len(__import__("sys").argv) > 1 and __import__("sys").argv[1] == "self-test":
        return self_test()
    return base.main()


if __name__ == "__main__":
    raise SystemExit(main())
