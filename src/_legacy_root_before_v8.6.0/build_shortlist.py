#!/usr/bin/env python3
"""Create an evidence-classed shortlist from the unified Dante catalog.

v8.0 adds a genome-wide uncertain-coding sweep:
- scans every carried PASS variant in the unified catalog, not a hand-picked gene list;
- retains ClinVar VUS coding/splice variants without requiring SIFT/PolyPhen damage calls;
- retains rare/unknown-frequency unclassified coding/splice variants (default max AF 1%);
- preserves generic ClinVar P/LP/conflict/pharmacogenomic and rare predicted-impact rules;
- enriches recognized gene symbols from the embedded CHESS 2.2 reference.

CHESS coordinates are metadata only; GRCh37 gene assignment remains the VEP annotation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
import traceback
from collections import Counter
from pathlib import Path
from typing import TextIO

RULES_SCHEMA = "deterministic-shortlist-v3-chess-genomewide"

PREFIX_COLUMNS = [
    "priority_rank", "reason_code", "evidence_class", "triage_reason",
    "genomewide_uncertain_sweep", "genomewide_reason_code",
    "clinical_status", "functional_class",
    "chess_gene_match", "chess_gene_symbols", "chess_gene_types",
    "chess_descriptions", "chess_gff_ids", "chess_source_locations",
    "chess_databases", "chess_refseq_geneids",
]

LOF_TERMS = (
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "start_lost",
)
SEVERE_OTHER_TERMS = ("stop_lost",)
INFRAME_TERMS = ("inframe_insertion", "inframe_deletion", "protein_altering_variant")
SPLICE_REGION_TERMS = (
    "splice_region_variant", "splice_donor_5th_base_variant",
    "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
)
UNCERTAIN_FUNCTIONAL_CLASSES = {
    "LOF", "MISSENSE", "INFRAME_OR_PROTEIN_ALTERING", "SPLICE_REGION", "SEVERE_OTHER"
}


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block_size):
            digest.update(chunk)
    return digest.hexdigest()


def max_frequency(value: str) -> float | None:
    numbers: list[float] = []
    for token in re.split(r"[,&|]", value or ""):
        try:
            numbers.append(float(token))
        except ValueError:
            pass
    return max(numbers) if numbers else None


def contains_any(value: str, terms: tuple[str, ...]) -> bool:
    lowered = (value or "").lower().replace("_", " ")
    return any(term in lowered for term in terms)


def significance_text(row: dict[str, str]) -> str:
    current = row.get("clinvar_current_significance", "")
    cached = row.get("clin_sig", "")
    return current if current not in {"", "."} else cached


def clinical_status(row: dict[str, str]) -> str:
    sig = significance_text(row)
    review = row.get("clinvar_current_review_status", "")
    pathogenic = contains_any(sig, ("pathogenic",))
    benign = contains_any(sig, ("benign",))
    conflicting = contains_any(sig, ("conflict",)) or contains_any(review, ("conflict",)) or (pathogenic and benign)
    if conflicting:
        return "CONFLICTING"
    if benign:
        return "BENIGN_OR_LIKELY_BENIGN"
    if pathogenic:
        return "PATHOGENIC_OR_LIKELY_PATHOGENIC"
    if contains_any(sig, ("uncertain significance", "vus")):
        return "VUS"
    normalized = (sig or "").strip().lower().replace("_", " ")
    if normalized in {"", ".", "not provided", "not specified", "no assertion provided"}:
        return "NO_ASSERTION"
    if contains_any(sig, ("drug response", "pharmacogen")):
        return "DRUG_RESPONSE"
    if contains_any(sig, ("risk factor", "association")):
        return "RISK_ASSOCIATION"
    if contains_any(sig, ("protective",)):
        return "PROTECTIVE"
    return "OTHER"


def consequence_terms(row: dict[str, str]) -> set[str]:
    value = row.get("all_consequences", row.get("consequence", ""))
    return {term.strip().lower() for term in re.split(r"[,;&|]", value or "") if term.strip() and term.strip() != "."}


def functional_class(row: dict[str, str]) -> str:
    terms = consequence_terms(row)
    if any(term in terms for term in LOF_TERMS):
        return "LOF"
    if "missense_variant" in terms or any("missense" in term for term in terms):
        return "MISSENSE"
    if any(term in terms for term in INFRAME_TERMS):
        return "INFRAME_OR_PROTEIN_ALTERING"
    if any(term in terms for term in SPLICE_REGION_TERMS):
        return "SPLICE_REGION"
    if any(term in terms for term in SEVERE_OTHER_TERMS):
        return "SEVERE_OTHER"
    return "OTHER"


def row_gene_symbols(row: dict[str, str]) -> set[str]:
    values = [row.get("symbol", ""), row.get("all_symbols", "")]
    genes: set[str] = set()
    for value in values:
        for token in re.split(r"[,;&|]", value or ""):
            token = token.strip().upper()
            if token and token != ".":
                genes.add(token)
    return genes


def load_chess_reference(path: Path) -> tuple[dict[str, list[dict[str, str]]], dict[str, int]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    index: dict[str, list[dict[str, str]]] = {}
    stats: Counter[str] = Counter()
    with path.open("rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Gene_Name", "Description", "Gene_Type", "GFF_ID", "Location", "Database", "RefSeq_GeneID"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise RuntimeError("CHESS gene reference has an unexpected header")
        for row in reader:
            stats["rows"] += 1
            symbol = (row.get("Gene_Name") or "").strip().upper()
            if not symbol or symbol == "-":
                stats["unnamed_rows"] += 1
                continue
            index.setdefault(symbol, []).append(dict(row))
            stats["named_rows"] += 1
            stats[f"type:{row.get('Gene_Type','')}"] += 1
    stats["unique_named_genes"] = len(index)
    return index, dict(stats)


def chess_metadata(symbols: set[str], reference: dict[str, list[dict[str, str]]]) -> dict[str, str]:
    matched = sorted(symbol for symbol in symbols if symbol in reference)
    records = [record for symbol in matched for record in reference[symbol]]
    def joined(field: str) -> str:
        values = []
        seen = set()
        for record in records:
            value = (record.get(field) or "").strip()
            if value and value != "-" and value not in seen:
                seen.add(value)
                values.append(value)
        return " | ".join(values)
    return {
        "chess_gene_match": "yes" if matched else "no",
        "chess_gene_symbols": ",".join(matched),
        "chess_gene_types": joined("Gene_Type"),
        "chess_descriptions": joined("Description"),
        "chess_gff_ids": joined("GFF_ID"),
        "chess_source_locations": joined("Location"),
        "chess_databases": joined("Database"),
        "chess_refseq_geneids": joined("RefSeq_GeneID"),
    }


def classify_generic(row: dict[str, str]) -> tuple[int, str, str, str] | None:
    if row.get("alt_carried") != "yes" or row.get("vcf_filter") != "PASS":
        return None
    sig = significance_text(row)
    review = row.get("clinvar_current_review_status", "")
    try:
        stars = int(row.get("clinvar_review_stars", "0") or 0)
    except ValueError:
        stars = 0
    pathogenic = contains_any(sig, ("pathogenic",))
    benign = contains_any(sig, ("benign",))
    conflicting = contains_any(sig, ("conflict",)) or contains_any(review, ("conflict",)) or (pathogenic and benign)
    if conflicting:
        return 10, "CLINVAR_CONFLICT", "clinical_database_conflict", "ClinVar interpretations conflict; manual evidence review required"
    if pathogenic and not benign and stars >= 2:
        return 1, "CLINVAR_PLP_STRONG", "clinical_database_assertion", f"ClinVar P/LP assertion with {stars}-star review status"
    if pathogenic and not benign:
        return 2, "CLINVAR_PLP_LIMITED", "clinical_database_assertion_limited", f"ClinVar P/LP assertion with limited review status ({stars} stars)"
    if contains_any(sig, ("drug response", "pharmacogen")):
        return 20, "CLINVAR_DRUG_RESPONSE", "pharmacogenomic_annotation", "ClinVar drug-response/pharmacogenomic annotation"
    if contains_any(sig, ("risk factor", "association")):
        return 30, "CLINVAR_RISK_ASSOCIATION", "risk_association", "ClinVar risk-factor/association annotation; not a causal diagnosis"
    if contains_any(sig, ("protective",)):
        return 31, "CLINVAR_PROTECTIVE", "protective_association", "ClinVar protective association"

    af = max_frequency(row.get("max_af", ""))
    rare = af is None or af <= 0.01
    impact = row.get("impact", "").upper()
    consequence = row.get("all_consequences", row.get("consequence", "")).lower()
    if impact == "HIGH" and rare and not benign:
        return 40, "RARE_PREDICTED_LOF", "predicted_molecular_consequence", "Rare/unknown-frequency HIGH-impact consequence; not a clinical classification"
    predicted_damaging = contains_any(row.get("sift", ""), ("deleterious",)) or contains_any(row.get("polyphen", ""), ("damaging", "probably damaging", "possibly damaging"))
    if impact == "MODERATE" and "missense" in consequence and rare and predicted_damaging and not benign:
        return 50, "RARE_MISSENSE_PREDICTED", "computational_support_only", "Rare/unknown-frequency missense with damaging computational prediction"
    return None


def classify_genomewide_uncertain(row: dict[str, str], unclassified_max_af: float) -> tuple[int, str, str, str] | None:
    """High-recall coding/splice uncertainty sweep across the entire annotated catalog."""
    if row.get("alt_carried") != "yes" or row.get("vcf_filter") != "PASS":
        return None
    status = clinical_status(row)
    if status in {"BENIGN_OR_LIKELY_BENIGN", "PATHOGENIC_OR_LIKELY_PATHOGENIC", "CONFLICTING", "DRUG_RESPONSE", "RISK_ASSOCIATION", "PROTECTIVE"}:
        return None
    func = functional_class(row)
    if func not in UNCERTAIN_FUNCTIONAL_CLASSES:
        return None

    explicit_vus = status == "VUS"
    if not explicit_vus:
        af = max_frequency(row.get("max_af", ""))
        if af is not None and af > unclassified_max_af:
            return None
        if status not in {"NO_ASSERTION", "OTHER"}:
            return None

    status_code = "VUS" if explicit_vus else "UNCLASSIFIED"
    ranks = {
        "LOF": (34, 36),
        "SEVERE_OTHER": (37, 38),
        "SPLICE_REGION": (44, 46),
        "MISSENSE": (54, 56),
        "INFRAME_OR_PROTEIN_ALTERING": (57, 58),
    }
    rank = ranks[func][0 if explicit_vus else 1]
    code = f"GENOME_{status_code}_{func}"
    evidence = "genomewide_clinvar_vus" if explicit_vus else "genomewide_rare_unclassified"
    label = {
        "LOF": "loss-of-function",
        "SEVERE_OTHER": "severe coding",
        "SPLICE_REGION": "splice-region",
        "MISSENSE": "missense",
        "INFRAME_OR_PROTEIN_ALTERING": "in-frame/protein-altering",
    }[func]
    if explicit_vus:
        reason = f"Genome-wide {label} variant classified as ClinVar VUS; included without requiring damaging SIFT/PolyPhen"
    else:
        reason = f"Genome-wide rare/unknown-frequency {label} variant without a definitive clinical classification; included without requiring damaging SIFT/PolyPhen"
    return rank, code, evidence, reason


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build deterministic + genome-wide uncertain-coding variant shortlist")
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genomewide-output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--gene-reference", type=Path, required=True)
    parser.add_argument("--unclassified-max-af", type=float, default=0.01)
    parser.add_argument("--progress-every", type=int, default=500_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = args.input.resolve()
    gene_reference_path = args.gene_reference.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if not 0 <= args.unclassified_max_af <= 1:
        raise ValueError("--unclassified-max-af must be between 0 and 1")
    chess, chess_stats = load_chess_reference(gene_reference_path)
    chess_hash = sha256_file(gene_reference_path)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.genomewide_output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    source_hash = sha256_file(source)
    state_path = args.output.with_name(args.output.name + ".state.json")
    if state_path.exists() and args.output.exists() and args.genomewide_output.exists():
        try:
            state = json.loads(state_path.read_text(encoding="utf-8"))
            if (
                state.get("catalog_sha256") == source_hash
                and state.get("status") == "complete"
                and state.get("rules_schema") == RULES_SCHEMA
                and state.get("gene_reference_sha256") == chess_hash
                and float(state.get("unclassified_max_af", -1)) == float(args.unclassified_max_af)
            ):
                print("Deterministic/genome-wide shortlist is current for the same catalog + CHESS reference; reusing it.")
                return 0
        except (OSError, ValueError, TypeError):
            pass

    selected: list[dict[str, str]] = []
    genomewide_rows: list[dict[str, str]] = []
    counts: Counter[str] = Counter()
    with open_text(source) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError("Catalog has no TSV header")
        source_columns = list(reader.fieldnames)
        for row in reader:
            counts["catalog_rows"] += 1
            if args.progress_every and counts["catalog_rows"] % args.progress_every == 0:
                print(f"  Reviewed {counts['catalog_rows']:,} catalog rows; shortlist {len(selected):,}; genome-wide uncertain {len(genomewide_rows):,}", flush=True)

            generic = classify_generic(row)
            genome = classify_genomewide_uncertain(row, args.unclassified_max_af)
            symbols = row_gene_symbols(row)
            metadata = chess_metadata(symbols, chess)
            if metadata["chess_gene_match"] == "yes":
                counts["catalog_rows_chess_match"] += 1

            if generic is None and genome is None:
                continue
            result = generic
            if genome is not None and (result is None or genome[0] < result[0]):
                result = genome
            assert result is not None
            rank, reason_code, evidence_class, reason = result
            base_row = {
                "priority_rank": str(rank),
                "reason_code": reason_code,
                "evidence_class": evidence_class,
                "triage_reason": reason,
                "genomewide_uncertain_sweep": "yes" if genome is not None else "no",
                "genomewide_reason_code": genome[1] if genome is not None else "",
                "clinical_status": clinical_status(row),
                "functional_class": functional_class(row),
                **metadata,
                **row,
            }
            selected.append(base_row)
            counts[reason_code] += 1
            if genome is not None:
                genomewide_rows.append(dict(base_row))
                counts["genomewide_rows"] += 1
                counts[f"genome_reason:{genome[1]}"] += 1
                counts[f"genome_func:{functional_class(row)}"] += 1
                if clinical_status(row) == "VUS":
                    counts["genome_explicit_vus"] += 1
                else:
                    counts["genome_unclassified"] += 1
                if metadata["chess_gene_match"] == "yes":
                    counts["genome_chess_match"] += 1

    sort_key = lambda row: (
        int(row["priority_rank"]), row.get("symbol", "."), row.get("chrom", ""),
        int(row.get("position", "0") or 0), row.get("alt", ""),
    )
    selected.sort(key=sort_key)
    genomewide_rows.sort(key=sort_key)
    fields = [*PREFIX_COLUMNS, *source_columns, "phenotype_or_question"]

    def publish(path: Path, records: list[dict[str, str]]) -> None:
        part = path.with_name(path.name + ".part")
        try:
            with part.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
                writer.writeheader()
                for record in records:
                    record["phenotype_or_question"] = ""
                    writer.writerow(record)
            os.replace(part, path)
        except Exception:
            part.unlink(missing_ok=True)
            raise

    publish(args.output, selected)
    publish(args.genomewide_output, genomewide_rows)

    report_lines = [
        "DANTE DETERMINISTIC + GENOME-WIDE UNCERTAIN-CODING SHORTLIST REPORT",
        "=" * 68,
        f"Rules schema: {RULES_SCHEMA}",
        f"Catalog rows reviewed: {counts['catalog_rows']:,}",
        f"Main shortlist rows: {len(selected):,}",
        f"Genome-wide uncertain-coding/splice rows: {len(genomewide_rows):,}",
        f"  explicit ClinVar VUS: {counts['genome_explicit_vus']:,}",
        f"  unclassified rare/unknown-frequency: {counts['genome_unclassified']:,}",
        "",
        "CHESS 2.2 GENE REFERENCE",
        f"Reference: {gene_reference_path}",
        f"Reference SHA256: {chess_hash}",
        f"Reference rows: {chess_stats.get('rows',0):,}",
        f"Unique named genes: {chess_stats.get('unique_named_genes',0):,}",
        f"Genome-wide candidates matching CHESS symbols: {counts['genome_chess_match']:,}",
        "CHESS source Location is metadata only; GRCh37 assignment is taken from VEP, so assemblies are not mixed.",
        "",
        "GENOME-WIDE UNCERTAIN-CODING RULE",
        "- No hand-picked gene list is used: every carried PASS catalog variant is evaluated.",
        "- ClinVar VUS: LoF/nonsense/frameshift/canonical splice, missense, in-frame/protein-altering, stop_lost and splice-region are retained regardless of AF.",
        f"- No assertion/other: the same molecular classes are retained when max AF is unknown or <= {args.unclassified_max_af:g}.",
        "- Benign/likely-benign variants are not promoted by this rule.",
        "- SIFT/PolyPhen are annotations/supporting evidence, not inclusion gates for this sweep.",
        "- Synonymous and ordinary intronic/intergenic variants are not included by this uncertainty sweep unless another generic rule selects them.",
        "",
        "COUNTS BY GENOME-WIDE MOLECULAR CLASS",
        *[f"{key.split(':',1)[1]}: {value:,}" for key,value in sorted(counts.items()) if key.startswith('genome_func:')],
        "",
        "COUNTS BY MAIN SHORTLIST REASON CODE",
        *[f"{key}: {value:,}" for key,value in sorted(counts.items()) if key.startswith(('CLINVAR_', 'RARE_', 'GENOME_'))],
        "",
        "COUNTS BY GENOME-WIDE SWEEP REASON CODE",
        *[f"{key.split(':',1)[1]}: {value:,}" for key,value in sorted(counts.items()) if key.startswith('genome_reason:')],
        "",
        "INTERPRETATION BOUNDARIES",
        "- VUS is a clinical-database classification; NO_ASSERTION/UNCLASSIFIED means no definitive classification was found.",
        "- Missense/nonsense/frameshift/splice are molecular consequences, not diagnoses.",
        "- CHESS membership or gene type is not treated as pathogenicity evidence.",
        "- This shortlist is triage for evidence review, not a diagnosis or treatment recommendation.",
        "",
    ]
    args.report.write_text("\n".join(report_lines), encoding="utf-8")
    state_path.write_text(json.dumps({
        "status": "complete",
        "rules_schema": RULES_SCHEMA,
        "catalog_sha256": source_hash,
        "gene_reference_sha256": chess_hash,
        "gene_reference_rows": chess_stats.get("rows", 0),
        "gene_reference_named_genes": chess_stats.get("unique_named_genes", 0),
        "unclassified_max_af": args.unclassified_max_af,
        "rows": len(selected),
        "genomewide_rows": len(genomewide_rows),
        "genomewide_output": args.genomewide_output.name,
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Created {args.output}: {len(selected):,} candidates; genome-wide uncertain-coding: {len(genomewide_rows):,}")
    return 0


def run_cli() -> None:
    try:
        code = main()
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 1
    except KeyboardInterrupt:
        print("\nERROR: cancelled by user.", file=sys.stderr)
        code = 130
    except BaseException as exc:
        print(f"\nERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        traceback.print_exc()
        code = 1
    raise SystemExit(code)


if __name__ == "__main__":
    run_cli()
