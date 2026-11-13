#!/usr/bin/env python3
"""Convert an Ensembl VEP VCF into practical Dante WGS outputs.

The script is deliberately deterministic and local. It:

* copies dbSNP rs identifiers reported by VEP Existing_variation into VCF ID;
* creates one compact TSV row per ALT allele while retaining the original VEP
  VCF with every transcript-level CSQ annotation;
* creates a conservative shortlist for later evidence review.

It expects VEP to have been run with --vcf --everything --allele_number.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import re
import sys
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from typing import BinaryIO, Iterable, TextIO
from urllib.parse import unquote


RSID_RE = re.compile(r"(?<![A-Za-z0-9])rs\d+", re.IGNORECASE)
CSQ_HEADER_RE = re.compile(r"Format:\s*([^\"]+)")
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}

SUMMARY_COLUMNS = [
    "variant_key",
    "source_layer",
    "chrom",
    "position",
    "ref",
    "alt",
    "allele_number",
    "rsid",
    "original_id",
    "vcf_filter",
    "GT",
    "genotype_bases",
    "alt_carried",
    "DP",
    "GQ",
    "AD",
    "allele_fraction",
    "impact",
    "consequence",
    "all_consequences",
    "symbol",
    "all_symbols",
    "gene",
    "feature_type",
    "transcript",
    "biotype",
    "canonical",
    "mane_select",
    "mane_plus_clinical",
    "hgvsc",
    "hgvsp",
    "protein_position",
    "amino_acids",
    "codons",
    "sift",
    "polyphen",
    "clin_sig",
    "phenotypes",
    "max_af",
    "max_af_pops",
    "pubmed",
    "variant_class",
    "csq_count",
    "source_file",
    "genome_build",
]

PRIORITY_COLUMNS = [
    "priority_tier",
    "priority_reason",
    *SUMMARY_COLUMNS,
    "phenotype_or_question",
]


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def open_atomic(path: Path) -> tuple[Path, TextIO, BinaryIO | None]:
    part = path.with_name(path.name + ".part")
    if path.name.lower().endswith(".gz"):
        raw_handle = part.open("wb")
        gzip_handle = gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0)
        handle = io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="")
    else:
        raw_handle = None
        handle = part.open("wt", encoding="utf-8", newline="")
    return part, handle, raw_handle


def finish_atomic(part: Path, handle: TextIO, raw_handle: BinaryIO | None, target: Path) -> None:
    handle.close()
    if raw_handle is not None:
        raw_handle.close()
    os.replace(part, target)


def parse_info(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    if not value or value == ".":
        return result
    for token in value.split(";"):
        if "=" in token:
            key, val = token.split("=", 1)
            result[key] = val
        elif token:
            result[token] = "true"
    return result


def parse_sample(fmt: str, value: str) -> dict[str, str]:
    if not fmt or fmt == ".":
        return {}
    keys = fmt.split(":")
    vals = value.split(":") if value else []
    return {key: vals[i] if i < len(vals) else "." for i, key in enumerate(keys)}


def genotype_bases(ref: str, alts: list[str], gt: str) -> str:
    if not gt or gt in {".", "./.", ".|."}:
        return "."
    sep = "|" if "|" in gt else "/"
    alleles = [ref, *alts]
    rendered: list[str] = []
    for token in re.split(r"[/|]", gt):
        if token == ".":
            rendered.append(".")
        elif token.isdigit() and int(token) < len(alleles):
            rendered.append(alleles[int(token)])
        else:
            rendered.append("?")
    return sep.join(rendered)


def called_alt_numbers(gt: str) -> set[int]:
    return {int(x) for x in re.split(r"[/|]", gt or "") if x.isdigit() and int(x) > 0}


def split_terms(value: str) -> list[str]:
    if not value:
        return []
    return [x for x in re.split(r"[&,]", value) if x]


def unique_join(values: Iterable[str], sep: str = "&") -> str:
    cleaned = sorted({unquote(v) for v in values if v and v not in {".", "-"}})
    return sep.join(cleaned) if cleaned else "."


def rsids_from(values: Iterable[str]) -> list[str]:
    found: set[str] = set()
    for value in values:
        for match in RSID_RE.findall(value or ""):
            found.add(match.lower())
    return sorted(found, key=lambda x: int(x[2:]))


def csq_record(raw: str, fields: list[str]) -> dict[str, str]:
    values = raw.split("|")
    if len(values) < len(fields):
        values.extend([""] * (len(fields) - len(values)))
    return {name: values[i] if i < len(values) else "" for i, name in enumerate(fields)}


def best_csq(records: list[dict[str, str]]) -> dict[str, str]:
    if not records:
        return {}

    def score(record: dict[str, str]) -> tuple[int, int, int, int, str]:
        mane = 2 if record.get("MANE_PLUS_CLINICAL") else (1 if record.get("MANE_SELECT") else 0)
        canonical = 1 if record.get("CANONICAL", "").upper() == "YES" else 0
        transcript = 1 if record.get("Feature_type") == "Transcript" else 0
        impact = IMPACT_RANK.get(record.get("IMPACT", "").upper(), 0)
        return mane, canonical, impact, transcript, record.get("Feature", "")

    return max(records, key=score)


def alt_specific(value: str, allele_number: int) -> str:
    """Select an allele-specific FORMAT value when possible.

    AD is Number=R (REF plus ALT), while AF is usually Number=A. Unknown shapes
    are preserved rather than guessed.
    """

    if not value or value == ".":
        return "."
    parts = value.split(",")
    if len(parts) > allele_number:
        return parts[allele_number]
    if len(parts) == 1:
        return parts[0]
    return value


def max_af_number(value: str) -> float | None:
    if not value or value in {".", "-"}:
        return None
    numbers: list[float] = []
    for token in re.split(r"[&,]", value):
        try:
            numbers.append(float(token))
        except ValueError:
            pass
    return max(numbers) if numbers else None


def priority_for(row: dict[str, str]) -> tuple[str, str] | None:
    if row["vcf_filter"] != "PASS" or row["alt_carried"] != "yes":
        return None
    clin = row["clin_sig"].lower()
    impact = row["impact"].upper()
    af = max_af_number(row["max_af"])
    pathogenic = "pathogenic" in clin
    benign = "benign" in clin
    if pathogenic and benign:
        return "T1_REVIEW", "conflicting ClinVar-style classifications"
    if pathogenic:
        return "T1", "pathogenic/likely pathogenic database classification"
    if impact == "HIGH":
        return "T1_REVIEW", "high-impact consequence; requires clinical database review"
    if any(term in clin for term in ("drug_response", "risk_factor", "protective", "association")):
        return "T2", "potentially actionable or phenotype-associated database annotation"
    if impact == "MODERATE" and (af is None or af <= 0.01):
        return "T2", "rare/unknown-frequency moderate-impact consequence"
    return None


def build_summary_row(
    fields: list[str],
    sample: dict[str, str],
    allele_number: int,
    alt: str,
    records: list[dict[str, str]],
    source: Path,
    layer: str,
) -> dict[str, str]:
    chrom, pos, original_id, ref, _, _, filt = fields[:7]
    best = best_csq(records)
    gt = sample.get("GT", ".")
    alts = fields[4].split(",")
    existing = [r.get("Existing_variation", "") for r in records]
    rsids = rsids_from(existing)
    consequences = [term for r in records for term in split_terms(r.get("Consequence", ""))]
    symbols = [r.get("SYMBOL", "") for r in records]
    clin = [term for r in records for term in split_terms(r.get("CLIN_SIG", ""))]
    phenotypes = [term for r in records for term in split_terms(r.get("PHENO", ""))]
    pubmed = [term for r in records for term in split_terms(r.get("PUBMED", ""))]
    return {
        "variant_key": f"GRCh37:{chrom.removeprefix('chr')}:{pos}:{ref}:{alt}",
        "source_layer": layer,
        "chrom": chrom,
        "position": pos,
        "ref": ref,
        "alt": alt,
        "allele_number": str(allele_number),
        "rsid": ";".join(rsids) if rsids else ".",
        "original_id": original_id,
        "vcf_filter": filt,
        "GT": gt,
        "genotype_bases": genotype_bases(ref, alts, gt),
        "alt_carried": "yes" if allele_number in called_alt_numbers(gt) else "no",
        "DP": sample.get("DP", "."),
        "GQ": sample.get("GQ", "."),
        "AD": alt_specific(sample.get("AD", "."), allele_number),
        "allele_fraction": alt_specific(sample.get("AF", "."), allele_number - 1),
        "impact": best.get("IMPACT", ".") or ".",
        "consequence": best.get("Consequence", ".") or ".",
        "all_consequences": unique_join(consequences),
        "symbol": best.get("SYMBOL", ".") or ".",
        "all_symbols": unique_join(symbols),
        "gene": best.get("Gene", ".") or ".",
        "feature_type": best.get("Feature_type", ".") or ".",
        "transcript": best.get("Feature", ".") or ".",
        "biotype": best.get("BIOTYPE", ".") or ".",
        "canonical": best.get("CANONICAL", ".") or ".",
        "mane_select": best.get("MANE_SELECT", ".") or ".",
        "mane_plus_clinical": best.get("MANE_PLUS_CLINICAL", ".") or ".",
        "hgvsc": unquote(best.get("HGVSc", ".") or "."),
        "hgvsp": unquote(best.get("HGVSp", ".") or "."),
        "protein_position": best.get("Protein_position", ".") or ".",
        "amino_acids": best.get("Amino_acids", ".") or ".",
        "codons": best.get("Codons", ".") or ".",
        "sift": best.get("SIFT", ".") or ".",
        "polyphen": best.get("PolyPhen", ".") or ".",
        "clin_sig": unique_join(clin),
        "phenotypes": unique_join(phenotypes),
        "max_af": best.get("MAX_AF", ".") or ".",
        "max_af_pops": best.get("MAX_AF_POPS", ".") or ".",
        "pubmed": unique_join(pubmed),
        "variant_class": best.get("VARIANT_CLASS", ".") or ".",
        "csq_count": str(len(records)),
        "source_file": source.name,
        "genome_build": "GRCh37/hg19",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Postprocess an Ensembl VEP VCF for Dante WGS")
    parser.add_argument("input", type=Path, help="VEP-annotated VCF or VCF.GZ")
    parser.add_argument("--rsid-vcf", type=Path, required=True, help="Output VCF.GZ with rsIDs copied into ID")
    parser.add_argument("--summary", type=Path, required=True, help="Compact allele-level TSV.GZ")
    parser.add_argument("--priority-output", type=Path, help="Optional legacy conservative candidate TSV")
    parser.add_argument("--append-priority", action="store_true", help="Append to an existing priority TSV")
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--layer", choices=["snp", "indel", "residual"], default="snp")
    parser.add_argument("--progress-every", type=int, default=500_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = args.input.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    required_paths = (args.rsid_vcf, args.summary, args.report)
    for path in required_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
    if args.priority_output:
        args.priority_output.parent.mkdir(parents=True, exist_ok=True)

    vcf_part, vcf_handle, vcf_raw_handle = open_atomic(args.rsid_vcf)
    summary_part, summary_handle, summary_raw_handle = open_atomic(args.summary)
    summary_writer = csv.DictWriter(summary_handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
    summary_writer.writeheader()

    priority_handle = None
    priority_writer = None
    if args.priority_output:
        priority_mode = "a" if args.append_priority and args.priority_output.exists() else "w"
        priority_handle = args.priority_output.open(priority_mode, encoding="utf-8", newline="")
        priority_writer = csv.DictWriter(
            priority_handle, fieldnames=PRIORITY_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        if priority_mode == "w" or args.priority_output.stat().st_size == 0:
            priority_writer.writeheader()

    csq_fields: list[str] = []
    counts: Counter[str] = Counter()
    impact_counts: Counter[str] = Counter()
    priority_counts: Counter[str] = Counter()

    try:
        with open_text(source) as input_handle:
            for raw in input_handle:
                if raw.startswith("##INFO=<ID=CSQ"):
                    match = CSQ_HEADER_RE.search(raw)
                    if not match:
                        raise RuntimeError("Could not parse the VEP CSQ Format header")
                    csq_fields = [x.strip() for x in match.group(1).strip().split("|")]
                    vcf_handle.write(raw)
                    continue
                if raw.startswith("#CHROM"):
                    if not csq_fields:
                        raise RuntimeError("No VEP CSQ header found; run VEP with --vcf first")
                    header_columns = raw.rstrip("\r\n").split("\t")
                    sample_count = max(0, len(header_columns) - 9)
                    if sample_count != 1:
                        raise RuntimeError(
                            f"Expected exactly one VCF sample, found {sample_count}; refusing to select one silently"
                        )
                    vcf_handle.write(
                        "##dante_postprocess=rsIDs copied from allele-matched VEP Existing_variation; full CSQ preserved\n"
                    )
                    vcf_handle.write(raw)
                    continue
                if raw.startswith("#"):
                    vcf_handle.write(raw)
                    continue
                if not raw.strip():
                    continue

                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 8:
                    raise RuntimeError(f"Malformed VCF row near record {counts['records'] + 1}")
                counts["records"] += 1
                info = parse_info(fields[7])
                csq_raw = info.get("CSQ", "")
                entries = [csq_record(x, csq_fields) for x in csq_raw.split(",") if x]
                by_allele: dict[int, list[dict[str, str]]] = defaultdict(list)
                for entry in entries:
                    token = entry.get("ALLELE_NUM", "")
                    if token.isdigit() and int(token) > 0:
                        by_allele[int(token)].append(entry)
                    elif len(fields[4].split(",")) == 1:
                        by_allele[1].append(entry)

                record_rsids = rsids_from(entry.get("Existing_variation", "") for entry in entries)
                original_id = fields[2]
                existing_ids = [] if original_id in {"", "."} else original_id.split(";")
                combined_ids = sorted(set(existing_ids + record_rsids))
                fields[2] = ";".join(combined_ids) if combined_ids else "."
                if record_rsids:
                    counts["records_with_rsid"] += 1
                    counts["rsids_total"] += len(record_rsids)
                else:
                    counts["records_without_rsid"] += 1
                vcf_handle.write("\t".join(fields) + "\n")

                sample = parse_sample(fields[8], fields[9]) if len(fields) > 9 else {}
                for allele_number, alt in enumerate(fields[4].split(","), start=1):
                    row = build_summary_row(
                        [fields[0], fields[1], original_id, fields[3], fields[4], fields[5], fields[6]],
                        sample,
                        allele_number,
                        alt,
                        by_allele.get(allele_number, []),
                        source,
                        args.layer,
                    )
                    summary_writer.writerow(row)
                    counts["allele_rows"] += 1
                    impact_counts[row["impact"]] += 1
                    priority = priority_for(row) if priority_writer else None
                    if priority and priority_writer:
                        tier, reason = priority
                        priority_writer.writerow(
                            {
                                "priority_tier": tier,
                                "priority_reason": reason,
                                **row,
                                "phenotype_or_question": "",
                            }
                        )
                        priority_counts[tier] += 1
                        counts["priority_rows"] += 1

                if args.progress_every and counts["records"] % args.progress_every == 0:
                    print(f"  {source.name}: {counts['records']:,} records", flush=True)

        finish_atomic(vcf_part, vcf_handle, vcf_raw_handle, args.rsid_vcf)
        finish_atomic(summary_part, summary_handle, summary_raw_handle, args.summary)
        if priority_handle:
            priority_handle.close()
    except Exception:
        for handle in (vcf_handle, summary_handle, priority_handle):
            if handle is None:
                continue
            try:
                handle.close()
            except Exception:
                pass
        for raw_handle in (vcf_raw_handle, summary_raw_handle):
            if raw_handle is None:
                continue
            try:
                raw_handle.close()
            except Exception:
                pass
        vcf_part.unlink(missing_ok=True)
        summary_part.unlink(missing_ok=True)
        raise

    pct = 100.0 * counts["records_with_rsid"] / max(1, counts["records"])
    lines = [
        "DANTE VEP POSTPROCESS REPORT",
        "=" * 40,
        f"Input: {source.name}",
        "Genome build: GRCh37/hg19",
        f"VCF records: {counts['records']:,}",
        f"Records with at least one rsID: {counts['records_with_rsid']:,} ({pct:.2f}%)",
        f"Records without rsID: {counts['records_without_rsid']:,}",
        f"Allele-level summary rows: {counts['allele_rows']:,}",
        f"Priority candidates emitted: {counts['priority_rows']:,}",
        "",
        "IMPACT COUNTS",
        *[f"{key}: {value:,}" for key, value in sorted(impact_counts.items())],
        "",
        "PRIORITY COUNTS",
        *[f"{key}: {value:,}" for key, value in sorted(priority_counts.items())],
        "",
        "NOTES",
        "- Full transcript-level VEP annotations remain in the source VEP VCF CSQ field.",
        "- The compact TSV selects one transcript deterministically and also preserves aggregate fields.",
        "- Priority tiers are triage rules, not clinical classifications.",
        "- Absence of an rsID does not make a variant false or clinically irrelevant.",
        "",
    ]
    args.report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Created {args.rsid_vcf.name}, {args.summary.name}, {args.report.name}")
    print(f"rsID coverage: {counts['records_with_rsid']:,}/{counts['records']:,} ({pct:.2f}%)")
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
