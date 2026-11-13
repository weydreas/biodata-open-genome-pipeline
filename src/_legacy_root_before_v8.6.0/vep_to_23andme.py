#!/usr/bin/env python3
"""Create a strict 23andMe-style genotype file from a minimal VEP VCF.

The input is expected to be the Dante filtered SNP VCF after Ensembl VEP was
run with --check_existing and --allele_number.  Only high-confidence canonical
chromosome SNP calls with a real allele-matched rsID are written.

The output has exactly the four data columns used by 23andMe raw files:

    rsid    chromosome    position    genotype

This is intentionally a *variant-only* export.  It cannot reconstruct
homozygous-reference positions that were absent from Dante's VCF.
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import sys
import traceback
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


RSID_RE = re.compile(r"(?<![A-Za-z0-9])rs\d+", re.IGNORECASE)
CSQ_HEADER_RE = re.compile(r"Format:\s*([^\"]+)")
CANONICAL = {*(str(i) for i in range(1, 23)), "X", "Y", "M", "MT"}
BASES = {"A", "C", "G", "T"}


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


@contextmanager
def open_anonymous_gzip_text(path: Path) -> Iterator[TextIO]:
    """Write gzip data without embedding the temporary .part filename."""
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gzip_handle:
            with io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="") as text_handle:
                yield text_handle


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


def normalise_chromosome(value: str) -> str | None:
    chrom = value[3:] if value.lower().startswith("chr") else value
    if chrom not in CANONICAL:
        return None
    return "MT" if chrom == "M" else chrom


def genotype_indices(gt: str) -> list[int] | None:
    if not gt or gt in {".", "./.", ".|."}:
        return None
    tokens = re.split(r"[/|]", gt)
    if not tokens or any(not token.isdigit() for token in tokens):
        return None
    return [int(token) for token in tokens]


def genotype_bases(ref: str, alts: list[str], indices: list[int]) -> str | None:
    alleles = [ref, *alts]
    if any(index >= len(alleles) for index in indices):
        return None
    called = [alleles[index].upper() for index in indices]
    if any(base not in BASES for base in called):
        return None
    return "".join(called)


def parse_csq_fields(header_line: str) -> list[str]:
    match = CSQ_HEADER_RE.search(header_line)
    if not match:
        raise RuntimeError("Could not parse the VEP CSQ Format header")
    fields = [field.strip() for field in match.group(1).strip().split("|")]
    required = {"Allele", "Existing_variation", "ALLELE_NUM"}
    missing = sorted(required.difference(fields))
    if missing:
        raise RuntimeError(
            "VEP CSQ is missing required field(s): " + ", ".join(missing)
            + ". Run the supplied minimal VEP command."
        )
    return fields


def rsids_by_alt(info: dict[str, str], csq_fields: list[str], alts: list[str]) -> dict[int, set[str]]:
    result: dict[int, set[str]] = {}
    raw_csq = info.get("CSQ", "")
    if not raw_csq:
        return result
    allele_index = csq_fields.index("Allele")
    existing_index = csq_fields.index("Existing_variation")
    number_index = csq_fields.index("ALLELE_NUM")
    for raw_record in raw_csq.split(","):
        values = raw_record.split("|")
        if len(values) < len(csq_fields):
            values.extend([""] * (len(csq_fields) - len(values)))
        number_text = values[number_index]
        alt_number: int | None = int(number_text) if number_text.isdigit() else None
        if alt_number is None:
            # Fallback is safe for SNPs because VEP's allele is a single base.
            allele = values[allele_index].upper()
            matches = [i + 1 for i, alt in enumerate(alts) if alt.upper() == allele]
            if len(matches) == 1:
                alt_number = matches[0]
        if alt_number is None or alt_number < 1 or alt_number > len(alts):
            continue
        found = {match.lower() for match in RSID_RE.findall(values[existing_index] or "")}
        if found:
            result.setdefault(alt_number, set()).update(found)
    return result


def rsid_sort_key(value: str) -> int:
    return int(value[2:])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create Dante_23andMe_style.txt.gz from a minimal VEP VCF")
    parser.add_argument("input", type=Path, help="VEP-annotated SNP VCF or VCF.GZ")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--progress-every", type=int, default=500_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    report = args.report.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    output.parent.mkdir(parents=True, exist_ok=True)
    report.parent.mkdir(parents=True, exist_ok=True)

    part = output.with_name(output.name + ".part")
    counts: Counter[str] = Counter()
    csq_fields: list[str] = []
    sample_name = ""

    try:
        with open_text(source) as input_handle, open_anonymous_gzip_text(part) as out:
            out.write("# Dante WGS variant-only export; GRCh37 forward-strand alleles\n")
            out.write("# Only FILTER=PASS canonical-chromosome SNPs with allele-matched rsIDs are included\n")
            out.write("# rsid\tchromosome\tposition\tgenotype\n")

            for raw in input_handle:
                if raw.startswith("##INFO=<ID=CSQ"):
                    csq_fields = parse_csq_fields(raw)
                    continue
                if raw.startswith("#CHROM"):
                    columns = raw.rstrip("\r\n").split("\t")
                    if len(columns) != 10:
                        raise RuntimeError(
                            f"Expected one Dante sample column, found {max(0, len(columns) - 9)}"
                        )
                    sample_name = columns[9]
                    continue
                if raw.startswith("#") or not raw.strip():
                    continue
                if not csq_fields:
                    raise RuntimeError("No VEP CSQ header found")

                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 10:
                    counts["malformed_or_wrong_sample_count"] += 1
                    continue
                counts["input_records"] += 1
                chrom_raw, pos, _, ref, alt_field, _, filt, info_text, fmt, sample_text = fields

                if filt != "PASS":
                    counts["skipped_non_pass"] += 1
                    continue
                chrom = normalise_chromosome(chrom_raw)
                if chrom is None:
                    counts["skipped_noncanonical_chromosome"] += 1
                    continue

                ref = ref.upper()
                alts = [alt.upper() for alt in alt_field.split(",")]
                if ref not in BASES or not alts or any(alt not in BASES for alt in alts):
                    counts["skipped_non_snp"] += 1
                    continue

                sample = parse_sample(fmt, sample_text)
                indices = genotype_indices(sample.get("GT", "."))
                if indices is None:
                    counts["skipped_missing_genotype"] += 1
                    continue
                genotype = genotype_bases(ref, alts, indices)
                if genotype is None:
                    counts["skipped_invalid_genotype"] += 1
                    continue

                called_alt_numbers = {index for index in indices if index > 0}
                if not called_alt_numbers:
                    counts["skipped_homozygous_reference"] += 1
                    continue

                matches = rsids_by_alt(parse_info(info_text), csq_fields, alts)
                row_rsids: set[str] = set()
                for alt_number in called_alt_numbers:
                    row_rsids.update(matches.get(alt_number, set()))
                if not row_rsids:
                    counts["skipped_without_allele_matched_rsid"] += 1
                    continue

                for rsid in sorted(row_rsids, key=rsid_sort_key):
                    out.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")
                    counts["output_rsid_rows"] += 1
                counts["output_variant_records"] += 1
                if len(alts) > 1:
                    counts["output_multiallelic_records"] += 1

                if args.progress_every and counts["input_records"] % args.progress_every == 0:
                    print(
                        f"  Read {counts['input_records']:,} SNP records; "
                        f"wrote {counts['output_rsid_rows']:,} rsID rows",
                        flush=True,
                    )

        os.replace(part, output)
    except Exception:
        part.unlink(missing_ok=True)
        raise

    if counts["output_rsid_rows"] == 0:
        output.unlink(missing_ok=True)
        raise RuntimeError("No allele-matched rsID rows were produced")

    report_lines = [
        "DANTE 23ANDME-STYLE EXPORT REPORT",
        "=" * 44,
        f"Input: {source.name}",
        f"Sample: {sample_name or '.'}",
        "Genome build: GRCh37/hg19",
        "Output is variant-only; absent positions are not assumed to be 0/0.",
        "Output alleles use the forward/reference orientation of the Dante VCF.",
        "",
    ]
    for key in (
        "input_records",
        "output_variant_records",
        "output_rsid_rows",
        "output_multiallelic_records",
        "skipped_non_pass",
        "skipped_noncanonical_chromosome",
        "skipped_non_snp",
        "skipped_missing_genotype",
        "skipped_invalid_genotype",
        "skipped_homozygous_reference",
        "skipped_without_allele_matched_rsid",
        "malformed_or_wrong_sample_count",
    ):
        report_lines.append(f"{key}: {counts[key]:,}")
    report_lines.append("")
    report.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"Created {output.name}: {counts['output_rsid_rows']:,} rsID rows")
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
