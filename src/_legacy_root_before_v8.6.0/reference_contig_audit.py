#!/usr/bin/env python3
"""Audit VCF header contigs against the installed reference target list."""

from __future__ import annotations

import argparse
import gzip
import os
import re
from pathlib import Path
from typing import TextIO


CONTIG_LINE = re.compile(r"^##contig=<(.*)>$")
CONTIG_ID = re.compile(r"(?:^|,)ID=([^,>]+)")


def open_vcf_text(path: Path) -> TextIO:
    with path.open("rb") as raw:
        is_gzip = raw.read(2) == b"\x1f\x8b"
    if is_gzip:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("rt", encoding="utf-8", errors="replace", newline="")


def load_reference_contigs(path: Path) -> tuple[list[str], set[str]]:
    ordered: list[str] = []
    seen: set[str] = set()
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            contig = line.split("\t", 1)[0].split(None, 1)[0]
            if contig and contig not in seen:
                seen.add(contig)
                ordered.append(contig)
    if not ordered:
        raise ValueError(f"Reference target list is empty: {path}")
    return ordered, seen


def load_vcf_header_contigs(path: Path) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    found_column_header = False
    with open_vcf_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("#CHROM\t") or line == "#CHROM":
                found_column_header = True
                break
            match = CONTIG_LINE.match(line)
            if not match:
                continue
            id_match = CONTIG_ID.search(match.group(1))
            if not id_match:
                continue
            contig = id_match.group(1)
            if contig not in seen:
                seen.add(contig)
                ordered.append(contig)
    if not found_column_header:
        raise ValueError(f"VCF header terminator #CHROM was not found: {path}")
    return ordered


def parse_vcf_spec(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--vcf must have the form LABEL=PATH")
    label, raw_path = value.split("=", 1)
    label = label.strip().upper()
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("--vcf must have the form LABEL=PATH")
    return label, Path(raw_path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--vcf", action="append", type=parse_vcf_spec, required=True)
    args = parser.parse_args()

    reference_order, reference = load_reference_contigs(args.targets)
    report_lines = [
        "Dante reference-contig audit",
        "Status: complete",
        f"Reference target file: {args.targets.name}",
        f"Reference FASTA contigs: {len(reference_order)}",
        (
            "Policy: records on contigs absent from the installed human GRCh37 "
            "FASTA are excluded before normalization."
        ),
        "",
    ]

    print(f"Reference contig audit: {len(reference_order)} sequences in the installed FASTA")
    for label, path in args.vcf:
        if not path.is_file():
            raise FileNotFoundError(path)
        header_contigs = load_vcf_header_contigs(path)
        absent = [contig for contig in header_contigs if contig not in reference]
        present_count = len(header_contigs) - len(absent)
        report_lines.extend(
            [
                f"[{label}]",
                f"Input: {path.name}",
                f"VCF header contigs: {len(header_contigs)}",
                f"Header contigs present in FASTA: {present_count}",
                f"Header contigs absent from FASTA: {len(absent)}",
                "Excluded contig definitions (records are excluded if present):",
            ]
        )
        report_lines.extend(absent or ["(none)"])
        report_lines.append("")

        if absent:
            preview = ", ".join(absent[:12])
            suffix = " ..." if len(absent) > 12 else ""
            print(
                f"{label}: {len(absent)} VCF header contig(s) absent from FASTA; "
                f"excluded if present: {preview}{suffix}"
            )
        else:
            print(f"{label}: all declared VCF contigs are present in the FASTA")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    part = args.report.with_name(args.report.name + ".part")
    with part.open("wt", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(report_lines).rstrip() + "\n")
    os.replace(part, args.report)
    print(f"Reference contig audit report: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
