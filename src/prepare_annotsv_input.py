#!/usr/bin/env python3
"""Prepare a high-confidence canonical-chromosome VCF for AnnotSV.

The original Dante CNV/SV VCF is never changed. For SV breakend pairs, an
entire EVENT/MATEID group is retained or rejected together so the preparation
step does not leave orphaned mates.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


BND_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


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


def canonical_chrom(value: str) -> bool:
    chrom = value.removeprefix("chr")
    return chrom in {*(str(x) for x in range(1, 23)), "X", "Y", "M", "MT"}


@dataclass
class Record:
    index: int
    raw: str
    fields: list[str]
    info: dict[str, str]
    sample: dict[str, str]

    @property
    def identifier(self) -> str:
        return self.fields[2]

    @property
    def event(self) -> str:
        return self.info.get("EVENT", "")

    @property
    def mates(self) -> list[str]:
        return [x for x in self.info.get("MATEID", "").split(",") if x and x != "."]

    def technical_pass(self) -> bool:
        if self.fields[6] != "PASS":
            return False
        ft = self.sample.get("FT", ".")
        return ft in {"", ".", "PASS"}

    def canonical(self) -> bool:
        if not canonical_chrom(self.fields[0]):
            return False
        match = BND_RE.search(self.fields[4])
        return not match or canonical_chrom(match.group(1))


def group_key(record: Record) -> str:
    if record.event:
        return "EVENT:" + record.event
    if record.identifier not in {"", "."} and record.mates:
        return "MATES:" + "|".join(sorted([record.identifier, *record.mates]))
    return f"ROW:{record.index}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare PASS canonical CNV/SV calls for AnnotSV")
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--kind", choices=["cnv", "sv"], required=True)
    parser.add_argument("--report", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source = args.input.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    headers: list[str] = []
    records: list[Record] = []
    with open_text(source) as handle:
        for raw in handle:
            if raw.startswith("#"):
                headers.append(raw)
                continue
            if not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise RuntimeError(f"Malformed VCF row {len(records) + 1}")
            info = parse_info(fields[7])
            sample = parse_sample(fields[8], fields[9]) if len(fields) > 9 else {}
            records.append(Record(len(records), raw, fields, info, sample))

    groups: dict[str, list[Record]] = defaultdict(list)
    if args.kind == "sv":
        for record in records:
            groups[group_key(record)].append(record)
    else:
        for record in records:
            groups[f"ROW:{record.index}"].append(record)

    counts: Counter[str] = Counter()
    selected_indices: set[int] = set()
    for members in groups.values():
        counts["groups_total"] += 1
        if not all(record.technical_pass() for record in members):
            counts["groups_rejected_filter"] += 1
            counts["rows_rejected_filter"] += len(members)
            continue
        if not all(record.canonical() for record in members):
            counts["groups_rejected_contig"] += 1
            counts["rows_rejected_contig"] += len(members)
            continue
        counts["groups_selected"] += 1
        counts["rows_selected"] += len(members)
        selected_indices.update(record.index for record in members)

    part = args.output.with_name(args.output.name + ".part")
    try:
        with part.open("wt", encoding="utf-8", newline="") as output:
            for raw in headers:
                if raw.startswith("#CHROM"):
                    output.write(
                        "##dante_annotsv_preparation=site FILTER PASS; sample FT PASS/absent; canonical chromosomes; event pairs retained together\n"
                    )
                output.write(raw)
            for record in records:
                if record.index in selected_indices:
                    output.write(record.raw if record.raw.endswith("\n") else record.raw + "\n")
        os.replace(part, args.output)
    except Exception:
        part.unlink(missing_ok=True)
        raise

    lines = [
        "DANTE ANNOTSV INPUT REPORT",
        "=" * 40,
        f"Input: {source.name}",
        f"Kind: {args.kind}",
        "Genome build: GRCh37/hg19",
        f"Input rows: {len(records):,}",
        f"Logical groups: {counts['groups_total']:,}",
        f"Selected groups: {counts['groups_selected']:,}",
        f"Selected rows: {counts['rows_selected']:,}",
        f"Rejected by site/sample filter: {counts['rows_rejected_filter']:,} rows",
        f"Rejected noncanonical/decoy: {counts['rows_rejected_contig']:,} rows",
        "",
        "The original VCF was not modified. AnnotSV receives only high-confidence canonical calls.",
        "For SVs, EVENT/MATEID-linked rows are retained or rejected together.",
        "",
    ]
    args.report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Created {args.output.name}: {counts['rows_selected']:,}/{len(records):,} rows")
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
