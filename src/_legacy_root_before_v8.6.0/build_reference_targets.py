#!/usr/bin/env python3
"""Build a deterministic GRCh37 rsID target panel from ClinVar and chip manifests."""

from __future__ import annotations

import argparse
import filecmp
import gzip
import os
import re
from dataclasses import dataclass
from pathlib import Path
from itertools import chain
from typing import Iterable, TextIO


RSID = re.compile(r"^rs(\d+)$", re.IGNORECASE)
CANONICAL = {str(value): value for value in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}
VALID_GENOTYPE = re.compile(r"^[ACGT]{1,2}$")


@dataclass
class Target:
    rsid: str
    chromosome: str
    position: int
    ref: str
    alts: str
    sources: set[str]
    reference_genotype: str = "--"
    reference_genotype_sources: str = ""
    reference_genotype_status: str = "MISSING"


def normalize_chromosome(value: str) -> str | None:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value == "M":
        value = "MT"
    return value if value in CANONICAL else None


def normalize_genotype(value: str) -> str | None:
    value = value.strip().upper().replace("/", "").replace("|", "")
    if not VALID_GENOTYPE.fullmatch(value):
        return None
    return "".join(sorted(value))


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("rt", encoding="utf-8", errors="replace", newline="")


def clinvar_targets(path: Path) -> Iterable[tuple[str, str, int, str, str, str, str | None]]:
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                continue
            chromosome = normalize_chromosome(fields[0])
            if chromosome is None or not fields[1].isdigit():
                continue
            ref = fields[3].upper()
            alts = [alt.upper() for alt in fields[4].split(",")]
            if len(ref) != 1 or not alts or any(len(alt) != 1 for alt in alts):
                continue
            info = fields[7]
            match = re.search(r"(?:^|;)RS=(\d+)(?:;|$)", info)
            if not match:
                continue
            yield f"rs{match.group(1)}", chromosome, int(fields[1]), ref, ",".join(alts), "ClinVar", None


def panel_targets(path: Path) -> Iterable[tuple[str, str, int, str, str, str, str | None]]:
    source_name = path.name
    with open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = re.split(r"[\t,; ]+", stripped)
            if len(fields) < 3 or not RSID.match(fields[0]):
                continue
            chromosome = normalize_chromosome(fields[1])
            if chromosome is None:
                continue
            try:
                position = int(fields[2])
            except ValueError:
                continue
            if position <= 0:
                continue
            genotype = normalize_genotype(fields[3]) if len(fields) >= 4 else None
            yield fields[0].lower(), chromosome, position, ".", ".", source_name, genotype


def append_source(existing: str, source: str) -> str:
    if not existing:
        return source
    values = existing.split(",")
    if source not in values:
        values.append(source)
    return ",".join(values)


def add_reference_genotype(
    target: Target,
    genotype: str | None,
    source: str,
    genotype_conflicts: dict[str, dict[str, set[str]]],
) -> None:
    if genotype is None:
        return
    if target.reference_genotype_status == "MISSING":
        target.reference_genotype = genotype
        target.reference_genotype_sources = source
        target.reference_genotype_status = "AVAILABLE"
        return
    if target.reference_genotype_status == "AVAILABLE":
        if target.reference_genotype == genotype:
            target.reference_genotype_sources = append_source(target.reference_genotype_sources, source)
            return
        details = genotype_conflicts.setdefault(target.rsid, {})
        details.setdefault(target.reference_genotype, set()).update(target.reference_genotype_sources.split(","))
        details.setdefault(genotype, set()).add(source)
        target.reference_genotype = "--"
        target.reference_genotype_sources = append_source(target.reference_genotype_sources, source)
        target.reference_genotype_status = "CONFLICT"
        return
    details = genotype_conflicts.setdefault(target.rsid, {})
    details.setdefault(genotype, set()).add(source)
    target.reference_genotype_sources = append_source(target.reference_genotype_sources, source)


def publish_text(path: Path, lines: Iterable[str]) -> None:
    """Atomically publish text, preserving the old file identity if content is unchanged."""
    part = Path(str(path) + ".part")
    path.parent.mkdir(parents=True, exist_ok=True)
    with part.open("wt", encoding="utf-8", newline="\n") as handle:
        for line in lines:
            handle.write(line)
    if path.exists() and filecmp.cmp(part, path, shallow=False):
        part.unlink()
        return
    os.replace(part, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clinvar", type=Path, required=True)
    parser.add_argument("--panel-dir", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--coordinate-conflicts", type=Path, required=True)
    parser.add_argument("--genotype-conflicts", type=Path, required=True)
    args = parser.parse_args()

    records: dict[str, Target] = {}
    coordinate_conflicts: dict[str, dict[tuple[str, int], set[str]]] = {}
    genotype_conflicts: dict[str, dict[str, set[str]]] = {}
    source_counts: dict[str, int] = {}
    genotype_counts: dict[str, int] = {}

    inputs: list[tuple[Path, Iterable[tuple[str, str, int, str, str, str, str | None]]]] = [
        (args.clinvar, clinvar_targets(args.clinvar))
    ]
    args.panel_dir.mkdir(parents=True, exist_ok=True)
    for path in sorted(args.panel_dir.iterdir()):
        if path.is_file() and re.search(r"(?i)\.(?:txt|tsv|csv)(?:\.gz)?$", path.name):
            inputs.append((path, panel_targets(path)))

    for _, iterator in inputs:
        for rsid, chromosome, position, ref, alts, source, genotype in iterator:
            source_counts[source] = source_counts.get(source, 0) + 1
            if genotype is not None:
                genotype_counts[source] = genotype_counts.get(source, 0) + 1
            existing = records.get(rsid)
            if existing is None:
                target = Target(rsid, chromosome, position, ref, alts, {source})
                add_reference_genotype(target, genotype, source, genotype_conflicts)
                records[rsid] = target
                continue

            if rsid in coordinate_conflicts:
                coordinate_conflicts[rsid].setdefault((chromosome, position), set()).add(source)
                if (existing.chromosome, existing.position) == (chromosome, position):
                    existing.sources.add(source)
                continue

            if (existing.chromosome, existing.position) != (chromosome, position):
                coordinate_conflicts[rsid] = {
                    (existing.chromosome, existing.position): set(existing.sources),
                    (chromosome, position): {source},
                }
                continue

            existing.sources.add(source)
            if existing.ref == "." and ref != ".":
                existing.ref = ref
                existing.alts = alts
            add_reference_genotype(existing, genotype, source, genotype_conflicts)

    for rsid in coordinate_conflicts:
        records.pop(rsid, None)
        genotype_conflicts.pop(rsid, None)

    ordered = sorted(records.values(), key=lambda row: (CANONICAL[row.chromosome], row.position, int(row.rsid[2:])))
    if not ordered:
        raise SystemExit("No valid GRCh37 rsID SNP targets were created.")

    publish_text(
        args.targets,
        chain(
            ["rsid\tchromosome\tposition\tref\talts\tsources\treference_genotype\treference_genotype_status\treference_genotype_sources\n"],
            (
                f"{row.rsid}\t{row.chromosome}\t{row.position}\t{row.ref}\t{row.alts}\t{','.join(sorted(row.sources))}\t"
                f"{row.reference_genotype}\t{row.reference_genotype_status}\t{row.reference_genotype_sources or '.'}\n"
                for row in ordered
            ),
        ),
    )

    unique_position_count = 0
    previous_position: tuple[str, int] | None = None

    def bed_lines() -> Iterable[str]:
        nonlocal unique_position_count, previous_position
        for row in ordered:
            current = (row.chromosome, row.position)
            if current == previous_position:
                continue
            previous_position = current
            unique_position_count += 1
            yield f"{row.chromosome}\t{row.position - 1}\t{row.position}\n"

    publish_text(args.bed, bed_lines())

    coordinate_lines = ["rsid\tchromosome\tposition\tsources\n"]
    for rsid in sorted(coordinate_conflicts, key=lambda item: int(item[2:])):
        locations = coordinate_conflicts[rsid]
        for (chromosome, position), sources in sorted(
            locations.items(), key=lambda item: (CANONICAL[item[0][0]], item[0][1])
        ):
            coordinate_lines.append(f"{rsid}\t{chromosome}\t{position}\t{','.join(sorted(sources))}\n")
    publish_text(args.coordinate_conflicts, coordinate_lines)

    genotype_lines = ["rsid\tchromosome\tposition\tgenotype\tsources\n"]
    for rsid in sorted(genotype_conflicts, key=lambda item: int(item[2:])):
        target = records.get(rsid)
        if target is None:
            continue
        for genotype, sources in sorted(genotype_conflicts[rsid].items()):
            genotype_lines.append(
                f"{rsid}\t{target.chromosome}\t{target.position}\t{genotype}\t{','.join(sorted(sources))}\n"
            )
    publish_text(args.genotype_conflicts, genotype_lines)

    available_genotypes = sum(1 for row in ordered if row.reference_genotype_status == "AVAILABLE")
    conflicting_genotypes = sum(1 for row in ordered if row.reference_genotype_status == "CONFLICT")
    report_lines = [
        "Dante reference target panel report\n",
        "Assembly: GRCh37\n",
        f"Unique rsIDs: {len(ordered)}\n",
        f"Unique positions: {unique_position_count}\n",
        f"Conflicting-coordinate rsIDs excluded: {len(coordinate_conflicts)}\n",
        f"Reference-source genotypes available for fallback: {available_genotypes}\n",
        f"Reference-source genotype conflicts disabled for fallback: {conflicting_genotypes}\n",
        f"Optional panel files: {max(0, len(inputs) - 1)}\n",
        "Input rows by source:\n",
    ]
    report_lines.extend(f"  {name}: {count}\n" for name, count in sorted(source_counts.items()))
    report_lines.append("Rows with usable genotype by source:\n")
    report_lines.extend(f"  {name}: {count}\n" for name, count in sorted(genotype_counts.items()))
    report_lines.extend(
        [
            f"Full coordinate-conflict table: {args.coordinate_conflicts.name}\n",
            f"Full reference-source genotype-conflict table: {args.genotype_conflicts.name}\n",
        ]
    )
    publish_text(args.report, report_lines)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
