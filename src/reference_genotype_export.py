#!/usr/bin/env python3
"""Create a tiered reference-complete 23andMe-style export with provenance."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import re
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


CHROMOSOME_RANK = {str(value): value for value in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}


@dataclass(frozen=True)
class RawCall:
    chromosome: str
    position: int
    ref: str
    alts: tuple[str, ...]
    sample: dict[str, str]
    info: dict[str, str]


@dataclass(frozen=True)
class CallEvaluation:
    genotype: str
    tier: str
    reason: str
    dp: int | None = None
    quality: int | None = None
    quality_metric: str = "."
    ad: str = "."
    allele_fractions: str = "."


@dataclass(frozen=True)
class ExportRow:
    rsid: str
    chromosome: str
    position: int
    genotype: str
    reason: str
    reference_genotype: str = "--"
    reference_status: str = "MISSING"
    reference_sources: str = "."
    fastq_candidate_genotype: str = "--"
    fastq_tier: str = "NONE"
    fastq_dp: int | None = None
    fastq_quality: int | None = None
    fastq_quality_metric: str = "."
    fastq_ad: str = "."
    fastq_allele_fractions: str = "."

    @property
    def key(self) -> tuple[int, int, int]:
        return CHROMOSOME_RANK[self.chromosome], self.position, int(self.rsid[2:])


def normalize_chromosome(value: str) -> str | None:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value == "M":
        value = "MT"
    return value if value in CHROMOSOME_RANK else None


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="")
    return path.open("rt", encoding="utf-8", errors="strict", newline="")


def parse_info(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in value.split(";"):
        if not field:
            continue
        key, separator, item = field.partition("=")
        result[key] = item if separator else "1"
    return result


def raw_call_rows(path: Path) -> Iterator[RawCall]:
    """Stream call VCF rows in canonical coordinate order without retaining the file in RAM."""
    previous: tuple[int, int] | None = None
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                continue
            chromosome = normalize_chromosome(fields[0])
            if chromosome is None:
                continue
            try:
                position = int(fields[1])
            except ValueError:
                continue
            key = (CHROMOSOME_RANK[chromosome], position)
            if previous is not None and key < previous:
                raise ValueError("Reference-call VCF is not in canonical coordinate order")
            previous = key
            ref = fields[3].upper()
            alts = tuple(alt.upper() for alt in fields[4].split(",") if alt and alt != ".")
            format_names = fields[8].split(":")
            format_values = fields[9].split(":")
            sample = dict(zip(format_names, format_values))
            yield RawCall(chromosome, position, ref, alts, sample, parse_info(fields[7]))


class CallStream:
    """Coordinate-merge accessor for a sorted call VCF; memory use stays O(rows at one position)."""

    def __init__(self, path: Path) -> None:
        self._iterator = iter(raw_call_rows(path))
        self._current = next(self._iterator, None)

    @staticmethod
    def _key(call: RawCall) -> tuple[int, int]:
        return CHROMOSOME_RANK[call.chromosome], call.position

    def take(self, chromosome: str, position: int) -> list[RawCall]:
        requested = (CHROMOSOME_RANK[chromosome], position)
        while self._current is not None and self._key(self._current) < requested:
            self._current = next(self._iterator, None)
        if self._current is None or self._key(self._current) != requested:
            return []
        result: list[RawCall] = []
        while self._current is not None and self._key(self._current) == requested:
            result.append(self._current)
            self._current = next(self._iterator, None)
        return result


def count_target_rows(path: Path) -> int:
    """Count target data rows with bounded memory for truthful merge progress."""
    count = 0
    with path.open("rb") as handle:
        for _ in handle:
            count += 1
    return max(0, count - 1)


class MergeProgress:
    def __init__(self, stage: int, total_targets: int) -> None:
        self.stage = stage
        self.total = max(1, total_targets)
        self.processed = 0
        self.started = time.monotonic()
        self.last_emit = 0.0

    def emit(self, percent: float, detail: str) -> None:
        bounded = min(100.0, max(0.0, percent))
        print(f"[OVERALL] Stage {self.stage} | {bounded:.2f}% | {detail}", file=sys.stderr, flush=True)

    def advance(self) -> None:
        self.processed += 1
        now = time.monotonic()
        if self.processed < self.total and now - self.last_emit < 5.0:
            return
        self.last_emit = now
        ratio = min(1.0, self.processed / self.total)
        elapsed = max(0.001, now - self.started)
        speed = self.processed / elapsed
        remaining = max(0, self.total - self.processed)
        eta = remaining / speed if speed > 0 else 0.0
        hours, rem = divmod(int(eta), 3600)
        minutes, seconds = divmod(rem, 60)
        self.emit(
            2.0 + 93.0 * ratio,
            f"MERGE {ratio * 100:.2f}% | {self.processed:,} / {self.total:,} target rows | "
            f"{speed:,.0f} rows/s | ETA {hours:02d}:{minutes:02d}:{seconds:02d}",
        )


def integer(value: str | None) -> int | None:
    if value is None or value in {"", "."}:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def pl_reference_quality(pl_value: str | None, allele_count: int, ploidy: int) -> int | None:
    """PL-derived quality used only for homozygous-reference calls.

    bcftools does not emit FORMAT/GQ on non-variant records in the multiallelic
    model, even when constrained ALT alleles are retained with -A. The retained
    PL vector still provides the evidence needed to compare 0/0 (or haploid 0)
    with the best non-reference genotype.
    """
    if not pl_value or pl_value in {"", "."} or allele_count < 2 or ploidy not in {1, 2}:
        return None
    try:
        values = [int(x) for x in pl_value.split(",")]
    except ValueError:
        return None
    if not values or values[0] < 0:
        return None

    if ploidy == 2:
        expected = allele_count * (allele_count + 1) // 2
        if len(values) < expected:
            return None
        competitors = values[1:expected]
    else:
        if len(values) == allele_count:
            competitors = values[1:]
        else:
            expected = allele_count * (allele_count + 1) // 2
            if len(values) < expected:
                return None
            competitors = []
            for allele in range(1, allele_count):
                idx = allele * (allele + 1) // 2 + allele
                if idx >= len(values):
                    return None
                competitors.append(values[idx])

    if not competitors:
        return None
    return max(0, min(competitors) - values[0])


def tier_rank(value: str) -> int:
    return {"NONE": 0, "C": 1, "B": 2, "A": 3}.get(value, 0)


def evaluate_call(
    call: RawCall,
    target_ref: str,
    target_alts: str,
    tier_a_min_dp: int = 10,
    tier_a_min_quality: int = 30,
    tier_a_min_allele_depth: int = 3,
    tier_a_hom_fraction: float = 0.90,
    tier_a_het_min_fraction: float = 0.20,
    tier_a_het_max_fraction: float = 0.80,
    tier_b_min_dp: int = 5,
    tier_b_min_quality: int = 10,
    tier_b_hom_min_support: int = 4,
    tier_b_hom_fraction: float = 0.80,
    tier_b_het_min_allele_depth: int = 2,
    tier_b_het_min_fraction: float = 0.20,
    tier_b_het_max_fraction: float = 0.80,
    tier_c_min_dp: int = 5,
    tier_c_min_quality: int = 5,
) -> CallEvaluation:
    if target_ref != "." and call.ref != target_ref:
        return CallEvaluation("--", "NONE", "REF_MISMATCH")
    gt = call.sample.get("GT", ".").replace("|", "/")
    if gt in {"", ".", "./."} or "." in gt:
        return CallEvaluation("--", "NONE", "MISSING_GT")
    try:
        allele_indexes = [int(item) for item in gt.split("/")]
    except ValueError:
        return CallEvaluation("--", "NONE", "INVALID_GT")
    if len(allele_indexes) not in {1, 2}:
        return CallEvaluation("--", "NONE", "UNSUPPORTED_PLOIDY")
    alleles = (call.ref,) + call.alts
    if any(index < 0 or index >= len(alleles) for index in allele_indexes):
        return CallEvaluation("--", "NONE", "ALLELE_INDEX_OUT_OF_RANGE")
    if any(len(alleles[index]) != 1 or alleles[index] not in {"A", "C", "G", "T"} for index in allele_indexes):
        return CallEvaluation("--", "NONE", "NON_SNP_ALLELE")
    allowed_alts = {
        alt.strip().upper()
        for alt in (target_alts or ".").split(",")
        if alt.strip() and alt.strip() != "."
    }
    if allowed_alts and any(alleles[index] != call.ref and alleles[index] not in allowed_alts for index in allele_indexes):
        return CallEvaluation("--", "NONE", "ALT_NOT_IN_TARGET")

    genotype = "".join(alleles[index] for index in allele_indexes)
    dp = integer(call.sample.get("DP"))
    if dp is None:
        dp = integer(call.info.get("DP"))
    if dp is None:
        return CallEvaluation(genotype, "NONE", "MISSING_DP")

    quality = integer(call.sample.get("GQ"))
    quality_metric = "GQ" if quality is not None else "."
    if quality is None and all(index == 0 for index in allele_indexes):
        quality = pl_reference_quality(call.sample.get("PL"), len(alleles), len(allele_indexes))
        if quality is not None:
            quality_metric = "REFQ_PL"

    ad_value = call.sample.get("AD", ".")
    try:
        ad = [int(item) for item in ad_value.split(",")]
    except ValueError:
        ad = []
    if not ad or any(index >= len(ad) for index in allele_indexes):
        return CallEvaluation(genotype, "NONE", "MISSING_AD", dp, quality, quality_metric, ad_value, ".")
    depth_sum = sum(ad)
    if depth_sum <= 0:
        return CallEvaluation(genotype, "NONE", "MISSING_AD", dp, quality, quality_metric, ad_value, ".")

    distinct = sorted(set(allele_indexes))
    called_depths = [ad[index] for index in distinct]
    fractions = [ad[index] / depth_sum for index in distinct]
    fraction_text = ",".join(f"{value:.3f}" for value in fractions)
    is_hom = len(distinct) == 1

    if quality is None:
        return CallEvaluation(genotype, "NONE", "MISSING_QUALITY", dp, None, quality_metric, ad_value, fraction_text)

    # Tier A — HIGH. This is the old strict contract, now named explicitly.
    if dp >= tier_a_min_dp and quality >= tier_a_min_quality:
        if is_hom:
            if called_depths[0] >= tier_a_min_allele_depth and fractions[0] >= tier_a_hom_fraction:
                return CallEvaluation(genotype, "A", "TIER_A_HIGH", dp, quality, quality_metric, ad_value, fraction_text)
        else:
            if (
                all(value >= tier_a_min_allele_depth for value in called_depths)
                and all(tier_a_het_min_fraction <= value <= tier_a_het_max_fraction for value in fractions)
            ):
                return CallEvaluation(genotype, "A", "TIER_A_HIGH", dp, quality, quality_metric, ad_value, fraction_text)

    # Tier B — SUPPORTIVE. Deliberately admits DP 5-9 when the likelihood and
    # observed allele support are coherent enough for hypothesis-driven clinical work.
    if dp >= tier_b_min_dp and quality >= tier_b_min_quality:
        if is_hom:
            if called_depths[0] >= tier_b_hom_min_support and fractions[0] >= tier_b_hom_fraction:
                return CallEvaluation(genotype, "B", "TIER_B_SUPPORTIVE", dp, quality, quality_metric, ad_value, fraction_text)
        else:
            if (
                all(value >= tier_b_het_min_allele_depth for value in called_depths)
                and all(tier_b_het_min_fraction <= value <= tier_b_het_max_fraction for value in fractions)
            ):
                return CallEvaluation(genotype, "B", "TIER_B_SUPPORTIVE", dp, quality, quality_metric, ad_value, fraction_text)

    # Tier C — HYPOTHESIS. Retain for research review, but never promote to the
    # reference-complete genotype export. Require at least one observed read for
    # every called allele so pure likelihood artefacts are not surfaced as candidates.
    if dp >= tier_c_min_dp and quality >= tier_c_min_quality and all(value >= 1 for value in called_depths):
        return CallEvaluation(genotype, "C", "TIER_C_HYPOTHESIS", dp, quality, quality_metric, ad_value, fraction_text)

    reason = "LOW_DP" if dp < tier_c_min_dp else "LOW_QUALITY" if quality < tier_c_min_quality else "LOW_ALLELE_SUPPORT"
    return CallEvaluation(genotype, "NONE", reason, dp, quality, quality_metric, ad_value, fraction_text)


def target_rows(
    targets: Path,
    calls: CallStream | None,
    args: argparse.Namespace,
    progress: MergeProgress | None = None,
) -> Iterator[ExportRow]:
    with targets.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "rsid", "chromosome", "position", "ref", "alts", "reference_genotype",
            "reference_genotype_status", "reference_genotype_sources",
        }
        if not required.issubset(set(reader.fieldnames or [])):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"Target panel is missing required columns: {', '.join(missing)}")
        for item in reader:
            if progress is not None:
                progress.advance()
            chromosome = normalize_chromosome(item["chromosome"])
            if chromosome is None:
                continue
            position = int(item["position"])
            reference_genotype = item.get("reference_genotype", "--").upper()
            reference_status = item.get("reference_genotype_status", "MISSING").upper()
            reference_sources = item.get("reference_genotype_sources", ".") or "."
            possible = calls.take(chromosome, position) if calls is not None else []
            if not possible:
                yield ExportRow(
                    item["rsid"].lower(), chromosome, position, "--", "NO_COVERAGE_OR_RECORD",
                    reference_genotype, reference_status, reference_sources,
                )
                continue
            evaluated = [
                evaluate_call(
                    call, item["ref"].upper(), item.get("alts", "."),
                    args.tier_a_min_dp, args.tier_a_min_quality, args.tier_a_min_allele_depth,
                    args.tier_a_hom_fraction, args.tier_a_het_min_fraction, args.tier_a_het_max_fraction,
                    args.tier_b_min_dp, args.tier_b_min_quality, args.tier_b_hom_min_support,
                    args.tier_b_hom_fraction, args.tier_b_het_min_allele_depth,
                    args.tier_b_het_min_fraction, args.tier_b_het_max_fraction,
                    args.tier_c_min_dp, args.tier_c_min_quality,
                )
                for call in possible
            ]
            candidates = [value for value in evaluated if value.tier in {"A", "B", "C"}]
            if candidates:
                genotypes = {value.genotype for value in candidates}
                if len(genotypes) == 1:
                    best = max(
                        candidates,
                        key=lambda value: (tier_rank(value.tier), value.quality or -1, value.dp or -1),
                    )
                    export_genotype = best.genotype if tier_rank(best.tier) >= tier_rank(args.export_min_tier) else "--"
                    yield ExportRow(
                        item["rsid"].lower(), chromosome, position, export_genotype, best.reason,
                        reference_genotype, reference_status, reference_sources,
                        best.genotype, best.tier, best.dp, best.quality, best.quality_metric, best.ad, best.allele_fractions,
                    )
                else:
                    yield ExportRow(
                        item["rsid"].lower(), chromosome, position, "--", "AMBIGUOUS_CALLS",
                        reference_genotype, reference_status, reference_sources,
                    )
            else:
                best = max(evaluated, key=lambda value: (value.quality or -1, value.dp or -1))
                yield ExportRow(
                    item["rsid"].lower(), chromosome, position, "--", best.reason,
                    reference_genotype, reference_status, reference_sources,
                    best.genotype, best.tier, best.dp, best.quality, best.quality_metric, best.ad, best.allele_fractions,
                )


def base_rows(path: Path) -> Iterator[ExportRow]:
    previous: tuple[int, int, int] | None = None
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 4 or not re.fullmatch(r"(?i)rs\d+", fields[0]):
                continue
            chromosome = normalize_chromosome(fields[1])
            if chromosome is None:
                continue
            row = ExportRow(fields[0].lower(), chromosome, int(fields[2]), fields[3].upper(), "BASE_VARIANT")
            if previous is not None and row.key < previous:
                raise ValueError("Base 23andMe-style file is not in canonical coordinate/rsID order")
            previous = row.key
            yield row


def canonical_genotype(value: str) -> str:
    if value in {"--", ".", ""}:
        return value
    return "".join(sorted(value))


def disagree(first: str, second: str) -> bool:
    return (
        first not in {"--", ".", ""}
        and second not in {"--", ".", ""}
        and canonical_genotype(first) != canonical_genotype(second)
    )


def choose_target(row: ExportRow, fastq_source: str) -> tuple[str, str]:
    if row.genotype != "--":
        return row.genotype, fastq_source
    if row.reference_status == "AVAILABLE" and row.reference_genotype not in {"--", ".", ""}:
        return row.reference_genotype, row.reference_sources
    return "--", "NO_CALL"


def atomic_replace_with_retry(source: Path, destination: Path, *, timeout_seconds: float = 15.0) -> None:
    """Atomically publish a completed artifact, tolerating short-lived Windows file locks."""
    deadline = time.monotonic() + timeout_seconds
    delay = 0.10
    while True:
        try:
            os.replace(source, destination)
            return
        except PermissionError as exc:
            # Windows Defender/indexers and other readers can hold a file briefly.
            # Retry WinError 5/32; persistent locks still fail with an actionable error.
            winerror = getattr(exc, "winerror", None)
            if os.name != "nt" or winerror not in {5, 32} or time.monotonic() >= deadline:
                raise
            time.sleep(delay)
            delay = min(1.0, delay * 1.7)


def atomic_text(path: Path, text: str) -> None:
    part = Path(str(path) + ".part")
    path.parent.mkdir(parents=True, exist_ok=True)
    part.write_text(text, encoding="utf-8", newline="\n")
    atomic_replace_with_retry(part, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--calls", type=Path)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--discrepancies", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--fastq-source", choices=["NONE", "R1", "R1+R2"], default="NONE")
    parser.add_argument("--evidence", type=Path)
    parser.add_argument("--export-min-tier", choices=["A", "B"], default="B")
    parser.add_argument("--tier-a-min-dp", type=int, default=10)
    parser.add_argument("--tier-a-min-quality", type=int, default=30)
    parser.add_argument("--tier-a-min-allele-depth", type=int, default=3)
    parser.add_argument("--tier-a-hom-fraction", type=float, default=0.90)
    parser.add_argument("--tier-a-het-min-fraction", type=float, default=0.20)
    parser.add_argument("--tier-a-het-max-fraction", type=float, default=0.80)
    parser.add_argument("--tier-b-min-dp", type=int, default=5)
    parser.add_argument("--tier-b-min-quality", type=int, default=10)
    parser.add_argument("--tier-b-hom-min-support", type=int, default=4)
    parser.add_argument("--tier-b-hom-fraction", type=float, default=0.80)
    parser.add_argument("--tier-b-het-min-allele-depth", type=int, default=2)
    parser.add_argument("--tier-b-het-min-fraction", type=float, default=0.20)
    parser.add_argument("--tier-b-het-max-fraction", type=float, default=0.80)
    parser.add_argument("--tier-c-min-dp", type=int, default=5)
    parser.add_argument("--tier-c-min-quality", type=int, default=5)
    args = parser.parse_args()

    stage_number = 17 if args.calls is not None else 14
    total_targets = count_target_rows(args.targets)
    progress = MergeProgress(stage_number, total_targets)
    progress.emit(1.0, f"prepare streaming merge | {total_targets:,} target rows | FASTQ source {args.fastq_source}")
    calls = CallStream(args.calls) if args.calls is not None else None
    targets = iter(target_rows(args.targets, calls, args, progress))
    base = iter(base_rows(args.base))
    target = next(targets, None)
    existing = next(base, None)
    counters: Counter[str] = Counter()
    reasons: Counter[str] = Counter()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.discrepancies.parent.mkdir(parents=True, exist_ok=True)
    if args.evidence is not None:
        args.evidence.parent.mkdir(parents=True, exist_ok=True)
    output_part = Path(str(args.output) + ".part")
    provenance_part = Path(str(args.provenance) + ".part")
    discrepancy_part = Path(str(args.discrepancies) + ".part")
    evidence_part = Path(str(args.evidence) + ".part") if args.evidence is not None else None

    evidence_raw = evidence_part.open("wb") if evidence_part is not None else None
    evidence_gzip = None
    evidence_text = None
    try:
        evidence_gzip = gzip.GzipFile(filename="", mode="wb", fileobj=evidence_raw, mtime=0) if evidence_raw is not None else None
        evidence_text = io.TextIOWrapper(evidence_gzip, encoding="utf-8", newline="\n") if evidence_gzip is not None else None
        if evidence_text is not None:
            evidence_text.write(
                "rsid\tchromosome\tposition\tgenotype\ttier\tdp\tquality\tquality_metric\tAD\tallele_fractions\tselected_in_export\n"
            )

        with provenance_part.open("wt", encoding="utf-8", newline="\n") as provenance, discrepancy_part.open(
            "wt", encoding="utf-8", newline="\n"
        ) as discrepancies:
            provenance.write(
                "rsid\tchromosome\tposition\tgenotype\tselected_source\tfastq_genotype\tfastq_tier\tfastq_reason\t"
                "fastq_dp\tfastq_quality\tfastq_quality_metric\tfastq_AD\tfastq_allele_fractions\t"
                "reference_genotype\treference_genotype_status\treference_sources\n"
            )
            discrepancies.write(
                "rsid\tchromosome\tposition\tconflict_type\tdante_vcf_genotype\tfastq_genotype\t"
                "reference_genotype\treference_sources\taction\n"
            )
            with output_part.open("wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                    with io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as output:
                        output.write("# Dante WGS reference-complete export; assembly GRCh37\n")
                        if args.calls is not None:
                            output.write(f"# FASTQ tiers included in genotype export: {args.export_min_tier} and higher; Tier C is research-only\n")
                            if args.evidence is not None:
                                output.write(f"# Tiered FASTQ evidence: {args.evidence.name}\n")
                        output.write("# rsid\tchromosome\tposition\tgenotype\n")
                        while target is not None or existing is not None:
                            target_context: ExportRow | None = None
                            if existing is None or (target is not None and target.key < existing.key):
                                assert target is not None
                                target_context = target
                                genotype, selected_source = choose_target(target, args.fastq_source)
                                chosen = ExportRow(target.rsid, target.chromosome, target.position, genotype, target.reason)
                                counters["target_only"] += 1
                                reasons[target.reason] += 1
                                if target.fastq_tier in {"A", "B", "C"}:
                                    counters[f"fastq_tier_{target.fastq_tier}"] += 1
                                if target.genotype != "--":
                                    counters["selected_fastq"] += 1
                                    counters[f"selected_fastq_tier_{target.fastq_tier}"] += 1
                                elif selected_source != "NO_CALL":
                                    counters["selected_reference_source"] += 1
                                if target.reference_status == "CONFLICT":
                                    counters["reference_source_conflict_no_fallback"] += 1
                                if disagree(target.genotype, target.reference_genotype) and target.reference_status == "AVAILABLE":
                                    counters["discrepancies"] += 1
                                    discrepancies.write(
                                        f"{target.rsid}\t{target.chromosome}\t{target.position}\tFASTQ_VS_REFERENCE_SOURCE\t.\t"
                                        f"{target.genotype}\t{target.reference_genotype}\t{target.reference_sources}\t{args.fastq_source}_KEPT\n"
                                    )
                                target = next(targets, None)
                            elif target is None or existing.key < target.key:
                                chosen = existing
                                assert chosen is not None
                                selected_source = "Dante_VCF"
                                counters["base_only"] += 1
                                counters["selected_dante_vcf"] += 1
                                existing = next(base, None)
                            else:
                                assert target is not None and existing is not None
                                target_context = target
                                chosen = existing
                                selected_source = "Dante_VCF"
                                counters["overlap"] += 1
                                counters["selected_dante_vcf"] += 1
                                reasons[target.reason] += 1
                                if target.fastq_tier in {"A", "B", "C"}:
                                    counters[f"fastq_tier_{target.fastq_tier}"] += 1
                                if disagree(existing.genotype, target.genotype):
                                    counters["discrepancies"] += 1
                                    discrepancies.write(
                                        f"{existing.rsid}\t{existing.chromosome}\t{existing.position}\tDANTE_VCF_VS_FASTQ\t"
                                        f"{existing.genotype}\t{target.genotype}\t{target.reference_genotype}\t{target.reference_sources}\tDANTE_VCF_KEPT\n"
                                    )
                                if target.reference_status == "AVAILABLE" and disagree(existing.genotype, target.reference_genotype):
                                    counters["discrepancies"] += 1
                                    discrepancies.write(
                                        f"{existing.rsid}\t{existing.chromosome}\t{existing.position}\tDANTE_VCF_VS_REFERENCE_SOURCE\t"
                                        f"{existing.genotype}\t{target.genotype}\t{target.reference_genotype}\t{target.reference_sources}\tDANTE_VCF_KEPT\n"
                                    )
                                if target.reference_status == "AVAILABLE" and disagree(target.genotype, target.reference_genotype):
                                    counters["discrepancies"] += 1
                                    discrepancies.write(
                                        f"{existing.rsid}\t{existing.chromosome}\t{existing.position}\tFASTQ_VS_REFERENCE_SOURCE\t"
                                        f"{existing.genotype}\t{target.genotype}\t{target.reference_genotype}\t{target.reference_sources}\tDANTE_VCF_KEPT\n"
                                    )
                                target = next(targets, None)
                                existing = next(base, None)

                            output.write(f"{chosen.rsid}\t{chosen.chromosome}\t{chosen.position}\t{chosen.genotype}\n")
                            counters["rows"] += 1
                            if chosen.genotype == "--":
                                counters["no_calls"] += 1
                            else:
                                counters["called"] += 1

                            if target_context is None:
                                provenance.write(
                                    f"{chosen.rsid}\t{chosen.chromosome}\t{chosen.position}\t{chosen.genotype}\t{selected_source}\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n"
                                )
                            else:
                                provenance.write(
                                    f"{chosen.rsid}\t{chosen.chromosome}\t{chosen.position}\t{chosen.genotype}\t{selected_source}\t"
                                    f"{target_context.fastq_candidate_genotype}\t{target_context.fastq_tier}\t{target_context.reason}\t"
                                    f"{target_context.fastq_dp if target_context.fastq_dp is not None else '.'}\t"
                                    f"{target_context.fastq_quality if target_context.fastq_quality is not None else '.'}\t"
                                    f"{target_context.fastq_quality_metric}\t{target_context.fastq_ad}\t{target_context.fastq_allele_fractions}\t"
                                    f"{target_context.reference_genotype}\t{target_context.reference_status}\t{target_context.reference_sources}\n"
                                )
                                if evidence_text is not None and target_context.fastq_tier in {"A", "B", "C"}:
                                    selected = "yes" if selected_source == args.fastq_source and target_context.genotype != "--" else "no"
                                    evidence_text.write(
                                        f"{target_context.rsid}\t{target_context.chromosome}\t{target_context.position}\t"
                                        f"{target_context.fastq_candidate_genotype}\t{target_context.fastq_tier}\t"
                                        f"{target_context.fastq_dp if target_context.fastq_dp is not None else '.'}\t"
                                        f"{target_context.fastq_quality if target_context.fastq_quality is not None else '.'}\t"
                                        f"{target_context.fastq_quality_metric}\t{target_context.fastq_ad}\t"
                                        f"{target_context.fastq_allele_fractions}\t{selected}\n"
                                    )
    finally:
        # Close every layer explicitly. TextIOWrapper closes GzipFile, but GzipFile
        # deliberately does NOT close an externally supplied fileobj. The previous
        # elif-chain therefore left evidence_raw open on Windows and made os.replace()
        # fail with WinError 32 at 97% after a successful merge.
        close_errors = []
        for handle in (evidence_text, evidence_gzip, evidence_raw):
            if handle is None:
                continue
            try:
                if not getattr(handle, "closed", False):
                    handle.close()
            except Exception as exc:
                close_errors.append(exc)
        if close_errors:
            raise close_errors[0]

    progress.emit(97.0, "merge complete | atomically publish export/provenance/conflicts/evidence")
    atomic_replace_with_retry(output_part, args.output)
    atomic_replace_with_retry(provenance_part, args.provenance)
    atomic_replace_with_retry(discrepancy_part, args.discrepancies)
    if evidence_part is not None and args.evidence is not None:
        atomic_replace_with_retry(evidence_part, args.evidence)

    report = [
        "Dante reference-complete export report\n",
        "Assembly: GRCh37\n",
        f"FASTQ evidence source: {args.fastq_source}\n",
        ("Priority: Dante VCF > reference_panel_sources genotype > explicit no-call\n"
         if args.fastq_source == "NONE" else
         "Priority: Dante VCF > FASTQ Tier A/B call > reference_panel_sources genotype > explicit no-call\n"),
        f"Output rows: {counters['rows']}\n",
        f"Called genotypes: {counters['called']}\n",
        f"Explicit no-calls (--): {counters['no_calls']}\n",
        f"Selected from Dante variant-only VCF export: {counters['selected_dante_vcf']}\n",
        f"Selected from tiered {args.fastq_source} calls (Tier {args.export_min_tier} or higher): {counters['selected_fastq']}\n",
        f"FASTQ Tier A — HIGH rows: {counters['fastq_tier_A']}\n",
        f"FASTQ Tier B — SUPPORTIVE rows: {counters['fastq_tier_B']}\n",
        f"FASTQ Tier C — HYPOTHESIS rows (research-only): {counters['fastq_tier_C']}\n",
        f"Selected FASTQ Tier A rows: {counters['selected_fastq_tier_A']}\n",
        f"Selected FASTQ Tier B rows: {counters['selected_fastq_tier_B']}\n",
        f"Selected as fallback from reference_panel_sources: {counters['selected_reference_source']}\n",
        f"Reference-source genotype conflicts not used as fallback: {counters['reference_source_conflict_no_fallback']}\n",
        f"Rows retained from variant-only base only: {counters['base_only']}\n",
        f"Target-only rows added: {counters['target_only']}\n",
        f"Overlapping base/target rows: {counters['overlap']}\n",
        f"Cross-source genotype discrepancies: {counters['discrepancies']}\n",
        f"Provenance table: {args.provenance.name}\n",
        f"Discrepancy table: {args.discrepancies.name}\n",
        f"Main export minimum FASTQ tier: {args.export_min_tier}\n",
        f"Tier A HIGH: DP>={args.tier_a_min_dp}; GQ/REFQ>={args.tier_a_min_quality}; hom AF>={args.tier_a_hom_fraction:.2f}; "
        f"het allele depth>={args.tier_a_min_allele_depth}; het AF={args.tier_a_het_min_fraction:.2f}-{args.tier_a_het_max_fraction:.2f}\n",
        f"Tier B SUPPORTIVE: DP>={args.tier_b_min_dp}; GQ/REFQ>={args.tier_b_min_quality}; hom support>={args.tier_b_hom_min_support}; "
        f"hom AF>={args.tier_b_hom_fraction:.2f}; het allele depth>={args.tier_b_het_min_allele_depth}; "
        f"het AF={args.tier_b_het_min_fraction:.2f}-{args.tier_b_het_max_fraction:.2f}\n",
        f"Tier C HYPOTHESIS: DP>={args.tier_c_min_dp}; GQ/REFQ>={args.tier_c_min_quality}; >=1 observed read for each called allele; research-only\n",
        (f"Tiered FASTQ evidence: {args.evidence.name}\n" if args.evidence is not None else ""),
        "Target call reasons:\n",
    ]
    report.extend(f"  {name}: {count}\n" for name, count in sorted(reasons.items()))
    atomic_text(args.report, "".join(report))
    progress.emit(100.0, f"DONE | {counters['rows']:,} output rows | {counters['called']:,} called | {counters['no_calls']:,} no-calls")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
