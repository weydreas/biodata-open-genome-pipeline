#!/usr/bin/env python3
"""Resumable reference/alternate genotype calling at the deterministic target panel.

Stage 16 checkpoints are independent, validated VCF+TBI chunks over contiguous
segments of the already-sorted target BED. Completed chunks survive crashes,
Docker/WSL transport failures and manual interruption. Final publication occurs
only after all chunks are concatenated, indexed and fully readable.
"""

from __future__ import annotations

import argparse
import atexit
import csv
import gzip
import hashlib
import mmap
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

CHECKPOINT_SCHEMA = "reference-call-checkpoints-v4-known-alleles-ploidy-coordinate-safe"
CHECKPOINT_CALLING_SCHEMA = "reference-calls-v4-known-alleles-ploidy-coordinate-safe"
CALLING_SCHEMA = "reference-calls-v5-known-alleles-ploidy-refq"
DEFAULT_CHUNK_SITES = 250_000
DEFAULT_DISK_RESERVE_GB = 5.0
PIPELINE_RELEASE = "8.4.1"
STAGE_LABEL = "biodata.dante.stage=16"
CANONICAL_CHROMS = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}
DNA_BASES = ("A", "C", "G", "T")


class TransientDockerError(RuntimeError):
    """Docker/WSL transport failure suitable for retrying only the current chunk."""


def mount(source: Path, target: str) -> list[str]:
    return ["--mount", f"type=bind,source={source.resolve()},target={target}"]


def atomic_json(path: Path, payload: dict) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(part, path)


def format_reads(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return str(value)


def format_duration(seconds: float | None) -> str:
    if seconds is None or seconds < 0:
        return "unknown"
    seconds = int(round(seconds))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def format_bytes(value: int) -> str:
    if value >= 1024**3:
        return f"{value / 1024**3:.1f} GiB"
    if value >= 1024**2:
        return f"{value / 1024**2:.1f} MiB"
    if value >= 1024:
        return f"{value / 1024:.1f} KiB"
    return f"{value} B"


def stage16_progress(percent: float, detail: str) -> None:
    bounded = min(100.0, max(0.0, percent))
    print(f"[OVERALL] Stage 16 | {bounded:.2f}% | {detail}", file=sys.stderr, flush=True)


def root_token(root: Path) -> str:
    return hashlib.sha1(str(root.resolve()).encode("utf-8")).hexdigest()[:12]


def stage16_root_label(root: Path) -> str:
    return f"biodata.dante.root={root_token(root)}"


def cleanup_stage16_containers(root: Path) -> int:
    """Remove only orphan containers belonging to this Stage 16 project root."""
    query = subprocess.run(
        [
            "docker", "ps", "-aq",
            "--filter", f"label={STAGE_LABEL}",
            "--filter", f"label={stage16_root_label(root)}",
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    ids = [line.strip() for line in query.stdout.splitlines() if line.strip()]
    if not ids:
        return 0
    subprocess.run(
        ["docker", "rm", "-f", *ids],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return len(ids)


def docker_engine_ready() -> bool:
    probe = subprocess.run(
        ["docker", "info", "--format", "{{.ServerVersion}}"],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return probe.returncode == 0


def wait_for_docker_engine(timeout_seconds: int = 120) -> bool:
    deadline = time.monotonic() + max(1, timeout_seconds)
    while time.monotonic() < deadline:
        if docker_engine_ready():
            return True
        time.sleep(5)
    return docker_engine_ready()


def run(command: list[str], *, transient_125: bool = False) -> None:
    completed = subprocess.run(command, check=False)
    if completed.returncode == 0:
        return
    if transient_125 and completed.returncode == 125:
        raise TransientDockerError(f"Docker command exited 125: {command[0]}")
    raise RuntimeError(f"command exited with {completed.returncode}: {command[0]}")



def capture(command: list[str], *, transient_125: bool = False) -> str:
    completed = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if completed.returncode == 0:
        return completed.stdout
    if transient_125 and completed.returncode == 125:
        raise TransientDockerError(f"Docker command exited 125: {command[0]}")
    tail = (completed.stderr or completed.stdout or "").strip().splitlines()[-5:]
    raise RuntimeError(
        f"command exited with {completed.returncode}: {command[0]}"
        + (f" | {' | '.join(tail)}" if tail else "")
    )


def normalize_chromosome(value: str) -> str | None:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value == "M":
        value = "MT"
    return value if value in CANONICAL_CHROMS else None

def docker_prefix(root: Path, checkpoint_dir: Path, image: str, name: str) -> list[str]:
    return [
        "docker", "run", "--rm",
        "--name", name,
        "--label", STAGE_LABEL,
        "--label", stage16_root_label(root),
        *mount(root, "/work"),
        *mount(checkpoint_dir, "/checkpoint"),
        image,
    ]


def count_target_sites(targets: Path) -> int:
    count = 0
    with targets.open("rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


def skip_target_sites(handle, count: int) -> int:
    skipped = 0
    while skipped < count:
        line = handle.readline()
        if not line:
            break
        if line.strip() and not line.startswith("#"):
            skipped += 1
    return skipped


def write_target_chunk(handle, output: Path, count: int) -> int:
    part = output.with_name(output.name + ".part")
    for path in (part, output):
        try:
            path.unlink()
        except OSError:
            pass
    written = 0
    with part.open("wt", encoding="utf-8", newline="\n") as out:
        while written < count:
            line = handle.readline()
            if not line:
                break
            if not line.strip() or line.startswith("#"):
                continue
            out.write(line.rstrip("\r\n") + "\n")
            written += 1
    if written == 0:
        try:
            part.unlink()
        except OSError:
            pass
        return 0
    os.replace(part, output)
    return written



@dataclass(frozen=True)
class TargetLocus:
    chromosome: str
    position: int
    ref: str
    alts: tuple[str, ...]
    fallback_all_snps: bool
    ref_conflict_rows: int
    callable: bool = True
    uncallable_reason: str = ""


class FastaIndex:
    """Memory-mapped single-base GRCh37 accessor backed by the FASTA .fai."""

    def __init__(self, fasta: Path) -> None:
        self.fasta = fasta
        self.fai = Path(str(fasta) + ".fai")
        if not self.fai.is_file():
            raise RuntimeError(f"FASTA index is missing: {self.fai}")
        self.entries: dict[str, tuple[int, int, int, int]] = {}
        with self.fai.open("rt", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 5:
                    continue
                chrom = normalize_chromosome(fields[0])
                if chrom is None:
                    continue
                self.entries[chrom] = (int(fields[1]), int(fields[2]), int(fields[3]), int(fields[4]))
        self._handle = self.fasta.open("rb")
        self._map = mmap.mmap(self._handle.fileno(), 0, access=mmap.ACCESS_READ)

    def close(self) -> None:
        self._map.close()
        self._handle.close()

    def __enter__(self) -> "FastaIndex":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # type: ignore[no-untyped-def]
        self.close()

    def symbol(self, chromosome: str, position: int) -> str:
        """Return the literal GRCh37 FASTA symbol at a 1-based coordinate.

        Assembly gaps legitimately contain N.  Such loci are not valid
        REF/ALT hypotheses for constrained variant calling, but they are a
        property of the reference assembly rather than a pipeline failure.
        """
        try:
            length, offset, line_bases, line_bytes = self.entries[chromosome]
        except KeyError as exc:
            raise RuntimeError(f"FASTA .fai lacks canonical chromosome {chromosome}") from exc
        if position < 1 or position > length:
            raise RuntimeError(f"Target coordinate outside FASTA contig: {chromosome}:{position}")
        zero = position - 1
        byte_offset = offset + (zero // line_bases) * line_bytes + (zero % line_bases)
        return chr(self._map[byte_offset]).upper()

    def base(self, chromosome: str, position: int) -> str:
        base = self.symbol(chromosome, position)
        if base not in DNA_BASES:
            raise RuntimeError(f"Non-ACGT GRCh37 reference base at {chromosome}:{position}: {base!r}")
        return base


def iter_target_loci(targets_tsv: Path, fasta: FastaIndex) -> Iterator[TargetLocus]:
    """Collapse rsID aliases to one coordinate and provide an explicit allele set.

    Known target alleles are used when they agree with GRCh37. If the panel has
    no trustworthy ALT definition for a SNP coordinate, all three non-reference
    SNP alleles are supplied to bcftools. This gives homozygous-reference calls a
    real alternative hypothesis and therefore a meaningful GQ without pretending
    that an unknown array allele is known.
    """
    with targets_tsv.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chromosome", "position", "ref", "alts"}
        if not required.issubset(set(reader.fieldnames or [])):
            missing = sorted(required - set(reader.fieldnames or []))
            raise RuntimeError(f"Target TSV lacks required columns: {', '.join(missing)}")

        current_key: tuple[int, int] | None = None
        current_chrom = ""
        current_pos = 0
        rows: list[tuple[str, tuple[str, ...]]] = []
        previous_key: tuple[int, int] | None = None

        def build_locus() -> TargetLocus:
            entry = fasta.entries.get(current_chrom)
            if entry is None:
                return TargetLocus(
                    current_chrom, current_pos, ".", tuple(), False, 0, False,
                    "MISSING_FASTA_CONTIG",
                )

            contig_length = int(entry[0])
            if current_pos < 1 or current_pos > contig_length:
                return TargetLocus(
                    current_chrom, current_pos, ".", tuple(), False, 0, False,
                    f"OUTSIDE_FASTA_CONTIG:length={contig_length}",
                )

            reference = fasta.symbol(current_chrom, current_pos)
            if reference not in DNA_BASES:
                # GRCh37 contains legitimate N/gap positions.  A constrained
                # caller cannot formulate a biologically meaningful REF,ALT
                # hypothesis there.  Preserve the panel locus as explicitly
                # uncallable instead of crashing or inventing an allele.
                return TargetLocus(
                    current_chrom,
                    current_pos,
                    reference,
                    tuple(),
                    False,
                    0,
                    False,
                    f"NON_ACGT_REFERENCE:{reference}",
                )

            alts: set[str] = set()
            conflicts = 0
            for row_ref, row_alts in rows:
                if row_ref in DNA_BASES and row_ref != reference:
                    conflicts += 1
                    continue
                for alt in row_alts:
                    if alt in DNA_BASES and alt != reference:
                        alts.add(alt)
            fallback = not alts
            if fallback:
                alts = set(DNA_BASES) - {reference}
            return TargetLocus(
                current_chrom,
                current_pos,
                reference,
                tuple(sorted(alts)),
                fallback,
                conflicts,
            )

        for row in reader:
            chromosome = normalize_chromosome(row.get("chromosome", ""))
            if chromosome is None:
                continue
            try:
                position = int(row.get("position", ""))
            except ValueError as exc:
                raise RuntimeError(f"Invalid target position: {row.get('position')!r}") from exc
            key = (CANONICAL_CHROMS[chromosome], position)
            if previous_key is not None and key < previous_key:
                raise RuntimeError("Target TSV is not in canonical GRCh37 coordinate order")
            previous_key = key
            row_ref = (row.get("ref") or ".").strip().upper()
            row_alts = tuple(
                alt.strip().upper()
                for alt in (row.get("alts") or ".").split(",")
                if alt.strip() and alt.strip() != "."
            )
            if current_key is None:
                current_key = key
                current_chrom = chromosome
                current_pos = position
            elif key != current_key:
                yield build_locus()
                current_key = key
                current_chrom = chromosome
                current_pos = position
                rows = []
            rows.append((row_ref, row_alts))

        if current_key is not None:
            yield build_locus()


def scan_target_reference_contract(
    targets_tsv: Path,
    fasta: FastaIndex,
    uncallable_output: Path,
) -> dict:
    """Audit every unique target coordinate against GRCh37 before calling.

    Non-ACGT reference positions (normally N in assembly gaps) are written to
    a deterministic QC report and excluded from bcftools known-alleles calls.
    They remain present in the master target panel, so Stage 17 will use an
    agreed external-source genotype when available and otherwise emit --.
    """
    part = uncallable_output.with_name(uncallable_output.name + ".part")
    uncallable_output.parent.mkdir(parents=True, exist_ok=True)
    total = callable_count = skipped = 0
    callable_by_chromosome: dict[str, int] = {}
    symbol_counts: dict[str, int] = {}
    reason_counts: dict[str, int] = {}

    with part.open("wt", encoding="utf-8", newline="\n") as out:
        out.write("chromosome\tposition\tgrch37_reference\treason\n")
        for locus in iter_target_loci(targets_tsv, fasta):
            total += 1
            if locus.callable:
                callable_count += 1
                callable_by_chromosome[locus.chromosome] = (
                    callable_by_chromosome.get(locus.chromosome, 0) + 1
                )
                continue
            skipped += 1
            symbol_counts[locus.ref] = symbol_counts.get(locus.ref, 0) + 1
            reason_class = locus.uncallable_reason.split(":", 1)[0]
            reason_counts[reason_class] = reason_counts.get(reason_class, 0) + 1
            out.write(
                f"{locus.chromosome}\t{locus.position}\t{locus.ref}\t"
                f"{locus.uncallable_reason}\n"
            )
    os.replace(part, uncallable_output)
    return {
        "panel_loci": total,
        "callable_loci": callable_count,
        "skipped_uncallable_reference_loci": skipped,
        "skipped_non_acgt_reference_loci": reason_counts.get("NON_ACGT_REFERENCE", 0),
        "skipped_outside_fasta_contig_loci": reason_counts.get("OUTSIDE_FASTA_CONTIG", 0),
        "skipped_missing_fasta_contig_loci": reason_counts.get("MISSING_FASTA_CONTIG", 0),
        "uncallable_reason_counts": reason_counts,
        "uncallable_reference_symbols": symbol_counts,
        "callable_counts_by_chromosome": callable_by_chromosome,
    }


def iter_callable_target_loci(loci: Iterator[TargetLocus]) -> Iterator[TargetLocus]:
    for locus in loci:
        if locus.callable:
            yield locus


def skip_target_loci(loci: Iterator[TargetLocus], count: int) -> int:
    skipped = 0
    while skipped < count:
        try:
            next(loci)
        except StopIteration:
            break
        skipped += 1
    return skipped


def write_call_target_chunk(
    loci: Iterator[TargetLocus],
    bed_output: Path,
    alleles_output: Path,
    count: int,
) -> dict[str, int]:
    bed_part = bed_output.with_name(bed_output.name + ".part")
    alleles_part = alleles_output.with_name(alleles_output.name + ".part")
    for path in (bed_part, alleles_part, bed_output, alleles_output):
        try:
            path.unlink()
        except OSError:
            pass
    written = 0
    fallback = 0
    ref_conflicts = 0
    with bed_part.open("wt", encoding="utf-8", newline="\n") as bed, alleles_part.open(
        "wt", encoding="utf-8", newline="\n"
    ) as alleles:
        while written < count:
            try:
                locus = next(loci)
            except StopIteration:
                break
            if not locus.callable:
                # Main uses iter_callable_target_loci(), but retain this guard so
                # future callers cannot accidentally send N-reference loci to bcftools.
                continue
            bed.write(f"{locus.chromosome}\t{locus.position - 1}\t{locus.position}\n")
            alleles.write(
                f"{locus.chromosome}\t{locus.position}\t{locus.ref},{','.join(locus.alts)}\n"
            )
            written += 1
            fallback += int(locus.fallback_all_snps)
            ref_conflicts += locus.ref_conflict_rows
    if written <= 0:
        for path in (bed_part, alleles_part):
            try:
                path.unlink()
            except OSError:
                pass
        return {"count": 0, "fallback": 0, "ref_conflicts": 0}
    os.replace(bed_part, bed_output)
    os.replace(alleles_part, alleles_output)
    return {"count": written, "fallback": fallback, "ref_conflicts": ref_conflicts}


def samtools_capture(root: Path, image: str, bam: Path, args: list[str], name_suffix: str) -> str:
    name = f"dante-stage16-{root_token(root)}-{os.getpid()}-{name_suffix}"
    return capture(
        [
            "docker", "run", "--rm",
            "--name", name,
            "--label", STAGE_LABEL,
            "--label", stage16_root_label(root),
            *mount(root, "/work"),
            "--entrypoint", "samtools",
            image,
            *args,
            "/work/" + bam.relative_to(root).as_posix(),
        ],
        transient_125=True,
    )


def bam_sample_name(root: Path, image: str, bam: Path) -> str:
    header = samtools_capture(root, image, bam, ["view", "-H"], "header")
    samples: set[str] = set()
    for line in header.splitlines():
        if not line.startswith("@RG\t"):
            continue
        for field in line.split("\t"):
            if field.startswith("SM:") and len(field) > 3:
                samples.add(field[3:])
    if len(samples) != 1:
        raise RuntimeError(
            f"Expected exactly one BAM sample (SM) in @RG headers, found {len(samples)}: {sorted(samples)}"
        )
    return next(iter(samples))


def target_counts_by_chromosome(targets_bed: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with targets_bed.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom = normalize_chromosome(line.split("\t", 1)[0])
            if chrom is not None:
                result[chrom] = result.get(chrom, 0) + 1
    return result


def mapped_reads_by_chromosome(root: Path, image: str, bam: Path) -> dict[str, int]:
    text = samtools_capture(root, image, bam, ["idxstats"], "idxstats")
    result: dict[str, int] = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        chrom = normalize_chromosome(fields[0])
        if chrom is None:
            continue
        try:
            result[chrom] = int(fields[2])
        except ValueError:
            continue
    return result


def infer_sex_from_normalized_depth(
    mapped: dict[str, int], targets: dict[str, int]
) -> dict[str, float | str]:
    auto_reads = sum(mapped.get(str(i), 0) for i in range(1, 23))
    auto_targets = sum(targets.get(str(i), 0) for i in range(1, 23))
    if auto_reads <= 0 or auto_targets <= 0:
        raise RuntimeError("Cannot infer sample sex: autosomal target/read denominator is empty")
    auto_rate = auto_reads / auto_targets
    x_targets = targets.get("X", 0)
    y_targets = targets.get("Y", 0)
    x_rate = mapped.get("X", 0) / max(1, x_targets)
    y_rate = mapped.get("Y", 0) / max(1, y_targets)
    x_ratio = x_rate / auto_rate
    y_ratio = y_rate / auto_rate
    if y_targets < 10:
        inferred = "AMBIGUOUS"
    else:
        inferred = "M" if y_ratio >= 0.15 else ("F" if y_ratio <= 0.03 else "AMBIGUOUS")
    return {
        "inferred": inferred,
        "autosome_reads_per_target": auto_rate,
        "x_to_autosome_ratio": x_ratio,
        "y_to_autosome_ratio": y_ratio,
    }


def normalize_sample_sex(value: str) -> str:
    normalized = value.strip().upper()
    aliases = {
        "": "AUTO", "AUTO": "AUTO",
        "M": "M", "MALE": "M", "XY": "M",
        "F": "F", "FEMALE": "F", "XX": "F",
    }
    if normalized not in aliases:
        raise RuntimeError("DANTE_SAMPLE_SEX must be AUTO, M/MALE/XY, or F/FEMALE/XX")
    return aliases[normalized]


def resolve_sample_sex(
    requested: str,
    mapped: dict[str, int],
    targets: dict[str, int],
) -> dict:
    requested = normalize_sample_sex(requested)
    inference = infer_sex_from_normalized_depth(mapped, targets)
    inferred = str(inference["inferred"])
    if requested == "AUTO":
        if inferred == "AMBIGUOUS":
            raise RuntimeError(
                "Automatic sex/ploidy inference is ambiguous. Set DANTE_SAMPLE_SEX=M or F in .env and rerun item 16."
            )
        resolved = inferred
        source = "auto_target_bam_depth"
    else:
        resolved = requested
        source = "explicit_env"
        if inferred in {"M", "F"} and inferred != resolved:
            raise RuntimeError(
                f"Explicit sample sex {resolved} conflicts with strong BAM inference {inferred}. "
                "Check sample identity or DANTE_SAMPLE_SEX before calling sex chromosomes."
            )
    return {
        "sex_requested": requested,
        "sex_resolved": resolved,
        "sex_source": source,
        **inference,
    }


def write_sample_sex_file(path: Path, sample_name: str, sex: str) -> None:
    part = path.with_name(path.name + ".part")
    part.write_text(f"{sample_name}\t{sex}\n", encoding="utf-8")
    os.replace(part, path)


def expected_ploidy(chromosome: str, position: int, sex: str) -> int:
    if chromosome == "MT":
        return 1
    if sex == "F":
        if chromosome == "Y":
            return 0
        return 2
    if chromosome == "Y":
        return 1
    if chromosome == "X" and (position <= 60000 or 2699521 <= position <= 154931043):
        return 1
    return 2


def pl_reference_quality(pl_value: str | None, allele_count: int, ploidy: int) -> int | None:
    """Return a PL-derived quality for a homozygous-reference genotype.

    bcftools multiallelic calling intentionally does not emit FORMAT/GQ when the
    selected site is non-variant. With -A, however, it preserves the constrained
    PL likelihoods. For a reference call we therefore compare the reference
    genotype likelihood with the best non-reference hypothesis.

    Diploid: compare PL(0/0) with every other genotype.
    Haploid: accept either a true haploid PL vector (Number=R) or the diploid
    likelihood layout used internally by mpileup; in the latter case compare
    0/0 only with alternate homozygotes and ignore heterozygous states.
    """
    if not pl_value or pl_value in {".", ""} or allele_count < 2 or ploidy not in {1, 2}:
        return None
    try:
        values = [int(x) for x in pl_value.split(",")]
    except ValueError:
        return None
    if not values or values[0] < 0:
        return None

    competitors: list[int] = []
    if ploidy == 2:
        expected = allele_count * (allele_count + 1) // 2
        if len(values) < expected:
            return None
        competitors = values[1:expected]
    else:
        # Native haploid layout: one PL per allele.
        if len(values) == allele_count:
            competitors = values[1:]
        else:
            # mpileup/call commonly retains diploid Number=G PLs even where the
            # final genotype is haploid. Use only homozygous allele states.
            expected = allele_count * (allele_count + 1) // 2
            if len(values) < expected:
                return None
            for allele in range(1, allele_count):
                idx = allele * (allele + 1) // 2 + allele
                if idx >= len(values):
                    return None
                competitors.append(values[idx])

    if not competitors:
        return None
    return max(0, min(competitors) - values[0])


def validate_call_contract(path: Path, sex: str) -> dict[str, int | float]:
    total = hom_ref = hom_ref_gq = hom_ref_pl_refq = 0
    alt_carried = alt_with_gq = missing = ploidy_mismatch = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = normalize_chromosome(fields[0])
            if chrom is None:
                continue
            position = int(fields[1])
            allele_count = 1 + (0 if fields[4] in {"", "."} else len(fields[4].split(",")))
            fmt = fields[8].split(":")
            values = fields[9].split(":")
            sample = dict(zip(fmt, values))
            gt = sample.get("GT", ".").replace("|", "/")
            total += 1
            if gt in {"", ".", "./."} or "." in gt:
                missing += 1
                continue
            parts = gt.split("/")
            try:
                indexes = [int(item) for item in parts]
            except ValueError:
                missing += 1
                continue
            expected = expected_ploidy(chrom, position, sex)
            if expected == 0:
                ploidy_mismatch += 1
            elif len(indexes) != expected:
                ploidy_mismatch += 1

            gq_raw = sample.get("GQ", ".")
            has_gq = gq_raw not in {"", "."}
            if all(index == 0 for index in indexes):
                hom_ref += 1
                if has_gq:
                    hom_ref_gq += 1
                else:
                    refq = pl_reference_quality(sample.get("PL"), allele_count, len(indexes))
                    if refq is not None:
                        hom_ref_pl_refq += 1
            elif any(index > 0 for index in indexes):
                alt_carried += 1
                if has_gq:
                    alt_with_gq += 1

    hom_ref_confident = hom_ref_gq + hom_ref_pl_refq
    hom_ref_confidence_fraction = hom_ref_confident / hom_ref if hom_ref else 1.0
    alt_gq_fraction = alt_with_gq / alt_carried if alt_carried else 1.0

    if hom_ref >= 100 and hom_ref_confidence_fraction < 0.95:
        raise RuntimeError(
            "Known-alleles call contract failed: only "
            f"{hom_ref_confidence_fraction:.2%} of homozygous-reference calls have "
            "either native GQ or a PL-derived reference quality."
        )
    if alt_carried >= 100 and alt_gq_fraction < 0.95:
        raise RuntimeError(
            f"Variant-call contract failed: only {alt_gq_fraction:.2%} of ALT-carried calls have native GQ."
        )
    if ploidy_mismatch:
        raise RuntimeError(
            f"Ploidy contract failed: {ploidy_mismatch:,} called records disagree with resolved sex/ploidy rules."
        )
    return {
        "total_records": total,
        "hom_ref_records": hom_ref,
        "hom_ref_with_native_gq": hom_ref_gq,
        "hom_ref_with_pl_derived_refq": hom_ref_pl_refq,
        "hom_ref_with_reference_quality": hom_ref_confident,
        "hom_ref_reference_quality_fraction": hom_ref_confidence_fraction,
        "alt_carried_records": alt_carried,
        "alt_carried_with_native_gq": alt_with_gq,
        "alt_carried_gq_fraction": alt_gq_fraction,
        "missing_genotype_records": missing,
        "ploidy_mismatch_records": ploidy_mismatch,
    }

def expected_identity(
    args,
    bam: Path,
    targets_bed: Path,
    targets_tsv: Path,
    fasta: Path,
    ploidy_file: Path,
    sample_name: str,
    resolved_sex: str,
) -> dict:
    return {
        "schema": CHECKPOINT_SCHEMA,
        "calling_schema": CHECKPOINT_CALLING_SCHEMA,
        "calling_model": "bcftools-m-known-alleles-gq-ploidy-coordinate-safe-v4",
        "reference_coordinate_policy": "skip-non-acgt-out-of-range-missing-contig-and-report-v2",
        "assembly": "GRCh37",
        "bam_name": bam.name,
        "bam_length": bam.stat().st_size,
        "bam_sha256": args.bam_sha256,
        "targets_bed_name": targets_bed.name,
        "targets_bed_length": targets_bed.stat().st_size,
        "targets_bed_sha256": args.targets_sha256,
        "targets_tsv_name": targets_tsv.name,
        "targets_tsv_length": targets_tsv.stat().st_size,
        "targets_tsv_sha256": args.targets_tsv_sha256,
        "fasta_name": fasta.name,
        "fasta_length": fasta.stat().st_size,
        "fasta_sha256": args.fasta_sha256,
        "ploidy_name": ploidy_file.name,
        "ploidy_sha256": args.ploidy_sha256,
        "sample_name": sample_name,
        "sample_sex": resolved_sex,
        "bcftools_contract": args.bcftools_contract,
        "samtools_image": args.samtools_image,
        "samtools_contract": args.samtools_contract,
        "min_mapq": max(0, args.min_mapq),
        "min_baseq": max(0, args.min_baseq),
        "checkpoint_format": "bgzip-vcf-tbi-known-alleles-coordinate-safe-v4",
    }

def manifest_matches(manifest: dict, identity: dict) -> bool:
    return manifest.get("schema") == CHECKPOINT_SCHEMA and manifest.get("identity") == identity


def valid_completed_chunks(manifest: dict, checkpoint_dir: Path) -> list[dict]:
    valid: list[dict] = []
    expected_index = 0
    expected_start = 0
    for row in manifest.get("completed_chunks", []):
        try:
            index = int(row["index"])
            start = int(row["target_start"])
            count = int(row["target_count"])
            vcf = checkpoint_dir / str(row["vcf_name"])
            tbi = Path(str(vcf) + ".tbi")
            if index != expected_index or start != expected_start or count <= 0:
                break
            if not vcf.is_file() or not tbi.is_file():
                break
            if vcf.stat().st_size != int(row["vcf_length"]):
                break
            if tbi.stat().st_size != int(row["tbi_length"]):
                break
            if vcf.stat().st_size <= 0 or tbi.stat().st_size <= 0:
                break
        except (KeyError, TypeError, ValueError, OSError):
            break
        valid.append(row)
        expected_index += 1
        expected_start += count
    return valid


def clean_transients(checkpoint_dir: Path) -> None:
    for pattern in ("*.part", "*.part.vcf.gz", "*.part.vcf.gz.tbi", "current_targets_*.bed", "current_alleles_*.tsv", "current_alleles_*.tsv.gz", "current_alleles_*.tsv.gz.tbi", "sample_sex.tsv"):
        for path in checkpoint_dir.glob(pattern):
            if path.is_file():
                try:
                    path.unlink()
                except OSError:
                    pass


def checkpoint_bytes(chunks: Iterable[dict], checkpoint_dir: Path) -> int:
    total = 0
    for row in chunks:
        vcf = checkpoint_dir / str(row["vcf_name"])
        tbi = Path(str(vcf) + ".tbi")
        for path in (vcf, tbi):
            try:
                total += path.stat().st_size
            except OSError:
                pass
    return total


def disk_guard(
    checkpoint_dir: Path,
    chunks: list[dict],
    completed_sites: int,
    total_sites: int,
    reserve_gb: float,
) -> None:
    usage = shutil.disk_usage(checkpoint_dir)
    used = checkpoint_bytes(chunks, checkpoint_dir)
    reserve = int(max(0.0, reserve_gb) * 1024**3)
    projected_remaining = 0
    if completed_sites > 0 and used > 0 and total_sites > completed_sites:
        projected_remaining = int((used / completed_sites) * (total_sites - completed_sites))
    stage16_progress(
        5.0 + 85.0 * (completed_sites / max(1, total_sites)),
        f"DISK {format_bytes(usage.free)} free | checkpoints {format_bytes(used)} | "
        f"reserve {format_bytes(reserve)} | projected remaining checkpoints {format_bytes(projected_remaining)}",
    )
    if usage.free <= reserve:
        raise RuntimeError(
            "SAFE STOP: free disk is already at/below the Stage 16 reserve. "
            "Completed call checkpoints were preserved; free/move storage and rerun item 16."
        )
    if projected_remaining and usage.free - reserve < projected_remaining:
        raise RuntimeError(
            "SAFE STOP: projected Stage 16 checkpoint growth exceeds free disk after reserve. "
            "Completed call checkpoints were preserved; free/move storage and rerun item 16."
        )


def call_chunk(
    *,
    root: Path,
    data_dir: Path,
    checkpoint_dir: Path,
    fasta_name: str,
    bam: Path,
    bed: Path,
    alleles_tsv: Path,
    ploidy_file: Path,
    sample_file: Path,
    output_vcf: Path,
    image: str,
    threads: int,
    min_mapq: int,
    min_baseq: int,
    chunk_index: int,
    attempt: int,
) -> None:
    part = output_vcf.with_name(output_vcf.name.replace(".vcf.gz", ".part.vcf.gz"))
    part_tbi = Path(str(part) + ".tbi")
    final_tbi = Path(str(output_vcf) + ".tbi")
    alleles_gz = Path(str(alleles_tsv) + ".gz")
    alleles_tbi = Path(str(alleles_gz) + ".tbi")
    for path in (part, part_tbi, output_vcf, final_tbi, alleles_gz, alleles_tbi):
        try:
            path.unlink()
        except OSError:
            pass

    # bcftools call -C alleles requires a bgzip-compressed and tabix-indexed
    # CHROM,POS,REF,ALT-list target file. Build it from the exact current chunk.
    prep_name = f"dante-stage16-{root_token(root)}-{os.getpid()}-als-{chunk_index:04d}-{attempt}"
    prep_pipeline = (
        "set -euo pipefail; "
        f"bgzip -c '/checkpoint/{alleles_tsv.name}' > '/checkpoint/{alleles_gz.name}'; "
        f"tabix -f -s1 -b2 -e2 '/checkpoint/{alleles_gz.name}'"
    )
    run(
        [
            "docker", "run", "--rm",
            "--name", prep_name,
            "--label", STAGE_LABEL,
            "--label", stage16_root_label(root),
            *mount(checkpoint_dir, "/checkpoint"),
            "--entrypoint", "bash",
            image, "-lc", prep_pipeline,
        ],
        transient_125=True,
    )
    if not alleles_gz.is_file() or not alleles_tbi.is_file():
        raise RuntimeError(f"Known-alleles target index was not created for chunk {chunk_index + 1}")

    bam_container = "/work/" + bam.relative_to(root).as_posix()
    bed_container = "/checkpoint/" + bed.name
    alleles_container = "/checkpoint/" + alleles_gz.name
    sample_container = "/checkpoint/" + sample_file.name
    ploidy_container = "/work/" + ploidy_file.relative_to(root).as_posix()
    part_container = "/checkpoint/" + part.name
    name = f"dante-stage16-{root_token(root)}-{os.getpid()}-{chunk_index:04d}-{attempt}"
    pipeline = (
        "set -o pipefail; "
        f"bcftools mpileup -Ou --threads {threads} -f '/data/fasta/{fasta_name}' "
        f"-R '{bed_container}' -q {min_mapq} -Q {min_baseq} "
        f"-a FORMAT/AD,FORMAT/DP '{bam_container}' | "
        f"bcftools call -m -C alleles -T '{alleles_container}' -A -f GQ "
        f"--ploidy-file '{ploidy_container}' -S '{sample_container}' "
        f"--threads {threads} -Oz -o '{part_container}'"
    )
    run(
        [
            "docker", "run", "--rm",
            "--name", name,
            "--label", STAGE_LABEL,
            "--label", stage16_root_label(root),
            *mount(root, "/work"),
            *mount(data_dir, "/data"),
            *mount(checkpoint_dir, "/checkpoint"),
            "--entrypoint", "bash",
            image, "-lc", pipeline,
        ],
        transient_125=True,
    )
    if not part.is_file() or part.stat().st_size <= 0:
        raise RuntimeError(f"bcftools did not create a non-empty checkpoint VCF for chunk {chunk_index + 1}")

    index_name = f"dante-stage16-{root_token(root)}-{os.getpid()}-idx-{chunk_index:04d}-{attempt}"
    run(
        [
            "docker", "run", "--rm",
            "--name", index_name,
            "--label", STAGE_LABEL,
            "--label", stage16_root_label(root),
            *mount(checkpoint_dir, "/checkpoint"),
            "--entrypoint", "bcftools",
            image, "index", "--tbi", "--force", part_container,
        ],
        transient_125=True,
    )
    if not part_tbi.is_file() or part_tbi.stat().st_size <= 0:
        raise RuntimeError(f"bcftools did not create a checkpoint tabix index for chunk {chunk_index + 1}")

    validate_name = f"dante-stage16-{root_token(root)}-{os.getpid()}-val-{chunk_index:04d}-{attempt}"
    run(
        [
            "docker", "run", "--rm",
            "--name", validate_name,
            "--label", STAGE_LABEL,
            "--label", stage16_root_label(root),
            *mount(checkpoint_dir, "/checkpoint"),
            "--entrypoint", "bcftools",
            image, "view", "--output-type", "v", "--output", "/dev/null", part_container,
        ],
        transient_125=True,
    )

    os.replace(part, output_vcf)
    os.replace(part_tbi, final_tbi)

def finalize_calls(
    *,
    root: Path,
    checkpoint_dir: Path,
    output: Path,
    image: str,
    chunks: list[dict],
) -> tuple[Path, Path]:
    part = output.with_name(output.name.replace(".vcf.gz", ".part.vcf.gz"))
    part_tbi = Path(str(part) + ".tbi")
    final_tbi = Path(str(output) + ".tbi")
    for path in (part, part_tbi):
        try:
            path.unlink()
        except OSError:
            pass

    chunk_paths = [f"/checkpoint/{row['vcf_name']}" for row in chunks]
    part_container = "/work/" + part.relative_to(root).as_posix()
    if len(chunks) == 1:
        shutil.copy2(checkpoint_dir / str(chunks[0]["vcf_name"]), part)
    else:
        concat_name = f"dante-stage16-{root_token(root)}-{os.getpid()}-concat"
        run(
            [
                "docker", "run", "--rm",
                "--name", concat_name,
                "--label", STAGE_LABEL,
                "--label", stage16_root_label(root),
                *mount(root, "/work"),
                *mount(checkpoint_dir, "/checkpoint"),
                "--entrypoint", "bcftools",
                image, "concat", "-Oz", "-o", part_container, *chunk_paths,
            ],
            transient_125=True,
        )
    if not part.is_file() or part.stat().st_size <= 0:
        raise RuntimeError("bcftools concat did not create a non-empty final partial VCF")

    index_name = f"dante-stage16-{root_token(root)}-{os.getpid()}-finalidx"
    run(
        [
            "docker", "run", "--rm",
            "--name", index_name,
            "--label", STAGE_LABEL,
            "--label", stage16_root_label(root),
            *mount(root, "/work"),
            "--entrypoint", "bcftools",
            image, "index", "--tbi", "--force", part_container,
        ],
        transient_125=True,
    )

    validate_name = f"dante-stage16-{root_token(root)}-{os.getpid()}-finalval"
    run(
        [
            "docker", "run", "--rm",
            "--name", validate_name,
            "--label", STAGE_LABEL,
            "--label", stage16_root_label(root),
            *mount(root, "/work"),
            "--entrypoint", "bcftools",
            image, "view", "--output-type", "v", "--output", "/dev/null", part_container,
        ],
        transient_125=True,
    )
    if not part_tbi.is_file() or part_tbi.stat().st_size <= 0:
        raise RuntimeError("bcftools did not create the final tabix index")

    return part, part_tbi


def retry_transient(root: Path, label: str, fn, max_attempts: int = 3) -> None:
    for attempt in range(1, max_attempts + 1):
        try:
            fn(attempt)
            return
        except TransientDockerError as exc:
            if attempt >= max_attempts:
                raise RuntimeError(
                    f"Docker transport failed repeatedly during {label}; completed Stage 16 checkpoints are preserved."
                ) from exc
            print(
                f"[RECOVERY] Stage 16 | {label} | retry {attempt + 1}/{max_attempts}",
                file=sys.stderr,
                flush=True,
            )
            cleanup_stage16_containers(root)
            if not wait_for_docker_engine(120):
                raise RuntimeError(
                    "Docker engine did not recover within 120 seconds; completed Stage 16 checkpoints are preserved."
                ) from exc


def main() -> int:
    if len(sys.argv) == 2 and sys.argv[1] == "--print-contract":
        print(f"{PIPELINE_RELEASE}|{CALLING_SCHEMA}|{CHECKPOINT_SCHEMA}")
        return 0

    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--fasta-name", required=True)
    parser.add_argument("--fasta-sha256", required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--bam-sha256", required=True)
    parser.add_argument("--targets-bed", type=Path, required=True)
    parser.add_argument("--targets-sha256", required=True)
    parser.add_argument("--targets-tsv", type=Path, required=True)
    parser.add_argument("--targets-tsv-sha256", required=True)
    parser.add_argument("--ploidy-file", type=Path, required=True)
    parser.add_argument("--ploidy-sha256", required=True)
    parser.add_argument("--sample-sex", default="AUTO")
    parser.add_argument("--samtools-image", required=True)
    parser.add_argument("--samtools-contract", required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    parser.add_argument("--uncallable-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    parser.add_argument("--bcftools-image", required=True)
    parser.add_argument("--bcftools-contract", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--chunk-sites", type=int, default=DEFAULT_CHUNK_SITES)
    parser.add_argument("--disk-reserve-gb", type=float, default=DEFAULT_DISK_RESERVE_GB)
    args = parser.parse_args()

    root = args.root.resolve()
    data_dir = args.data_dir.resolve()
    bam = args.bam.resolve()
    targets_bed = args.targets_bed.resolve()
    targets_tsv = args.targets_tsv.resolve()
    ploidy_file = args.ploidy_file.resolve()
    metadata_output = args.metadata_output.resolve()
    uncallable_output = args.uncallable_output.resolve()
    output = args.output.resolve()
    checkpoint_dir = args.checkpoint_dir.resolve()
    fasta = data_dir / "fasta" / args.fasta_name
    for path in (
        bam, Path(str(bam) + ".bai"), targets_bed, targets_tsv,
        fasta, Path(str(fasta) + ".fai"), ploidy_file,
    ):
        if not path.is_file():
            raise SystemExit(f"Required file is missing: {path}")
    for path in (bam, targets_bed, targets_tsv, ploidy_file, output.parent, metadata_output.parent, uncallable_output.parent):
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise SystemExit(f"Path must be inside pipeline root: {path}") from exc

    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata_output.parent.mkdir(parents=True, exist_ok=True)
    uncallable_output.parent.mkdir(parents=True, exist_ok=True)
    threads = max(1, args.threads)
    chunk_sites = max(1, args.chunk_sites)
    min_mapq = max(0, args.min_mapq)
    min_baseq = max(0, args.min_baseq)

    cleanup_stage16_containers(root)
    atexit.register(cleanup_stage16_containers, root)
    clean_transients(checkpoint_dir)

    legacy_part = output.with_name(output.name.replace(".vcf.gz", ".part.vcf.gz"))
    legacy_part_tbi = Path(str(legacy_part) + ".tbi")
    for path in (legacy_part, legacy_part_tbi):
        if path.exists():
            path.unlink()
            print(f"[MIGRATION] Removed legacy non-resumable Stage 16 partial: {path.name}", file=sys.stderr, flush=True)

    total_panel_sites = count_target_sites(targets_bed)
    if total_panel_sites <= 0:
        raise SystemExit("Target BED contains no sites")

    # Audit the reference contract before any bcftools work.  GRCh37 N/gap
    # coordinates are valid panel metadata but cannot be used as REF alleles.
    with FastaIndex(fasta) as fasta_audit:
        reference_contract = scan_target_reference_contract(
            targets_tsv, fasta_audit, uncallable_output
        )
    if int(reference_contract["panel_loci"]) != total_panel_sites:
        raise RuntimeError(
            "Target BED/TSV unique-locus count mismatch during GRCh37 reference audit: "
            f"BED={total_panel_sites:,}, TSV={int(reference_contract['panel_loci']):,}"
        )
    total_sites = int(reference_contract["callable_loci"])
    skipped_uncallable = int(reference_contract["skipped_uncallable_reference_loci"])
    skipped_non_acgt = int(reference_contract["skipped_non_acgt_reference_loci"])
    skipped_outside = int(reference_contract["skipped_outside_fasta_contig_loci"])
    skipped_missing_contig = int(reference_contract["skipped_missing_fasta_contig_loci"])
    if total_sites <= 0:
        raise SystemExit("No callable A/C/G/T GRCh37 target loci remain after reference audit")
    print(
        f"[REFERENCE] panel={format_reads(total_panel_sites)} loci | "
        f"callable A/C/G/T={format_reads(total_sites)} | "
        f"skipped uncallable={skipped_uncallable:,} "
        f"(non-ACGT={skipped_non_acgt:,}, outside-contig={skipped_outside:,}, "
        f"missing-contig={skipped_missing_contig:,}) | report={uncallable_output.name}",
        file=sys.stderr,
        flush=True,
    )

    sample_name = bam_sample_name(root, args.samtools_image, bam)
    # Sex inference is based on the target-filtered BAM, which was produced
    # from the full target BED.  Keep the denominator aligned to that same BED
    # rather than mixing callable-only loci with full-BAM read counts.
    target_chrom_counts = target_counts_by_chromosome(targets_bed)
    mapped_chrom_counts = mapped_reads_by_chromosome(root, args.samtools_image, bam)
    sex_info = resolve_sample_sex(args.sample_sex, mapped_chrom_counts, target_chrom_counts)
    resolved_sex = str(sex_info["sex_resolved"])
    sample_file = checkpoint_dir / "sample_sex.tsv"
    write_sample_sex_file(sample_file, sample_name, resolved_sex)
    print(
        f"[PLOIDY] sample={sample_name} | sex={resolved_sex} ({sex_info['sex_source']}) | "
        f"X/autosome={float(sex_info['x_to_autosome_ratio']):.3f} | "
        f"Y/autosome={float(sex_info['y_to_autosome_ratio']):.3f}",
        file=sys.stderr,
        flush=True,
    )

    identity = expected_identity(
        args, bam, targets_bed, targets_tsv, fasta, ploidy_file, sample_name, resolved_sex
    )
    manifest_path = checkpoint_dir / "checkpoint_manifest.json"
    manifest: dict
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            manifest = {}
    else:
        manifest = {}

    if not manifest_matches(manifest, identity):
        for pattern in (
            "call_chunk_*.vcf.gz", "call_chunk_*.vcf.gz.tbi", "current_targets_*.bed",
            "current_alleles_*.tsv", "current_alleles_*.tsv.gz", "current_alleles_*.tsv.gz.tbi",
        ):
            for path in checkpoint_dir.glob(pattern):
                if path.is_file():
                    try:
                        path.unlink()
                    except OSError:
                        pass
        manifest = {
            "schema": CHECKPOINT_SCHEMA,
            "identity": identity,
            "status": "in_progress",
            "total_target_sites": total_sites,
            "total_panel_loci": total_panel_sites,
            "skipped_uncallable_reference_loci": skipped_uncallable,
            "uncallable_reason_counts": reference_contract["uncallable_reason_counts"],
            "completed_chunks": [],
        }
        atomic_json(manifest_path, manifest)

    chunks = valid_completed_chunks(manifest, checkpoint_dir)
    if len(chunks) != len(manifest.get("completed_chunks", [])):
        manifest["completed_chunks"] = chunks
        manifest["status"] = "in_progress"
        atomic_json(manifest_path, manifest)

    completed_sites = sum(int(row["target_count"]) for row in chunks)
    if completed_sites > total_sites:
        raise RuntimeError("Stage 16 checkpoint manifest exceeds current target BED size")

    if completed_sites:
        stage16_progress(
            5.0 + 85.0 * (completed_sites / total_sites),
            f"RESUME | {len(chunks)} completed known-alleles checkpoints | "
            f"{format_reads(completed_sites)} / {format_reads(total_sites)} target loci preserved",
        )
    else:
        stage16_progress(
            5.0,
            f"bcftools known-alleles/ploidy caller ready | {format_reads(total_sites)} callable loci | "
            f"{skipped_uncallable:,} uncallable GRCh37 loci skipped/reported | resumable chunks",
        )

    run_started = time.monotonic()
    initial_completed_sites = completed_sites
    total_fallback = sum(int(row.get("fallback_all_snp_loci", 0)) for row in chunks)
    total_ref_conflicts = sum(int(row.get("target_ref_conflict_rows", 0)) for row in chunks)

    with FastaIndex(fasta) as fasta_index:
        loci = iter_callable_target_loci(iter_target_loci(targets_tsv, fasta_index))
        skipped = skip_target_loci(loci, completed_sites)
        if skipped != completed_sites:
            raise RuntimeError("Target TSV ended while scanning to the Stage 16 checkpoint boundary")

        chunk_index = len(chunks)
        while completed_sites < total_sites:
            wanted = min(chunk_sites, total_sites - completed_sites)
            current_bed = checkpoint_dir / f"current_targets_{chunk_index:05d}.bed"
            current_alleles = checkpoint_dir / f"current_alleles_{chunk_index:05d}.tsv"
            chunk_meta = write_call_target_chunk(loci, current_bed, current_alleles, wanted)
            current_count = int(chunk_meta["count"])
            if current_count <= 0:
                raise RuntimeError("Target TSV ended before the expected number of Stage 16 loci")
            total_fallback += int(chunk_meta["fallback"])
            total_ref_conflicts += int(chunk_meta["ref_conflicts"])
            final_vcf = checkpoint_dir / f"call_chunk_{chunk_index:05d}.vcf.gz"

            disk_guard(checkpoint_dir, chunks, completed_sites, total_sites, args.disk_reserve_gb)
            stage16_progress(
                5.0 + 85.0 * (completed_sites / total_sites),
                f"CHUNK {chunk_index + 1} | call {format_reads(current_count)} loci against explicit REF/ALT hypotheses | checkpoint after validation",
            )

            def attempt_chunk(attempt: int) -> None:
                call_chunk(
                    root=root,
                    data_dir=data_dir,
                    checkpoint_dir=checkpoint_dir,
                    fasta_name=args.fasta_name,
                    bam=bam,
                    bed=current_bed,
                    alleles_tsv=current_alleles,
                    ploidy_file=ploidy_file,
                    sample_file=sample_file,
                    output_vcf=final_vcf,
                    image=args.bcftools_image,
                    threads=threads,
                    min_mapq=min_mapq,
                    min_baseq=min_baseq,
                    chunk_index=chunk_index,
                    attempt=attempt,
                )

            retry_transient(root, f"Stage 16 known-alleles chunk {chunk_index + 1}", attempt_chunk)

            tbi = Path(str(final_vcf) + ".tbi")
            row = {
                "index": chunk_index,
                "target_start": completed_sites,
                "target_count": current_count,
                "vcf_name": final_vcf.name,
                "vcf_length": final_vcf.stat().st_size,
                "tbi_length": tbi.stat().st_size,
                "fallback_all_snp_loci": int(chunk_meta["fallback"]),
                "target_ref_conflict_rows": int(chunk_meta["ref_conflicts"]),
            }
            chunks.append(row)
            completed_sites += current_count
            manifest["completed_chunks"] = chunks
            manifest["status"] = "in_progress"
            manifest["total_target_sites"] = total_sites
            manifest["current_chunk_sites"] = chunk_sites
            atomic_json(manifest_path, manifest)
            for path in (
                current_bed,
                current_alleles,
                Path(str(current_alleles) + ".gz"),
                Path(str(current_alleles) + ".gz.tbi"),
            ):
                try:
                    path.unlink()
                except OSError:
                    pass

            elapsed = max(0.001, time.monotonic() - run_started)
            done_this_run = max(0, completed_sites - initial_completed_sites)
            speed = done_this_run / elapsed if elapsed > 0 else 0.0
            remaining = max(0, total_sites - completed_sites)
            eta = remaining / speed if speed > 0 else None
            stage16_progress(
                5.0 + 85.0 * (completed_sites / total_sites),
                f"CHECKPOINT SAVED | chunk {chunk_index + 1} | "
                f"{format_reads(completed_sites)} / {format_reads(total_sites)} loci | "
                f"{format_reads(int(speed))} loci/s | ETA {format_duration(eta)}",
            )
            chunk_index += 1

        # The target TSV should collapse to exactly the BED's unique coordinate count.
        try:
            extra = next(loci)
        except StopIteration:
            extra = None
        if extra is not None:
            raise RuntimeError(
                f"Target TSV contains additional unique locus after BED exhausted: {extra.chromosome}:{extra.position}"
            )

    if completed_sites != total_sites:
        raise RuntimeError(f"Stage 16 completed {completed_sites:,} of {total_sites:,} target loci")

    stage16_progress(92.0, f"ALL CALL CHECKPOINTS COMPLETE | {len(chunks)} chunks | concatenate VCFs")
    finalized: dict[str, Path] = {}

    def attempt_finalize(_attempt: int) -> None:
        part, part_tbi = finalize_calls(
            root=root, checkpoint_dir=checkpoint_dir, output=output, image=args.bcftools_image, chunks=chunks
        )
        finalized["vcf"] = part
        finalized["tbi"] = part_tbi

    try:
        retry_transient(root, "Stage 16 final concat/index", attempt_finalize)
        part = finalized["vcf"]
        part_tbi = finalized["tbi"]
        stage16_progress(96.0, "SEMANTIC QC | verify native ALT GQ + PL-derived reference quality + sex-aware ploidy")
        call_qc = validate_call_contract(part, resolved_sex)
        os.replace(part, output)
        os.replace(part_tbi, Path(str(output) + ".tbi"))
    finally:
        cleanup_stage16_containers(root)

    metadata = {
        "status": "complete",
        "schema": CALLING_SCHEMA,
        "calling_model": "bcftools-m-C-alleles-native-GQ-plus-PL-refQ-sex-aware-v5",
        "reference_confidence_model": "hom-ref: native GQ if present, otherwise PL-derived best-nonref likelihood gap; ALT-carried: native bcftools GQ",
        "reference_coordinate_policy": "skip-non-acgt-out-of-range-missing-contig-and-report-v2",
        "sample": {"name": sample_name, **sex_info},
        "known_alleles": {
            "panel_loci": total_panel_sites,
            "callable_target_loci": total_sites,
            "skipped_uncallable_reference_loci": skipped_uncallable,
            "skipped_non_acgt_reference_loci": skipped_non_acgt,
            "skipped_outside_fasta_contig_loci": skipped_outside,
            "skipped_missing_fasta_contig_loci": skipped_missing_contig,
            "uncallable_reason_counts": reference_contract["uncallable_reason_counts"],
            "uncallable_reference_symbols": reference_contract["uncallable_reference_symbols"],
            "uncallable_report": uncallable_output.name,
            "fallback_all_snp_loci_this_run": total_fallback,
            "target_ref_conflict_rows_this_run": total_ref_conflicts,
        },
        "call_qc": call_qc,
    }
    atomic_json(metadata_output, metadata)
    manifest["status"] = "complete_waiting_for_outer_state_validation"
    atomic_json(manifest_path, manifest)
    stage16_progress(
        100.0,
        f"DONE | known-alleles VCF finalized | hom-ref reference-quality coverage "
        f"{float(call_qc['hom_ref_reference_quality_fraction']):.2%} | ALT native-GQ "
        f"{float(call_qc['alt_carried_gq_fraction']):.2%} | checkpoints retained until state validation",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
