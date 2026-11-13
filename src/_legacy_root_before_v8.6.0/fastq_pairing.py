#!/usr/bin/env python3
"""Discover and validate one single-end R1 or one paired-end FASTQ input."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import time
import os
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


FASTQ_SUFFIX = re.compile(r"(?i)\.(?:fastq|fq)(?:\.gz)?$")
PAIR_PATTERNS = (
    re.compile(
        r"(?i)^(?P<prefix>.+?)(?P<sep>[_\.-])R(?P<read>[12])"
        r"(?P<tail>(?:[_\.-].*)?)(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$"
    ),
    re.compile(
        r"(?i)^(?P<prefix>.+?)(?P<sep>[_\.-])(?P<read>[12])"
        r"(?P<tail>(?:[_\.-].*)?)(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$"
    ),
)


@dataclass(frozen=True)
class Candidate:
    path: Path
    read: int
    key: tuple[str, str, str, str]


def classify(path: Path) -> Candidate | None:
    if not path.is_file() or not FASTQ_SUFFIX.search(path.name):
        return None
    for index, pattern in enumerate(PAIR_PATTERNS):
        match = pattern.match(path.name)
        if match:
            key = (
                str(index),
                match.group("prefix").lower(),
                match.group("tail").lower(),
                match.group("ext").lower(),
            )
            return Candidate(path=path, read=int(match.group("read")), key=key)
    return None


def open_fastq(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="strict", newline="")
    return path.open("rt", encoding="ascii", errors="strict", newline="")


def normalized_read_name(header: str) -> str:
    token = header[1:].strip().split()[0]
    return re.sub(r"/[12]$", "", token)


def validate_fastq(path: Path, expected_read: int, records: int = 16) -> list[str]:
    errors: list[str] = []
    seen = 0
    try:
        with open_fastq(path) as handle:
            while seen < records:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().rstrip("\r\n")
                plus = handle.readline()
                quality = handle.readline().rstrip("\r\n")
                seen += 1
                if not header.startswith("@"):
                    errors.append(f"record {seen}: header does not start with @")
                    break
                if not plus.startswith("+"):
                    errors.append(f"record {seen}: separator does not start with +")
                    break
                if not sequence or len(sequence) != len(quality):
                    errors.append(f"record {seen}: sequence/quality lengths differ")
                    break
                if not re.fullmatch(r"[A-Za-z.*=-]+", sequence):
                    errors.append(f"record {seen}: sequence contains invalid characters")
                    break
                fields = header[1:].strip().split()
                if len(fields) > 1 and re.match(r"^[12]:", fields[1]):
                    observed = int(fields[1][0])
                    if observed != expected_read:
                        errors.append(
                            f"record {seen}: header declares read {observed}, expected {expected_read}"
                        )
                        break
    except (OSError, EOFError, UnicodeError) as exc:
        errors.append(f"cannot read FASTQ: {exc}")
    if seen == 0 and not errors:
        errors.append("FASTQ is empty")
    return errors


def validate_pair_names(r1: Path, r2: Path, records: int = 16) -> list[str]:
    errors: list[str] = []
    try:
        with open_fastq(r1) as first, open_fastq(r2) as second:
            for index in range(1, records + 1):
                h1 = first.readline()
                h2 = second.readline()
                if not h1 and not h2:
                    break
                if not h1 or not h2:
                    errors.append(f"record {index}: one mate ended before the other")
                    break
                if normalized_read_name(h1) != normalized_read_name(h2):
                    errors.append(
                        f"record {index}: mate names differ ({normalized_read_name(h1)} != "
                        f"{normalized_read_name(h2)})"
                    )
                    break
                for handle in (first, second):
                    for _ in range(3):
                        if not handle.readline():
                            errors.append(f"record {index}: truncated FASTQ record")
                            return errors
    except (OSError, EOFError, UnicodeError) as exc:
        errors.append(f"cannot compare mate names: {exc}")
    return errors


def full_stream_check(
    path: Path,
    *,
    overall_start: float = 0.0,
    overall_span: float = 5.0,
    progress_stage: int = 15,
) -> tuple[str | None, str | None, dict[str, int] | None]:
    """Validate the complete stream, hash decompressed bytes, and count FASTQ records.

    The exact record count is persisted so downstream stages can report progress as
    processed reads / total reads instead of guessing from output size.
    """
    digest = hashlib.sha256()
    try:
        source_size = path.stat().st_size
        processed = 0
        newline_count = 0
        last_byte: int | None = None
        next_progress = time.monotonic()
        handle_context = gzip.open(path, "rb") if path.name.lower().endswith(".gz") else path.open("rb")
        with handle_context as handle:
            while chunk := handle.read(8 * 1024 * 1024):
                digest.update(chunk)
                processed += len(chunk)
                newline_count += chunk.count(b"\n")
                last_byte = chunk[-1]
                if source_size >= 1024 * 1024 and time.monotonic() >= next_progress:
                    if path.name.lower().endswith(".gz"):
                        source_position = handle.fileobj.tell()  # type: ignore[attr-defined]
                    else:
                        source_position = handle.tell()
                    percent = 100.0 * source_position / source_size if source_size else 100.0
                    percent = min(100.0, percent)
                    overall = overall_start + overall_span * percent / 100.0
                    print(
                        f"FASTQ integrity {path.name}: {percent:.2f}% | "
                        f"{processed / (1024 ** 3):.2f} GiB decompressed",
                        file=sys.stderr,
                        flush=True,
                    )
                    print(
                        f"[OVERALL] Stage {progress_stage} | {overall:.2f}% | VALIDATION {percent:.2f}% | {path.name}",
                        file=sys.stderr,
                        flush=True,
                    )
                    next_progress = time.monotonic() + 5.0
        logical_lines = newline_count + (1 if processed and last_byte != 0x0A else 0)
        if logical_lines == 0:
            return "full stream integrity check failed: FASTQ is empty", None, None
        if logical_lines % 4 != 0:
            return (
                f"full stream integrity check failed: {logical_lines:,} logical lines is not divisible by 4",
                None,
                None,
            )
        records = logical_lines // 4
        if source_size >= 1024 * 1024:
            print(
                f"FASTQ integrity {path.name}: 100.00% | {processed / (1024 ** 3):.2f} GiB decompressed | "
                f"{records:,} records",
                file=sys.stderr,
                flush=True,
            )
            print(
                f"[OVERALL] Stage {progress_stage} | {overall_start + overall_span:.2f}% | VALIDATION 100.00% | "
                f"{path.name} | {records:,} records",
                file=sys.stderr,
                flush=True,
            )
    except (OSError, EOFError) as exc:
        return f"full stream integrity check failed: {exc}", None, None
    return None, digest.hexdigest(), {
        "record_count": records,
        "decompressed_bytes": processed,
        "logical_lines": logical_lines,
    }


def discover(root: Path, do_full_gzip_check: bool = False, progress_stage: int = 15, progress_start: float = 0.0, progress_span: float = 5.0) -> dict[str, object]:
    candidates = [item for item in (classify(path) for path in root.iterdir()) if item]
    grouped: dict[tuple[str, str, str, str], dict[int, list[Path]]] = {}
    for candidate in candidates:
        grouped.setdefault(candidate.key, {1: [], 2: []})[candidate.read].append(candidate.path)

    complete: list[tuple[Path, Path]] = []
    incomplete: list[str] = []
    for reads in grouped.values():
        r1s = sorted(reads[1])
        r2s = sorted(reads[2])
        if len(r1s) == 1 and len(r2s) == 1:
            complete.append((r1s[0], r2s[0]))
        else:
            names = [path.name for path in r1s + r2s]
            incomplete.append(", ".join(names))

    result: dict[str, object] = {
        "status": "MISSING",
        "read_mode": None,
        "r1": None,
        "r2": None,
        "r1_content_sha256": None,
        "r2_content_sha256": None,
        "r1_record_count": None,
        "r2_record_count": None,
        "r1_decompressed_bytes": None,
        "r2_decompressed_bytes": None,
        "candidate_count": len(candidates),
        "complete_pair_count": len(complete),
        "details": [],
    }
    details: list[str] = result["details"]  # type: ignore[assignment]

    if len(complete) > 1:
        result["status"] = "MULTIPLE_PAIRS"
        details.append("More than one complete FASTQ pair was found; lane merging is not implicit.")
        details.extend(f"{r1.name} + {r2.name}" for r1, r2 in complete)
        return result
    if len(complete) == 0:
        r1s = sorted(candidate.path for candidate in candidates if candidate.read == 1)
        r2s = sorted(candidate.path for candidate in candidates if candidate.read == 2)
        if len(r1s) == 1 and not r2s:
            r1 = r1s[0]
            result["read_mode"] = "single"
            result["r1"] = str(r1s[0].resolve())
            errors = validate_fastq(r1, 1)
            if do_full_gzip_check:
                error, digest, stats = full_stream_check(r1, overall_start=progress_start, overall_span=progress_span, progress_stage=progress_stage)
                if error:
                    errors.append(f"{r1.name}: {error}")
                else:
                    result["r1_content_sha256"] = digest
                    assert stats is not None
                    result["r1_record_count"] = stats["record_count"]
                    result["r1_decompressed_bytes"] = stats["decompressed_bytes"]
            if errors:
                result["status"] = "INVALID"
                details.extend(errors)
            else:
                result["status"] = "SINGLE_END_READY"
                details.append(f"Single-end R1: {r1.name}; matching R2 is not present.")
        elif len(r2s) == 1 and not r1s:
            result["status"] = "R1_MISSING"
            result["r2"] = str(r2s[0].resolve())
            details.append(f"R2 found: {r2s[0].name}; matching R1 is missing.")
        elif candidates:
            result["status"] = "PAIRING_ERROR"
            details.append("FASTQ files were found but could not be paired unambiguously.")
            details.extend(incomplete)
        return result

    r1, r2 = complete[0]
    result["read_mode"] = "paired"
    result["r1"] = str(r1.resolve())
    result["r2"] = str(r2.resolve())
    errors = validate_fastq(r1, 1) + validate_fastq(r2, 2) + validate_pair_names(r1, r2)
    if do_full_gzip_check:
        pair_stats: dict[int, dict[str, int]] = {}
        for index, (path, digest_name, record_name, bytes_name) in enumerate(
            (
                (r1, "r1_content_sha256", "r1_record_count", "r1_decompressed_bytes"),
                (r2, "r2_content_sha256", "r2_record_count", "r2_decompressed_bytes"),
            )
        ):
            half_span = progress_span / 2.0
            error, digest, stats = full_stream_check(
                path,
                overall_start=progress_start + half_span * index,
                overall_span=half_span,
                progress_stage=progress_stage,
            )
            if error:
                errors.append(f"{path.name}: {error}")
            else:
                result[digest_name] = digest
                assert stats is not None
                result[record_name] = stats["record_count"]
                result[bytes_name] = stats["decompressed_bytes"]
                pair_stats[index + 1] = stats
        if 1 in pair_stats and 2 in pair_stats and pair_stats[1]["record_count"] != pair_stats[2]["record_count"]:
            errors.append(
                f"paired FASTQ record counts differ: R1={pair_stats[1]['record_count']:,}, "
                f"R2={pair_stats[2]['record_count']:,}"
            )
    if errors:
        result["status"] = "INVALID"
        details.extend(errors)
    else:
        result["status"] = "READY"
        details.append(f"Paired FASTQ: {r1.name} + {r2.name}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--full-gzip-check", action="store_true")
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--progress-stage", type=int, default=15)
    parser.add_argument("--progress-start", type=float, default=0.0)
    parser.add_argument("--progress-span", type=float, default=5.0)
    args = parser.parse_args()
    result = discover(
        args.root.resolve(),
        args.full_gzip_check,
        progress_stage=args.progress_stage,
        progress_start=args.progress_start,
        progress_span=args.progress_span,
    )
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        part = Path(str(args.output_json) + ".part")
        part.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(part, args.output_json)
    else:
        json.dump(result, sys.stdout, ensure_ascii=False)
        sys.stdout.write("\n")
    if args.require_complete and result["status"] not in {"READY", "SINGLE_END_READY"}:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
