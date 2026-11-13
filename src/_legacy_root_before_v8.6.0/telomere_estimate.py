#!/usr/bin/env python3
"""Estimate mean telomere length from Dante Illumina WGS FASTQ.

The primary estimate is a transparent TelSeq-compatible read-count estimator
implemented directly on FASTQ so Stage 4 can run on R1-only data without
requiring a whole-genome BAM.  It counts canonical TTAGGG/CCCTAA repeats and
normalizes by reads in the 48-52% GC window using the published TelSeq human
constants.  Because raw FASTQ does not carry duplicate flags, the report makes
that difference from the original BAM-based TelSeq explicit.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import sys
import time
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, Optional

SCHEMA = "dante-telomere-metrics-v1"
METHOD = "Dante TelSeq-compatible FASTQ estimator v1"
MOTIF_FWD = b"TTAGGG"
MOTIF_REV = b"CCCTAA"
DEFAULT_K = 7
GENOME_LENGTH_AT_TEL_GC = 332_720_800
TELOMERE_ENDS = 46
LENGTH_UNIT = 1000  # output in kb
GC_LOW = 0.48
GC_HIGH = 0.52
SENSITIVITY_K = tuple(range(5, 11))


def open_fastq(path: Path) -> BinaryIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rb")
    return path.open("rb")


def human_duration(seconds: float) -> str:
    if not math.isfinite(seconds) or seconds < 0:
        return "unknown"
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = Path(str(path) + ".part")
    part.write_text(text, encoding="utf-8", newline="\n")
    os.replace(part, path)


def atomic_json(path: Path, payload: dict) -> None:
    atomic_text(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def telomeric_motif_count(seq: bytes) -> int:
    # bytes.count is non-overlapping, matching the original TelSeq motif-count logic.
    return max(seq.count(MOTIF_FWD), seq.count(MOTIF_REV))


def empty_counts() -> Dict[str, object]:
    return {
        "records": 0,
        "bases": 0,
        "gc_normalizer_reads": 0,
        "telomeric_reads_k7": 0,
        "threshold_counts": {str(k): 0 for k in SENSITIVITY_K},
        "min_read_length": None,
        "max_read_length": 0,
        "first_header": None,
    }


def estimate_kb(tel_reads: int, gc_reads: int) -> Optional[float]:
    if gc_reads <= 0:
        return None
    return (tel_reads / gc_reads) * GENOME_LENGTH_AT_TEL_GC / LENGTH_UNIT / TELOMERE_ENDS


def scan_fastq(path: Path, expected_records: int, already_processed: int, total_expected: int, label: str, overall_start: float = 5.0, overall_span: float = 95.0) -> dict:
    counts = empty_counts()
    started = time.monotonic()
    last_print = 0.0
    with open_fastq(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not seq or not plus or not qual:
                raise RuntimeError(f"Truncated FASTQ record in {path.name} near record {int(counts['records']) + 1:,}")
            if not header.startswith(b"@") or not plus.startswith(b"+"):
                raise RuntimeError(f"Invalid FASTQ structure in {path.name} near record {int(counts['records']) + 1:,}")
            seq = seq.strip().upper()
            qual = qual.rstrip(b"\r\n")
            if len(seq) != len(qual):
                raise RuntimeError(f"Sequence/quality length mismatch in {path.name} near record {int(counts['records']) + 1:,}")
            if counts["first_header"] is None:
                counts["first_header"] = header.decode("ascii", "replace").strip()
            length = len(seq)
            counts["records"] = int(counts["records"]) + 1
            counts["bases"] = int(counts["bases"]) + length
            counts["min_read_length"] = length if counts["min_read_length"] is None else min(int(counts["min_read_length"]), length)
            counts["max_read_length"] = max(int(counts["max_read_length"]), length)

            if length:
                gc = (seq.count(b"G") + seq.count(b"C")) / length
                if GC_LOW <= gc < GC_HIGH:
                    counts["gc_normalizer_reads"] = int(counts["gc_normalizer_reads"]) + 1
            motifs = telomeric_motif_count(seq)
            thresholds = counts["threshold_counts"]
            assert isinstance(thresholds, dict)
            for k in SENSITIVITY_K:
                if motifs >= k:
                    thresholds[str(k)] = int(thresholds[str(k)]) + 1
            if motifs >= DEFAULT_K:
                counts["telomeric_reads_k7"] = int(counts["telomeric_reads_k7"]) + 1

            now = time.monotonic()
            if now - last_print >= 5.0:
                processed = already_processed + int(counts["records"])
                scan_pct = (100.0 * processed / total_expected) if total_expected > 0 else 0.0
                pct = overall_start + overall_span * scan_pct / 100.0
                elapsed = now - started
                local_rate = int(counts["records"]) / elapsed if elapsed > 0 else 0.0
                remaining = max(0, total_expected - processed)
                eta = human_duration(remaining / local_rate) if local_rate > 0 else "unknown"
                print(
                    f"[OVERALL] Stage 4 | {pct:6.2f}% | TELOMERE SCAN {label} | "
                    f"{processed/1_000_000:.1f}M / {total_expected/1_000_000:.1f}M reads | "
                    f"{local_rate/1000:.1f}k reads/s | ETA {eta}",
                    flush=True,
                )
                last_print = now
    if expected_records > 0 and int(counts["records"]) != expected_records:
        raise RuntimeError(
            f"FASTQ record count changed for {path.name}: expected {expected_records:,}, observed {int(counts['records']):,}"
        )
    counts["mean_read_length"] = (int(counts["bases"]) / int(counts["records"])) if int(counts["records"]) else 0.0
    counts["estimate_kb_k7"] = estimate_kb(int(counts["telomeric_reads_k7"]), int(counts["gc_normalizer_reads"]))
    counts["sensitivity_estimate_kb"] = {
        str(k): estimate_kb(int(counts["threshold_counts"][str(k)]), int(counts["gc_normalizer_reads"]))
        for k in SENSITIVITY_K
    }
    return counts


def combined_metrics(parts: Iterable[dict]) -> dict:
    parts = list(parts)
    total = empty_counts()
    total["min_read_length"] = None
    for part in parts:
        total["records"] = int(total["records"]) + int(part["records"])
        total["bases"] = int(total["bases"]) + int(part["bases"])
        total["gc_normalizer_reads"] = int(total["gc_normalizer_reads"]) + int(part["gc_normalizer_reads"])
        total["telomeric_reads_k7"] = int(total["telomeric_reads_k7"]) + int(part["telomeric_reads_k7"])
        for k in SENSITIVITY_K:
            total["threshold_counts"][str(k)] = int(total["threshold_counts"][str(k)]) + int(part["threshold_counts"][str(k)])
        if part.get("min_read_length") is not None:
            total["min_read_length"] = int(part["min_read_length"]) if total["min_read_length"] is None else min(int(total["min_read_length"]), int(part["min_read_length"]))
        total["max_read_length"] = max(int(total["max_read_length"]), int(part["max_read_length"]))
    total["mean_read_length"] = (int(total["bases"]) / int(total["records"])) if int(total["records"]) else 0.0
    total["estimate_kb_k7"] = estimate_kb(int(total["telomeric_reads_k7"]), int(total["gc_normalizer_reads"]))
    total["sensitivity_estimate_kb"] = {
        str(k): estimate_kb(int(total["threshold_counts"][str(k)]), int(total["gc_normalizer_reads"]))
        for k in SENSITIVITY_K
    }
    total.pop("first_header", None)
    return total


def read_vcf_metadata(root: Path) -> list[str]:
    clues: list[str] = []
    for path in sorted(root.glob("*.vcf*")):
        if not path.is_file():
            continue
        try:
            opener = gzip.open if path.name.lower().endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[arg-type]
                for _ in range(250):
                    line = handle.readline()
                    if not line or line.startswith("#CHROM"):
                        break
                    if re.match(r"^##(?:fileDate|date|runDate|analysisDate)=", line, flags=re.I):
                        clues.append(f"{path.name}: {line.strip()}")
        except Exception as exc:
            clues.append(f"{path.name}: metadata header unreadable ({exc.__class__.__name__})")
    return clues


def format_estimate(value: Optional[float]) -> str:
    return "UNKNOWN" if value is None else f"{value:.3f} kb"


def render_report(payload: dict, sample_date: str, sample_year: str, sample_material: str, metadata_root: Path) -> str:
    combined = payload["combined"]
    lines: list[str] = []
    lines.append("DANTE WGS - TELOMERE LENGTH ESTIMATE")
    lines.append("=" * 72)
    lines.append("")
    lines.append("PRIMARY RESULT")
    lines.append(f"Estimated mean telomere length: {format_estimate(combined.get('estimate_kb_k7'))}")
    lines.append("Estimator: TelSeq-compatible canonical-repeat count from raw FASTQ, k>=7")
    lines.append(f"Read mode: {payload['read_mode'].upper()}")
    lines.append("")
    lines.append("SAMPLE / DATE")
    if sample_date:
        lines.append(f"Sample date (user-supplied): {sample_date}")
    elif sample_year:
        lines.append(f"Sample year (user-supplied): {sample_year}")
    else:
        lines.append("Sample date: UNKNOWN")
        lines.append("The standard Illumina FASTQ read identifier does not encode a reliable specimen-collection date.")
    lines.append(f"Sample material: {sample_material if sample_material else 'UNKNOWN'}")
    lines.append("")
    lines.append("INPUTS")
    for key in ("r1", "r2"):
        item = payload.get(key)
        if not item:
            continue
        lines.append(f"{key.upper()}: {item['name']}")
        lines.append(f"  SHA256 (validated decompressed FASTQ stream): {item.get('content_sha256') or 'UNKNOWN'}")
        lines.append(f"  records: {int(item['metrics']['records']):,}")
        lines.append(f"  bases: {int(item['metrics']['bases']):,}")
        lines.append(f"  mean read length: {float(item['metrics']['mean_read_length']):.2f} bp")
        lines.append(f"  first FASTQ header: {item['metrics'].get('first_header') or 'UNKNOWN'}")
        lines.append(f"  individual k=7 estimate: {format_estimate(item['metrics'].get('estimate_kb_k7'))}")
    lines.append("")
    lines.append("COUNTING / NORMALIZATION QC")
    lines.append(f"Total reads scanned: {int(combined['records']):,}")
    lines.append(f"Total sequenced bases scanned: {int(combined['bases']):,}")
    lines.append(f"GC-normalizer reads (48-52% GC): {int(combined['gc_normalizer_reads']):,}")
    lines.append(f"Canonical telomeric reads (>=7 non-overlapping TTAGGG or CCCTAA motifs): {int(combined['telomeric_reads_k7']):,}")
    if int(combined["gc_normalizer_reads"]) > 0:
        lines.append(f"Telomeric/GC-normalizer ratio: {int(combined['telomeric_reads_k7']) / int(combined['gc_normalizer_reads']):.8f}")
    lines.append("")
    lines.append("THRESHOLD SENSITIVITY (QC ONLY; k=7 is the primary TelSeq threshold)")
    for k in SENSITIVITY_K:
        count = int(combined["threshold_counts"][str(k)])
        est = combined["sensitivity_estimate_kb"][str(k)]
        lines.append(f"k>={k}: {count:,} reads; calibrated-formula value {format_estimate(est)}")
    lines.append("")
    lines.append("METHOD")
    lines.append("Formula: (telomeric reads / 48-52% GC reads) x 332,720,800 / 1000 / 46.")
    lines.append("These constants and the default k=7 canonical-repeat threshold reproduce the human TelSeq calculation contract.")
    lines.append("This implementation scans raw FASTQ directly. Unlike original BAM-based TelSeq, raw FASTQ has no duplicate flags, so PCR/optical duplicates cannot be excluded here.")
    lines.append("This is therefore a TelSeq-compatible historical WGS estimate, not a laboratory telomere assay and not a claim of chromosome-specific telomere length.")
    lines.append("Do not convert this value directly into a biological-age number. Longitudinal comparison is most useful when tissue/source and sequencing workflow are comparable.")
    lines.append("")
    lines.append("REFERENCES")
    lines.append("Ding Z et al. Estimating telomere length from whole genome sequence data. Nucleic Acids Res. 2014;42:e75. DOI: 10.1093/nar/gku181")
    lines.append("TelSeq source constants: https://github.com/zd1/telseq")
    lines.append("Nersisyan L, Arakelyan A. Computel. PLoS ONE. 2015;10:e0125201. DOI: 10.1371/journal.pone.0125201")
    lines.append("")
    lines.append("DATE / METADATA CLUES (NOT AUTOMATICALLY TREATED AS SAMPLE DATE)")
    clues = read_vcf_metadata(metadata_root)
    if clues:
        lines.extend(f"- {clue}" for clue in clues)
    else:
        lines.append("- No explicit date-like VCF header fields were found in the root VCF files.")
    for key in ("r1", "r2"):
        item = payload.get(key)
        if item:
            lines.append(f"- {key.upper()} filesystem mtime: {item.get('filesystem_mtime_utc', 'UNKNOWN')} (weak provenance only)")
    lines.append("")
    lines.append("REPRODUCIBILITY")
    lines.append(f"Metrics schema: {payload['schema']}")
    lines.append(f"Method implementation: {payload['method']}")
    lines.append(f"Generated UTC: {payload['generated_utc']}")
    lines.append("")
    return "\n".join(lines)


def cached_matches(cached: dict, read_mode: str, r1_sha: str, r2_sha: str) -> bool:
    if cached.get("schema") != SCHEMA or cached.get("read_mode") != read_mode:
        return False
    if (cached.get("r1") or {}).get("content_sha256") != r1_sha:
        return False
    if read_mode == "paired" and (cached.get("r2") or {}).get("content_sha256") != r2_sha:
        return False
    if read_mode == "single" and cached.get("r2"):
        return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path)
    parser.add_argument("--r1-sha256", required=True)
    parser.add_argument("--r2-sha256", default="")
    parser.add_argument("--r1-records", type=int, required=True)
    parser.add_argument("--r2-records", type=int, default=0)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--metadata-root", type=Path, required=True)
    parser.add_argument("--sample-date", nargs="?", const="", default="")
    parser.add_argument("--sample-year", nargs="?", const="", default="")
    parser.add_argument("--sample-material", nargs="?", const="", default="")
    args = parser.parse_args()

    read_mode = "paired" if args.r2 else "single"
    cached = None
    if args.metrics.exists():
        try:
            cached = json.loads(args.metrics.read_text(encoding="utf-8"))
        except Exception:
            cached = None
    if cached and cached_matches(cached, read_mode, args.r1_sha256, args.r2_sha256):
        payload = cached
        print("Telomere FASTQ metrics are current; reusing cached scan and regenerating the TXT report.", flush=True)
    else:
        total_expected = args.r1_records + (args.r2_records if args.r2 else 0)
        r1_metrics = scan_fastq(args.r1, args.r1_records, 0, total_expected, "R1")
        r1 = {
            "name": args.r1.name,
            "content_sha256": args.r1_sha256,
            "filesystem_mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(args.r1.stat().st_mtime)),
            "metrics": r1_metrics,
        }
        r2 = None
        if args.r2:
            r2_metrics = scan_fastq(args.r2, args.r2_records, args.r1_records, total_expected, "R2")
            r2 = {
                "name": args.r2.name,
                "content_sha256": args.r2_sha256,
                "filesystem_mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(args.r2.stat().st_mtime)),
                "metrics": r2_metrics,
            }
        payload = {
            "schema": SCHEMA,
            "method": METHOD,
            "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "read_mode": read_mode,
            "r1": r1,
            "r2": r2,
            "combined": combined_metrics([r1_metrics] + ([r2["metrics"]] if r2 else [])),
        }
        atomic_json(args.metrics, payload)
        print(f"Saved reusable telomere metrics: {args.metrics}", flush=True)

    report_text = render_report(payload, args.sample_date.strip(), args.sample_year.strip(), args.sample_material.strip(), args.metadata_root)
    atomic_text(args.report, report_text)
    print(f"Saved telomere report: {args.report}", flush=True)
    estimate = payload["combined"].get("estimate_kb_k7")
    if estimate is None:
        print("Primary telomere estimate: UNKNOWN", flush=True)
    else:
        print(f"Primary telomere estimate: {float(estimate):.3f} kb", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
