#!/usr/bin/env python3
"""Build one allele-level Dante catalog enriched with current ClinVar GRCh37.

The input VEP summaries are deterministic local annotations. ClinVar is joined
by the stable GRCh37 chromosome/position/REF/ALT key. The script never uploads
genomic data and keeps cache-derived and current-ClinVar fields separate.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import sqlite3
import sys
import tempfile
import time
import traceback
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO
from urllib.parse import unquote


CURRENT_CLINVAR_COLUMNS = [
    "clinvar_current_variation_id",
    "clinvar_current_rsid",
    "clinvar_current_significance",
    "clinvar_current_review_status",
    "clinvar_review_stars",
    "clinvar_current_conditions",
    "clinvar_current_disease_ids",
    "clinvar_current_origin",
    "clinvar_current_variant_type",
    "clinvar_current_geneinfo",
    "clinvar_exact_match",
]


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


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block_size):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "name": path.name,
        "size": stat.st_size,
        "sha256": sha256_file(path),
    }


def parse_info(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    if not value or value == ".":
        return result
    for token in value.split(";"):
        if "=" in token:
            key, val = token.split("=", 1)
            result[key] = unquote(val)
        elif token:
            result[token] = "true"
    return result


def normalize_chrom(value: str) -> str:
    chrom = value.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom


def variant_key(chrom: str, pos: str, ref: str, alt: str) -> str:
    return f"GRCh37:{normalize_chrom(chrom)}:{pos}:{ref}:{alt}"


def review_stars(status: str) -> int:
    normalized = status.lower().replace("_", " ")
    if "practice guideline" in normalized:
        return 4
    if "reviewed by expert panel" in normalized:
        return 3
    if "multiple submitters" in normalized and "no conflicts" in normalized:
        return 2
    if "criteria provided" in normalized:
        return 1
    return 0


def clinvar_file_date(path: Path) -> str:
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("##fileDate="):
                return raw.strip().split("=", 1)[1]
            if raw.startswith("#CHROM"):
                break
    return "unknown"


def create_database(parent: Path) -> tuple[sqlite3.Connection, Path]:
    fd, name = tempfile.mkstemp(prefix="dante_clinvar_", suffix=".sqlite3", dir=parent)
    os.close(fd)
    path = Path(name)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    conn.execute(
        """CREATE TABLE clinvar (
            variant_key TEXT PRIMARY KEY,
            variation_id TEXT,
            rsid TEXT,
            significance TEXT,
            review_status TEXT,
            conditions TEXT,
            disease_ids TEXT,
            origin TEXT,
            variant_type TEXT,
            geneinfo TEXT
        ) WITHOUT ROWID"""
    )
    conn.execute("CREATE TABLE seen (variant_key TEXT PRIMARY KEY) WITHOUT ROWID")
    return conn, path


def index_clinvar(conn: sqlite3.Connection, path: Path, progress_every: int) -> Counter[str]:
    counts: Counter[str] = Counter()
    with open_text(path) as handle:
        conn.execute("BEGIN")
        for raw in handle:
            if not raw or raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                counts["malformed"] += 1
                continue
            chrom, pos, ident, ref, alt_field = fields[:5]
            info = parse_info(fields[7])
            for alt in alt_field.split(","):
                key = variant_key(chrom, pos, ref, alt)
                values = (
                    key,
                    ident,
                    info.get("RS", "."),
                    info.get("CLNSIG", "."),
                    info.get("CLNREVSTAT", "."),
                    info.get("CLNDN", "."),
                    info.get("CLNDISDB", "."),
                    info.get("CLNORIGIN", "."),
                    info.get("CLNVC", info.get("CLNVCSO", ".")),
                    info.get("GENEINFO", "."),
                )
                cursor = conn.execute(
                    "INSERT OR IGNORE INTO clinvar VALUES (?,?,?,?,?,?,?,?,?,?)", values
                )
                if cursor.rowcount == 0:
                    counts["duplicate_keys"] += 1
                else:
                    counts["alleles_indexed"] += 1
            counts["records"] += 1
            if progress_every and counts["records"] % progress_every == 0:
                print(f"  ClinVar: {counts['records']:,} records indexed", flush=True)
                conn.commit()
                conn.execute("BEGIN")
        conn.commit()
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a unified Dante VEP + current ClinVar catalog")
    parser.add_argument("--summary", type=Path, action="append", required=True)
    parser.add_argument("--clinvar", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--vep-version", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--progress-every", type=int, default=500_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sources = [path.resolve() for path in args.summary]
    clinvar = args.clinvar.resolve()
    for path in [*sources, clinvar]:
        if not path.is_file():
            raise FileNotFoundError(path)
    for target in (args.output, args.report, args.manifest):
        target.parent.mkdir(parents=True, exist_ok=True)

    print("Calculating SHA256 input identities...", flush=True)
    inputs = {
        "vep_summaries": [file_identity(path) for path in sources],
        "clinvar": file_identity(clinvar),
    }
    if args.manifest.exists() and args.output.exists():
        try:
            previous = json.loads(args.manifest.read_text(encoding="utf-8"))
            if previous.get("inputs") == inputs and previous.get("catalog", {}).get("status") == "complete":
                print("Unified catalog is current for the same SHA256 inputs; reusing it.")
                return 0
        except (OSError, ValueError, TypeError):
            pass

    started = time.time()
    conn, db_path = create_database(args.output.parent)
    output_part = args.output.with_name(args.output.name + ".part")
    counts: Counter[str] = Counter()
    layer_counts: Counter[str] = Counter()
    clinvar_counts: Counter[str] = Counter()
    header: list[str] | None = None
    try:
        print("Indexing current ClinVar GRCh37...", flush=True)
        clinvar_counts = index_clinvar(conn, clinvar, args.progress_every)
        with open_anonymous_gzip_text(output_part) as output_handle:
            writer: csv.DictWriter | None = None
            for source in sources:
                print(f"Merging {source.name}...", flush=True)
                with open_text(source) as input_handle:
                    reader = csv.DictReader(input_handle, delimiter="\t")
                    if reader.fieldnames is None:
                        raise RuntimeError(f"No TSV header found in {source}")
                    if header is None:
                        header = list(reader.fieldnames)
                        writer = csv.DictWriter(
                            output_handle,
                            fieldnames=[*header, *CURRENT_CLINVAR_COLUMNS],
                            delimiter="\t",
                            lineterminator="\n",
                            extrasaction="ignore",
                        )
                        writer.writeheader()
                    elif list(reader.fieldnames) != header:
                        raise RuntimeError(f"VEP summary columns differ in {source.name}")
                    assert writer is not None
                    for row in reader:
                        counts["summary_allele_rows"] += 1
                        if row.get("alt_carried") != "yes":
                            counts["uncarried_alt_rows_skipped"] += 1
                            continue
                        key = row.get("variant_key") or variant_key(
                            row.get("chrom", ""), row.get("position", ""), row.get("ref", ""), row.get("alt", "")
                        )
                        row["variant_key"] = key
                        cursor = conn.execute("INSERT OR IGNORE INTO seen VALUES (?)", (key,))
                        if cursor.rowcount == 0:
                            counts["duplicate_carried_alleles_skipped"] += 1
                            continue
                        match = conn.execute(
                            "SELECT variation_id,rsid,significance,review_status,conditions,disease_ids,origin,variant_type,geneinfo "
                            "FROM clinvar WHERE variant_key=?",
                            (key,),
                        ).fetchone()
                        if match:
                            review = match[3] or "."
                            enrichment = {
                                "clinvar_current_variation_id": match[0] or ".",
                                "clinvar_current_rsid": f"rs{match[1]}" if str(match[1]).isdigit() else (match[1] or "."),
                                "clinvar_current_significance": match[2] or ".",
                                "clinvar_current_review_status": review,
                                "clinvar_review_stars": str(review_stars(review)),
                                "clinvar_current_conditions": match[4] or ".",
                                "clinvar_current_disease_ids": match[5] or ".",
                                "clinvar_current_origin": match[6] or ".",
                                "clinvar_current_variant_type": match[7] or ".",
                                "clinvar_current_geneinfo": match[8] or ".",
                                "clinvar_exact_match": "yes",
                            }
                            counts["current_clinvar_exact_matches"] += 1
                        else:
                            enrichment = {column: "." for column in CURRENT_CLINVAR_COLUMNS}
                            enrichment["clinvar_review_stars"] = "0"
                            enrichment["clinvar_exact_match"] = "no"
                        writer.writerow({**row, **enrichment})
                        counts["catalog_rows"] += 1
                        layer_counts[row.get("source_layer", "unknown")] += 1
                        if args.progress_every and counts["catalog_rows"] % args.progress_every == 0:
                            print(f"  Catalog: {counts['catalog_rows']:,} carried alleles written", flush=True)
        os.replace(output_part, args.output)
    except Exception:
        output_part.unlink(missing_ok=True)
        raise
    finally:
        conn.close()
        db_path.unlink(missing_ok=True)

    elapsed = time.time() - started
    file_date = clinvar_file_date(clinvar)
    report_lines = [
        "DANTE UNIFIED SMALL-VARIANT CATALOG REPORT",
        "=" * 48,
        f"Assembly: {args.assembly}",
        f"VEP cache/tool version: {args.vep_version}",
        f"ClinVar file date: {file_date}",
        f"Summary allele rows read: {counts['summary_allele_rows']:,}",
        f"Uncarried ALT rows skipped: {counts['uncarried_alt_rows_skipped']:,}",
        f"Duplicate carried alleles skipped: {counts['duplicate_carried_alleles_skipped']:,}",
        f"Carried allele catalog rows: {counts['catalog_rows']:,}",
        f"Current ClinVar exact matches: {counts['current_clinvar_exact_matches']:,}",
        f"ClinVar records indexed: {clinvar_counts['records']:,}",
        f"ClinVar alleles indexed: {clinvar_counts['alleles_indexed']:,}",
        "",
        "ROWS BY SOURCE LAYER",
        *[f"{key}: {value:,}" for key, value in sorted(layer_counts.items())],
        "",
        "NOTES",
        "- One row is retained per actually carried ALT allele.",
        "- VEP-cache clinical fields and current ClinVar fields are kept separate.",
        "- ClinVar is joined by exact normalized GRCh37 chromosome/position/REF/ALT.",
        "- Review stars summarize ClinVar review status; they do not replace variant interpretation.",
        f"- Elapsed seconds: {elapsed:.1f}",
        "",
    ]
    args.report.write_text("\n".join(report_lines), encoding="utf-8")
    manifest = {
        "schema_version": 1,
        "pipeline": "Dante WGS local pipeline v7.5.2",
        "assembly": args.assembly,
        "vep_version": args.vep_version,
        "clinvar_file_date": file_date,
        "inputs": inputs,
        "catalog": {
            "status": "complete",
            "file": args.output.name,
            "sha256": sha256_file(args.output),
            "rows": counts["catalog_rows"],
        },
        "counts": dict(counts),
    }
    args.manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Created {args.output}: {counts['catalog_rows']:,} carried alleles")
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
