#!/usr/bin/env python3
"""Extract Dante raw small variants absent from filtered SNP and indel splits."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import os
import sqlite3
import sys
import tempfile
import time
import traceback
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


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


def iter_records(path: Path) -> Iterator[str]:
    with open_text(path) as handle:
        for raw in handle:
            if raw and not raw.startswith("#") and raw.strip():
                yield raw


def key_from_raw(raw: str) -> bytes:
    fields = raw.rstrip("\r\n").split("\t", 5)
    if len(fields) < 5:
        raise RuntimeError("VCF record has fewer than five columns")
    value = "\x1f".join((fields[0], fields[1], fields[3], fields[4])).encode("utf-8")
    return hashlib.blake2b(value, digest_size=16).digest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract raw Dante variants missing from SNP/indel splits")
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--snp", type=Path, required=True)
    parser.add_argument("--indel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--progress-every", type=int, default=500_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for source in (args.raw, args.snp, args.indel):
        if not source.is_file():
            raise FileNotFoundError(source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    fd, db_name = tempfile.mkstemp(prefix="dante_residual_", suffix=".sqlite3", dir=args.output.parent)
    os.close(fd)
    db_path = Path(db_name)
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    conn.execute("CREATE TABLE seen (k BLOB PRIMARY KEY) WITHOUT ROWID")
    counts = {"snp": 0, "indel": 0, "raw": 0, "residual": 0}
    started = time.time()

    part = args.output.with_name(args.output.name + ".part")
    try:
        for label, source in (("snp", args.snp), ("indel", args.indel)):
            conn.execute("BEGIN")
            for raw in iter_records(source):
                conn.execute("INSERT OR IGNORE INTO seen(k) VALUES (?)", (key_from_raw(raw),))
                counts[label] += 1
                if args.progress_every and counts[label] % args.progress_every == 0:
                    print(f"  {source.name}: {counts[label]:,} indexed", flush=True)
            conn.commit()

        with open_anonymous_gzip_text(part) as output:
            with open_text(args.raw) as raw_handle:
                header_done = False
                for raw in raw_handle:
                    if raw.startswith("#"):
                        if raw.startswith("#CHROM") and not header_done:
                            output.write(
                                "##dante_raw_residual=records absent from filtered.snp and filtered.indel by CHROM,POS,REF,ALT\n"
                            )
                            header_done = True
                        output.write(raw)
                        continue
                    if not raw.strip():
                        continue
                    counts["raw"] += 1
                    if conn.execute("SELECT 1 FROM seen WHERE k=?", (key_from_raw(raw),)).fetchone() is None:
                        output.write(raw if raw.endswith("\n") else raw + "\n")
                        counts["residual"] += 1
                    if args.progress_every and counts["raw"] % args.progress_every == 0:
                        print(f"  {args.raw.name}: {counts['raw']:,} checked", flush=True)
        os.replace(part, args.output)
    except Exception:
        part.unlink(missing_ok=True)
        raise
    finally:
        conn.close()
        db_path.unlink(missing_ok=True)

    lines = [
        "DANTE RAW RESIDUAL EXTRACTION REPORT",
        "=" * 40,
        "Match key: CHROM, POS, REF, ALT",
        f"SNP split rows indexed: {counts['snp']:,}",
        f"Indel split rows indexed: {counts['indel']:,}",
        f"Raw rows checked: {counts['raw']:,}",
        f"Raw residual rows written: {counts['residual']:,}",
        f"Elapsed seconds: {time.time() - started:.1f}",
        "",
        "The original VCF files were not modified.",
        "Residual records are typically mixed/complex variants excluded from pure SNP and indel splits.",
        "",
    ]
    args.report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Created {args.output.name}: {counts['residual']:,} residual records")
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
