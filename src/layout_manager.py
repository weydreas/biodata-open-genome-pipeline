#!/usr/bin/env python3
"""One-time v8.6.0 migration from a flat code layout to src/ + docs/.

This migrator is intentionally conservative:
- patient data, annotation_data and reference_data are never touched;
- root *.py/*.ps1 files are moved under src/_legacy_root_before_v8.6.0/;
- known root documentation files are moved under docs/_legacy_root_before_v8.6.0/;
- the canonical v8.6.0 code already shipped under src/ is never overwritten;
- migration is idempotent and leaves an audit JSON under src/.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import sys
from typing import Any

PIPELINE_RELEASE = "8.6.0"
SCHEMA = "dante-code-layout-migration-v1"
MARKER_NAME = ".layout_migration_v8.6.0.json"
LEGACY_CODE_DIR = "_legacy_root_before_v8.6.0"
LEGACY_DOCS_DIR = "_legacy_root_before_v8.6.0"
KNOWN_DOCS = {
    "BROWSER_THIRD_PARTY_NOTICE.md",
    "CHANGELOG.md",
    "DANTE_WGS_V7_README_RU.md",
    "DANTE_WGS_V8_README_RU.md",
    "DANTE_WGS_V9_README_RU.md",
    "STAGE24_KNOWLEDGE_ENRICHMENT_CONTRACT_RU.md",
}
CODE_SUFFIXES = {".py", ".ps1"}
CODE_EXACT_NAMES = {"requirements-browser.txt"}

def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()

def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(part, path)

def unique_destination(directory: Path, name: str) -> Path:
    candidate = directory / name
    if not candidate.exists():
        return candidate
    stem = Path(name).stem
    suffix = Path(name).suffix
    i = 2
    while True:
        candidate = directory / f"{stem}__legacy_{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1

def move_preserving(src: Path, dst_dir: Path, records: list[dict[str, Any]]) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = unique_destination(dst_dir, src.name)
    stat = src.stat()
    os.replace(src, dst)
    records.append({
        "from": str(src),
        "to": str(dst),
        "length": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    })

def candidates(root: Path) -> tuple[list[Path], list[Path]]:
    code: list[Path] = []
    docs: list[Path] = []
    for p in root.iterdir():
        if not p.is_file():
            continue
        if p.suffix.lower() in CODE_SUFFIXES or p.name in CODE_EXACT_NAMES:
            code.append(p)
        elif p.name in KNOWN_DOCS:
            docs.append(p)
    return sorted(code), sorted(docs)

def migrate(root: Path) -> dict[str, Any]:
    root = root.resolve()
    src = root / "src"
    docs = root / "docs"
    marker = src / MARKER_NAME

    if not src.is_dir():
        raise RuntimeError(
            f"Canonical src directory is missing: {src}. "
            "Extract the full v9.0.0 archive over the pipeline root before launching."
        )
    canonical = src / "dante_wgs_pipeline.ps1"
    if not canonical.is_file():
        raise RuntimeError(f"Canonical pipeline script is missing: {canonical}")

    # Even if a marker already exists, re-check root for leftovers from an older overlay.
    # This makes the migration self-healing without touching data directories.
    code_files, doc_files = candidates(root)
    moved_code: list[dict[str, Any]] = []
    moved_docs: list[dict[str, Any]] = []

    code_backup = src / LEGACY_CODE_DIR
    docs_backup = docs / LEGACY_DOCS_DIR
    for p in code_files:
        move_preserving(p, code_backup, moved_code)
    for p in doc_files:
        move_preserving(p, docs_backup, moved_docs)

    previous = {}
    if marker.is_file():
        try:
            previous = json.loads(marker.read_text(encoding="utf-8"))
        except Exception:
            previous = {}

    history = list(previous.get("runs", []))
    history.append({
        "utc": utcnow(),
        "moved_code_files": moved_code,
        "moved_doc_files": moved_docs,
    })
    payload = {
        "schema": SCHEMA,
        "status": "complete",
        "pipeline_release": PIPELINE_RELEASE,
        "root": str(root),
        "canonical_src": str(src),
        "canonical_docs": str(docs),
        "runs": history[-20:],
        "last_run_moved_code": len(moved_code),
        "last_run_moved_docs": len(moved_docs),
        "protected_data_directories": [
            "patients",
            "annotation_data",
            "reference_data",
            "reference_panel_sources",
            "tests",
        ],
    }
    write_json_atomic(marker, payload)
    return payload

def self_test() -> int:
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        src = root / "src"
        docs = root / "docs"
        src.mkdir()
        docs.mkdir()
        (src / "dante_wgs_pipeline.ps1").write_text("# canonical\n", encoding="utf-8")
        (src / "patient_manager.py").write_text("# canonical\n", encoding="utf-8")

        # Legacy flat code + docs + data that MUST remain untouched.
        (root / "old_a.py").write_text("print('old')\n", encoding="utf-8")
        (root / "old_b.ps1").write_text("Write-Host old\n", encoding="utf-8")
        (root / "requirements-browser.txt").write_text("Flask\n", encoding="utf-8")
        (root / "CHANGELOG.md").write_text("legacy\n", encoding="utf-8")
        (root / "user_notes.txt").write_text("do not move\n", encoding="utf-8")
        (root / "patients").mkdir()
        (root / "patients" / "keep.dat").write_bytes(b"x")
        (root / "reference_data").mkdir()
        (root / "reference_data" / "keep.dat").write_bytes(b"x")

        result = migrate(root)
        assert result["status"] == "complete"
        assert not (root / "old_a.py").exists()
        assert not (root / "old_b.ps1").exists()
        assert not (root / "requirements-browser.txt").exists()
        assert not (root / "CHANGELOG.md").exists()
        assert (src / LEGACY_CODE_DIR / "old_a.py").is_file()
        assert (src / LEGACY_CODE_DIR / "old_b.ps1").is_file()
        assert (docs / LEGACY_DOCS_DIR / "CHANGELOG.md").is_file()
        assert (root / "user_notes.txt").is_file()
        assert (root / "patients" / "keep.dat").is_file()
        assert (root / "reference_data" / "keep.dat").is_file()
        assert (src / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8") == "# canonical\n"

        # Idempotent second run.
        second = migrate(root)
        assert second["last_run_moved_code"] == 0
        assert second["last_run_moved_docs"] == 0

    print("SELF-TEST PASS | flat-code migration + canonical protection + data-directory protection")
    return 0

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("migrate", "self-test"))
    ap.add_argument("--root", type=Path, default=Path.cwd())
    args = ap.parse_args()
    if args.mode == "self-test":
        return self_test()
    payload = migrate(args.root)
    if payload["last_run_moved_code"] or payload["last_run_moved_docs"]:
        print("")
        print("v8.6.0 code-layout migration completed:")
        print(f"  old root code moved: {payload['last_run_moved_code']}")
        print(f"  old root docs moved: {payload['last_run_moved_docs']}")
        print(f"  canonical code: {payload['canonical_src']}")
        print(f"  canonical docs: {payload['canonical_docs']}")
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Code-layout migration cancelled.", file=sys.stderr)
        raise SystemExit(10)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
