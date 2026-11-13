#!/usr/bin/env python3
"""Patient workspace selector and one-time v8.5 legacy-layout migrator.

Shared installation data stay in the pipeline root:
  annotation_data/
  reference_data/
  scripts / docs

Patient-specific data live under:
  patients/<Patient Name>/
      input/
          original source VCF/FASTQ files
          reference_panel_sources/
      Dante_work/
      Dante_qc/
      Dante_results/
      .env
      patient.json

The v8.5.0 bootstrap migrates the pre-patient-layout project into
patients/Demo 1 without recomputing genomic outputs.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import shutil
import sys
from typing import Iterable

SCHEMA = "dante-patient-workspace-v1"
ACTIVE_SCHEMA = "dante-active-patient-v1"
MIGRATION_SCHEMA = "dante-patient-layout-migration-v1"
PIPELINE_RELEASE = "8.5.2"
DEFAULT_PATIENT = "Demo 1"

PATIENT_DIRS = ("Dante_work", "Dante_qc", "Dante_results")
SAMPLE_ENV_KEYS = {
    "DANTE_SAMPLE_DATE",
    "DANTE_SAMPLE_YEAR",
    "DANTE_SAMPLE_MATERIAL",
    "DANTE_SAMPLE_SEX",
}
INVALID_WINDOWS_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
RESERVED_WINDOWS_NAMES = {
    "CON", "PRN", "AUX", "NUL",
    *(f"COM{i}" for i in range(1, 10)),
    *(f"LPT{i}" for i in range(1, 10)),
}

SOURCE_PATTERNS = (
    re.compile(r"(?i)\.raw\.vcf(?:\.gz)?$"),
    re.compile(r"(?i)\.filtered\.snp\.vcf(?:\.gz)?$"),
    re.compile(r"(?i)\.filtered\.indel\.vcf(?:\.gz)?$"),
    re.compile(r"(?i)\.cnv\.vcf(?:\.gz)?$"),
    re.compile(r"(?i)\.sv\.vcf(?:\.gz)?$"),
    re.compile(r"(?i)\.(?:fastq|fq)(?:\.gz)?$"),
)
SOURCE_COMPANION_PATTERNS = (
    re.compile(r"(?i)\.(?:vcf|fastq|fq)(?:\.gz)?\.(?:tbi|csi|md5|sha256)$"),
    re.compile(r"(?i)\.(?:fastq|fq)(?:\.gz)?\.md5$"),
)

def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()

def configure_console() -> None:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

def validate_patient_name(name: str) -> str:
    value = name.strip()
    if not value:
        raise ValueError("Patient name cannot be empty.")
    if value in {".", ".."}:
        raise ValueError("Invalid patient name.")
    if INVALID_WINDOWS_CHARS.search(value):
        raise ValueError('Patient name contains a Windows-invalid character: <>:"/\\|?*')
    if value.endswith((" ", ".")):
        raise ValueError("Patient name cannot end with a space or period on Windows.")
    stem = value.split(".", 1)[0].upper()
    if stem in RESERVED_WINDOWS_NAMES:
        raise ValueError(f"Reserved Windows name: {value}")
    return value

def patient_root(root: Path, name: str) -> Path:
    return root / "patients" / validate_patient_name(name)

def patient_input(root: Path, name: str) -> Path:
    return patient_root(root, name) / "input"

def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(part, path)

def ensure_patient(root: Path, name: str) -> Path:
    name = validate_patient_name(name)
    p = patient_root(root, name)
    inp = p / "input"
    panels = inp / "reference_panel_sources"
    for directory in (p, inp, panels):
        directory.mkdir(parents=True, exist_ok=True)
    meta = p / "patient.json"
    if not meta.is_file():
        write_json_atomic(meta, {
            "schema": SCHEMA,
            "name": name,
            "created_utc": utcnow(),
            "pipeline_release_created": PIPELINE_RELEASE,
        })
    template = root / "reference_panel_sources" / "README.md"
    patient_readme = panels / "README.md"
    if template.is_file() and not patient_readme.exists():
        shutil.copy2(template, patient_readme)
    return p

def is_source_file(path: Path) -> bool:
    if not path.is_file():
        return False
    return any(p.search(path.name) for p in SOURCE_PATTERNS + SOURCE_COMPANION_PATTERNS)

def list_legacy_source_files(root: Path) -> list[Path]:
    return sorted((p for p in root.iterdir() if is_source_file(p)), key=lambda p: p.name.lower())

def move_file_no_overwrite(src: Path, dst: Path, moved: list[dict]) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        raise RuntimeError(
            f"Migration conflict: destination already exists: {dst}. "
            "Nothing was overwritten; resolve the duplicate manually."
        )
    stat = src.stat()
    os.replace(src, dst)
    moved.append({
        "kind": "file",
        "from": str(src),
        "to": str(dst),
        "length": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    })

def move_tree_no_overwrite(src: Path, dst: Path, moved: list[dict]) -> None:
    if not src.exists():
        return
    if not dst.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.replace(src, dst)
        moved.append({"kind": "directory", "from": str(src), "to": str(dst)})
        return
    if not src.is_dir() or not dst.is_dir():
        raise RuntimeError(f"Migration conflict: cannot merge {src} into {dst}")
    for child in sorted(src.iterdir(), key=lambda p: p.name.lower()):
        target = dst / child.name
        if child.is_dir():
            move_tree_no_overwrite(child, target, moved)
        else:
            move_file_no_overwrite(child, target, moved)
    try:
        src.rmdir()
    except OSError:
        pass

def migrate_reference_panel_sources(root: Path, target: Path, moved: list[dict]) -> None:
    legacy = root / "reference_panel_sources"
    target.mkdir(parents=True, exist_ok=True)
    if not legacy.is_dir():
        return
    for child in sorted(legacy.iterdir(), key=lambda p: p.name.lower()):
        if child.name.lower() == "readme.md":
            continue
        destination = target / child.name
        if child.is_dir():
            move_tree_no_overwrite(child, destination, moved)
        else:
            move_file_no_overwrite(child, destination, moved)

def migrate_patient_env(root: Path, patient: Path, moved: list[dict]) -> None:
    global_env = root / ".env"
    if not global_env.is_file():
        return
    lines = global_env.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    keep: list[str] = []
    patient_lines: list[str] = []
    for line in lines:
        stripped = line.strip()
        match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=", stripped)
        if match and match.group(1) in SAMPLE_ENV_KEYS:
            patient_lines.append(line)
        else:
            keep.append(line)
    if not patient_lines:
        return

    patient_env = patient / ".env"
    existing = patient_env.read_text(encoding="utf-8-sig", errors="replace").splitlines() if patient_env.exists() else []
    existing_keys = set()
    for line in existing:
        m = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
        if m:
            existing_keys.add(m.group(1))
    additions = [
        line for line in patient_lines
        if (re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", line).group(1) not in existing_keys)
    ]
    if additions:
        out = list(existing)
        if out and out[-1].strip():
            out.append("")
        out.append("# Migrated from root .env by Dante WGS v8.5.0 patient-layout migration")
        out.extend(additions)
        patient_env.write_text("\n".join(out).rstrip() + "\n", encoding="utf-8")
        moved.append({"kind": "env-patient-keys", "to": str(patient_env), "keys": sorted(SAMPLE_ENV_KEYS)})

    backup = root / ".env.pre-v8.5.0.bak"
    if not backup.exists():
        shutil.copy2(global_env, backup)
    global_env.write_text("\n".join(keep).rstrip() + ("\n" if keep else ""), encoding="utf-8")
    moved.append({
        "kind": "env-global-cleanup",
        "from": str(global_env),
        "backup": str(backup),
        "removed_keys": sorted(
            re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", line).group(1)
            for line in patient_lines
        ),
    })

def legacy_layout_present(root: Path) -> bool:
    if any((root / name).exists() for name in PATIENT_DIRS):
        return True
    if list_legacy_source_files(root):
        return True
    panel_dir = root / "reference_panel_sources"
    if panel_dir.is_dir():
        if any(p.name.lower() != "readme.md" for p in panel_dir.iterdir()):
            return True
    return False

def bootstrap_legacy_patient(root: Path) -> dict:
    """One-time conversion from the pre-v8.5 root layout."""
    patients = root / "patients"
    patients.mkdir(parents=True, exist_ok=True)
    # An upgrade must reuse the previous migration workspace, whatever its name.
    # Never rename local patient folders or rewrite patient metadata.
    markers = sorted(
        p / "migration_v8.5.0.json" for p in patients.iterdir()
        if p.is_dir() and (p / "migration_v8.5.0.json").is_file()
    )
    if len(markers) > 1:
        raise RuntimeError("Multiple legacy migration markers; manual review required.")
    if markers:
        return json.loads(markers[0].read_text(encoding="utf-8"))

    # An already initialized patient layout is not a new patient.
    # In particular, never create an empty Demo 1 beside an existing workspace.
    if not legacy_layout_present(root) and any(
        p.is_dir() and (p / "patient.json").is_file() for p in patients.iterdir()
    ):
        return {"schema": MIGRATION_SCHEMA, "status": "already_initialized", "moved": []}

    patient = ensure_patient(root, DEFAULT_PATIENT)
    marker = patient / "migration_v8.5.0.json"
    if marker.is_file():
        return json.loads(marker.read_text(encoding="utf-8"))

    moved: list[dict] = []
    had_legacy = legacy_layout_present(root)

    # Generated patient-specific directories.
    for name in PATIENT_DIRS:
        move_tree_no_overwrite(root / name, patient / name, moved)

    # Original source genome files from which the project started.
    inp = patient / "input"
    for src in list_legacy_source_files(root):
        move_file_no_overwrite(src, inp / src.name, moved)

    # Original array/reference-panel inputs are patient-specific too.
    migrate_reference_panel_sources(root, inp / "reference_panel_sources", moved)

    # Sample-specific metadata must not leak to patient #2.
    migrate_patient_env(root, patient, moved)

    payload = {
        "schema": MIGRATION_SCHEMA,
        "status": "complete",
        "pipeline_release": PIPELINE_RELEASE,
        "patient": DEFAULT_PATIENT,
        "performed_utc": utcnow(),
        "legacy_layout_detected": had_legacy,
        "moved": moved,
        "notes": [
            "annotation_data and reference_data remain shared at pipeline root",
            "root reference_panel_sources/README.md remains as a template only",
            "existing genomic outputs were moved, not recomputed",
        ],
    }
    write_json_atomic(marker, payload)
    return payload

def list_patients(root: Path) -> list[str]:
    base = root / "patients"
    if not base.is_dir():
        return []
    names = []
    for p in base.iterdir():
        if p.is_dir() and (p / "patient.json").is_file():
            try:
                meta = json.loads((p / "patient.json").read_text(encoding="utf-8"))
                names.append(str(meta.get("name") or p.name))
            except Exception:
                names.append(p.name)
    return sorted(set(names), key=str.casefold)

def patient_summary(root: Path, name: str) -> str:
    p = patient_root(root, name)
    inp = p / "input"
    source_count = sum(1 for x in inp.iterdir() if is_source_file(x)) if inp.is_dir() else 0
    results = p / "Dante_results"
    qc = p / "Dante_qc"
    flags = []
    if source_count:
        flags.append(f"{source_count} source files")
    if results.is_dir():
        flags.append("results")
    if qc.is_dir():
        flags.append("QC")
    return ", ".join(flags) if flags else "empty workspace"

def activate(root: Path, name: str) -> None:
    p = ensure_patient(root, name)
    write_json_atomic(root / ".active_patient.json", {
        "schema": ACTIVE_SCHEMA,
        "name": name,
        "relative_path": str(p.relative_to(root)).replace("\\", "/"),
        "selected_utc": utcnow(),
    })

def create_interactive(root: Path) -> str:
    while True:
        value = input("New patient name: ").strip()
        try:
            name = validate_patient_name(value)
            if patient_root(root, name).exists():
                print(f"Patient already exists: {name}")
                continue
            ensure_patient(root, name)
            return name
        except ValueError as exc:
            print(f"Invalid name: {exc}")

def select_interactive(root: Path, forced_patient: str | None = None) -> int:
    migration = bootstrap_legacy_patient(root)
    if migration.get("moved"):
        print("")
        print("v8.5.0 patient-layout migration completed:")
        print(f"  Patient: {DEFAULT_PATIENT}")
        print(f"  Moved entries: {len(migration['moved'])}")
        print(f"  Folder: {patient_root(root, DEFAULT_PATIENT)}")

    if forced_patient:
        name = validate_patient_name(forced_patient)
        ensure_patient(root, name)
        activate(root, name)
        print(f"Selected patient: {name}")
        return 0

    while True:
        patients = list_patients(root)
        if not patients:
            ensure_patient(root, DEFAULT_PATIENT)
            patients = list_patients(root)

        print("")
        print("=" * 78)
        print("SELECT PATIENT")
        print("=" * 78)
        for i, name in enumerate(patients, 1):
            print(f"{i:2}. {name}  [{patient_summary(root, name)}]")
        print(" N. Create new patient")
        print(" 0. Exit")
        choice = input("Select patient: ").strip()

        if choice == "0":
            return 10
        if choice.lower() in {"n", "new", "+"}:
            name = create_interactive(root)
            activate(root, name)
            print(f"Selected patient: {name}")
            return 0
        if choice.isdigit():
            index = int(choice)
            if 1 <= index <= len(patients):
                name = patients[index - 1]
                activate(root, name)
                print(f"Selected patient: {name}")
                return 0
        print("Unknown selection.")

def self_test() -> int:
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "reference_panel_sources").mkdir()
        (root / "reference_panel_sources" / "README.md").write_text("template\n", encoding="utf-8")
        (root / "reference_panel_sources" / "23andme.tsv").write_text("rs1\t1\t1\tAA\n", encoding="utf-8")
        (root / "sample.raw.vcf.gz").write_bytes(b"raw")
        (root / "sample.filtered.snp.vcf.gz").write_bytes(b"snp")
        (root / "sample_R1_001.fastq.gz").write_bytes(b"fq")
        (root / "Dante_results").mkdir()
        (root / "Dante_results" / "out.txt").write_text("done", encoding="utf-8")
        (root / "Dante_qc").mkdir()
        (root / "Dante_qc" / "qc.txt").write_text("qc", encoding="utf-8")
        (root / "Dante_work").mkdir()
        (root / "Dante_work" / "state.json").write_text("{}", encoding="utf-8")
        (root / ".env").write_text(
            "OPENAI_API_KEY=x\nDANTE_SAMPLE_SEX=M\nDANTE_SAMPLE_YEAR=2019\n",
            encoding="utf-8",
        )

        result = bootstrap_legacy_patient(root)
        p = root / "patients" / DEFAULT_PATIENT
        assert result["status"] == "complete"
        assert (p / "input" / "sample.raw.vcf.gz").is_file()
        assert (p / "input" / "sample.filtered.snp.vcf.gz").is_file()
        assert (p / "input" / "sample_R1_001.fastq.gz").is_file()
        assert (p / "input" / "reference_panel_sources" / "23andme.tsv").is_file()
        assert (p / "Dante_results" / "out.txt").is_file()
        assert (p / "Dante_qc" / "qc.txt").is_file()
        assert (p / "Dante_work" / "state.json").is_file()
        assert (root / "reference_panel_sources" / "README.md").is_file()
        assert "DANTE_SAMPLE_SEX" not in (root / ".env").read_text(encoding="utf-8")
        assert "OPENAI_API_KEY=x" in (root / ".env").read_text(encoding="utf-8")
        assert "DANTE_SAMPLE_SEX=M" in (p / ".env").read_text(encoding="utf-8")
        activate(root, DEFAULT_PATIENT)
        active = json.loads((root / ".active_patient.json").read_text(encoding="utf-8"))
        assert active["name"] == DEFAULT_PATIENT
        assert bootstrap_legacy_patient(root)["status"] == "complete"
    print("SELF-TEST PASS | patient migration + source isolation + env isolation + active selection")
    return 0

def main() -> int:
    configure_console()
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("select", "migrate", "list", "self-test"))
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--patient")
    args = parser.parse_args()
    root = args.root.resolve()

    if args.mode == "self-test":
        return self_test()
    if args.mode == "migrate":
        payload = bootstrap_legacy_patient(root)
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0
    if args.mode == "list":
        bootstrap_legacy_patient(root)
        for name in list_patients(root):
            print(name)
        return 0
    return select_interactive(root, args.patient)

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nPatient selection cancelled.", file=sys.stderr)
        raise SystemExit(10)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
