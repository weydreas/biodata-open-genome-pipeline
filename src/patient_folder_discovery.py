#!/usr/bin/env python3
"""Non-destructive patient-folder discovery and input-profile refresh.

Immediate valid folders under patients/ are adopted without moving or overwriting
existing files. After adoption, every valid patient workspace gets a conservative
input_profile.json so the interactive pipeline can expose only analytically valid
stages for the data actually present (WGS bundle, genotyping-array VCF, or unknown).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import patient_input_profile as input_profile
import patient_manager as pm

DISCOVERY_SCHEMA = "dante-manual-patient-discovery-v2"


def discover_manual_patient_folders(root: Path) -> tuple[list[str], list[str]]:
    base = root / "patients"
    base.mkdir(parents=True, exist_ok=True)
    adopted: list[str] = []
    ignored: list[str] = []

    for folder in sorted(base.iterdir(), key=lambda p: p.name.casefold()):
        if not folder.is_dir() or folder.is_symlink():
            continue
        if folder.name.startswith("."):
            continue
        try:
            name = pm.validate_patient_name(folder.name)
        except ValueError:
            ignored.append(folder.name)
            continue

        if not (folder / "patient.json").is_file():
            # ensure_patient is additive/idempotent: it creates patient.json and
            # standard subdirectories but never moves/replaces existing inputs.
            pm.ensure_patient(root, name)
            adopted.append(name)

    return adopted, ignored


def refresh_input_profiles(root: Path) -> list[tuple[str, str]]:
    refreshed: list[tuple[str, str]] = []
    base = root / "patients"
    if not base.is_dir():
        return refreshed
    for folder in sorted(base.iterdir(), key=lambda p: p.name.casefold()):
        if (
            not folder.is_dir()
            or folder.is_symlink()
            or folder.name.startswith(".")
            or not (folder / "patient.json").is_file()
        ):
            continue
        profile = input_profile.detect(folder)
        input_profile.save_atomic(folder / "input_profile.json", profile)
        refreshed.append((folder.name, str(profile.get("kind") or "UNKNOWN")))
    return refreshed


def main() -> int:
    pm.configure_console()
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    root = args.root.resolve()

    adopted, ignored = discover_manual_patient_folders(root)
    profiles = refresh_input_profiles(root)
    for name in adopted:
        print(f"[PATIENTS] adopted manual folder: {name}")
    for name, kind in profiles:
        print(f"[PATIENTS] input profile: {name} -> {kind}")
    for name in ignored:
        print(f"[PATIENTS] ignored invalid folder name: {name}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
