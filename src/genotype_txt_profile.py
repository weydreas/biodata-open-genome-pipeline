#!/usr/bin/env python3
"""Detect and annotate 23andMe-style genotype-table patient inputs.

This is intentionally a third patient-input contract.  It does not reclassify
Dante WGS bundles or array VCFs.  The PowerShell runtime calls this helper only
when the existing input-profile resolver returned UNKNOWN.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
from collections import Counter
from pathlib import Path
from typing import Any, TextIO

SCHEMA = "dante-patient-input-profile-v1"
KIND = "GENOTYPE_TABLE_TXT"
UNKNOWN_KIND = "UNKNOWN"
CANONICAL_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
GENOTYPE_RE = re.compile(r"^(?:--|[ACGT]{1,2}|[ID]{1,2})$", re.IGNORECASE)
TXT_RE = re.compile(r"\.(?:txt|tsv)(?:\.gz)?$", re.IGNORECASE)


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def quick_identity(path: Path) -> dict[str, Any]:
    st = path.stat()
    return {
        "name": path.name,
        "length": int(st.st_size),
        "modified_utc_ticks": int(st.st_mtime_ns // 100),
    }


def normalize_chrom(value: str) -> str:
    text = value.strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.upper() == "M":
        return "MT"
    return text.upper() if text.upper() in {"X", "Y", "MT"} else text


def inspect_table(path: Path, max_records: int = 10_000) -> dict[str, Any]:
    comments: list[str] = []
    counts: Counter[str] = Counter()
    header_contract = False
    source = "UNKNOWN"
    assembly = "UNKNOWN"

    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("#"):
                lower = line.lower()
                if len(comments) < 80:
                    comments.append(line)
                if "23andme" in lower:
                    source = "23andMe"
                if (
                    "build 37" in lower
                    or "grch37" in lower
                    or "hg19" in lower
                    or "annotation release 104" in lower
                ):
                    assembly = "GRCh37"
                if all(token in lower for token in ("rsid", "chromosome", "position", "genotype")):
                    header_contract = True
                continue
            if not line.strip():
                continue

            fields = line.split("\t")
            counts["sampled"] += 1
            if len(fields) != 4:
                counts["malformed"] += 1
                if counts["sampled"] >= max_records:
                    break
                continue
            marker, chrom_raw, pos_raw, genotype_raw = fields
            chrom = normalize_chrom(chrom_raw)
            genotype = genotype_raw.strip().upper()
            if chrom in CANONICAL_CHROMS:
                counts["canonical_chrom"] += 1
            try:
                if int(pos_raw) > 0:
                    counts["valid_position"] += 1
            except ValueError:
                pass
            if GENOTYPE_RE.fullmatch(genotype):
                counts["recognized_genotype"] += 1
            if re.fullmatch(r"rs\d+", marker, re.IGNORECASE):
                counts["rsid"] += 1
            elif re.fullmatch(r"i\d+", marker, re.IGNORECASE):
                counts["internal_id"] += 1
            if genotype == "--":
                counts["no_call"] += 1
            elif set(genotype) <= set("ACGT"):
                if len(genotype) == 1:
                    counts["haploid_acgt"] += 1
                elif len(genotype) == 2:
                    counts["diploid_acgt"] += 1
            elif set(genotype) <= set("ID"):
                counts["symbolic_indel"] += 1

            if counts["sampled"] >= max_records:
                break

    sampled = max(1, counts["sampled"])
    if source == "23andMe" and assembly == "UNKNOWN":
        # 23andMe raw files must still explicitly expose the build contract;
        # source branding alone is not enough to infer GRCh37.
        assembly = "UNKNOWN"

    return {
        "source": source,
        "assembly": assembly,
        "header_contract": bool(header_contract),
        "sampled_records": int(counts["sampled"]),
        "canonical_chrom_fraction": counts["canonical_chrom"] / sampled,
        "valid_position_fraction": counts["valid_position"] / sampled,
        "recognized_genotype_fraction": counts["recognized_genotype"] / sampled,
        "id_fraction": (counts["rsid"] + counts["internal_id"]) / sampled,
        "counts": dict(counts),
        "comments": comments,
    }


def classify(path: Path) -> tuple[str, dict[str, Any]]:
    detail = inspect_table(path)
    supported = (
        detail["assembly"] == "GRCh37"
        and detail["header_contract"]
        and detail["sampled_records"] >= 10
        and detail["canonical_chrom_fraction"] >= 0.95
        and detail["valid_position_fraction"] >= 0.95
        and detail["recognized_genotype_fraction"] >= 0.95
        and detail["id_fraction"] >= 0.50
    )
    return (KIND if supported else UNKNOWN_KIND), detail


def detect(patient_root: Path) -> dict[str, Any]:
    patient_root = patient_root.resolve()
    input_dir = patient_root / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    candidates = sorted(
        (p for p in input_dir.iterdir() if p.is_file() and TXT_RE.search(p.name)),
        key=lambda p: p.name.casefold(),
    )
    matches: list[tuple[Path, dict[str, Any]]] = []
    inspected: list[dict[str, Any]] = []
    for path in candidates:
        kind, detail = classify(path)
        inspected.append({"file": quick_identity(path), "kind": kind, "inspection": detail})
        if kind == KIND:
            matches.append((path, detail))

    if len(matches) == 1:
        source_file, detail = matches[0]
        source = detail.get("source") or "GENOTYPE_TABLE_TXT"
        if source == "UNKNOWN":
            source = "23andMe-style genotype table"
        return {
            "schema": SCHEMA,
            "kind": KIND,
            "assembly": "GRCh37",
            "source": source,
            "scope": (
                "Interrogated genotype-table loci only. Calls are accepted only at rows present in the source; "
                "-- is a no-call; absence outside the table is UNKNOWN. Symbolic I/D indel codes do not encode "
                "sequence alleles and are excluded from sequence-level VEP annotation."
            ),
            "capabilities": {
                "direct_rsid_export": True,
                "small_variant_annotation": True,
                "structural_variant_annotation": False,
                "fastq_alignment": False,
                "reference_genome_completion": False,
                "telomere_from_fastq": False,
                "local_bam_browser": False,
                "knowledge_enrichment": True,
                "final_reports": True,
                "whole_genome_negative_evidence": False,
                "sequence_resolved_indels": False,
            },
            "source_file": quick_identity(source_file),
            "inspection": detail,
            "txt_candidates": inspected,
        }

    return {
        "schema": SCHEMA,
        "kind": UNKNOWN_KIND,
        "assembly": "UNKNOWN",
        "source": "UNKNOWN",
        "scope": "No unique supported 23andMe-style GRCh37 genotype table was detected.",
        "capabilities": {
            "direct_rsid_export": False,
            "small_variant_annotation": False,
            "structural_variant_annotation": False,
            "fastq_alignment": False,
            "reference_genome_completion": False,
            "telomere_from_fastq": False,
            "local_bam_browser": False,
            "knowledge_enrichment": False,
            "final_reports": False,
            "whole_genome_negative_evidence": False,
            "sequence_resolved_indels": False,
        },
        "source_file": None,
        "txt_candidates": inspected,
    }


def save_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(part, path)


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def annotate_stage24(patient_root: Path) -> dict[str, Any]:
    patient_root = patient_root.resolve()
    profile_path = patient_root / "input_profile.json"
    profile = detect(patient_root)
    if profile.get("kind") != KIND:
        raise RuntimeError("active patient is not a supported genotype-table TXT profile")
    save_atomic(profile_path, profile)

    output_dir = patient_root / "Dante_results" / "knowledge_enrichment"
    state_path = output_dir / "Dante_knowledge_enrichment.state.json"
    report_path = output_dir / "Dante_knowledge_enrichment_report.txt"
    if not state_path.is_file() or not report_path.is_file():
        raise FileNotFoundError("Stage 24 state/report is missing")

    state = json.loads(state_path.read_text(encoding="utf-8"))
    scope_block = [
        "INPUT SCOPE",
        f"Input profile: {KIND}",
        f"Source: {profile.get('source', 'UNKNOWN')}",
        "Assembly: GRCh37",
        f"Scope: {profile.get('scope', '')}",
        "Whole-genome negative evidence: NO",
        "Structural variants: NOT TESTED by this source; CNV/SV Stage-24 inputs are explicit N/A sentinels only.",
        "Technical confidence: genotype-table calls have no read-level DP/GQ/QUAL/FILTER evidence contract.",
        "Symbolic I/D indel calls: NOT sequence-resolved; retained in direct genotype export/QC but excluded from VEP and Stage 24.",
    ]
    marker = "INPUT SCOPE\n"
    report_text = report_path.read_text(encoding="utf-8")
    if marker in report_text:
        before, after = report_text.split(marker, 1)
        if "\n\n" in after:
            _, tail = after.split("\n\n", 1)
            report_text = before + "\n".join(scope_block) + "\n\n" + tail
        else:
            report_text = before + "\n".join(scope_block) + "\n"
    else:
        report_text = "\n".join(scope_block) + "\n\n" + report_text
    report_path.write_text(report_text, encoding="utf-8")

    state["input_profile"] = {
        "schema": SCHEMA,
        "kind": KIND,
        "source": profile.get("source"),
        "assembly": "GRCh37",
        "scope": profile.get("scope"),
        "whole_genome_negative_evidence": False,
        "sequence_resolved_indels": False,
        "profile_sha256": sha256_file(profile_path),
    }
    for identity in (state.get("outputs") or {}).values():
        if identity.get("name") == report_path.name:
            identity["length"] = report_path.stat().st_size
            identity["sha256"] = sha256_file(report_path)
            break
    save_atomic(state_path, state)
    return state


def _synthetic_text() -> str:
    rows = [
        "# This data file generated by 23andMe at: synthetic",
        "# We are using reference human assembly build 37 (also known as Annotation Release 104).",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    for i in range(1, 21):
        gt = "AG" if i % 3 == 0 else "AA"
        rows.append(f"rs{i}\t1\t{i}\t{gt}")
    rows.extend(["rs100\tX\t1\tA", "rs101\tY\t1\tG", "rs102\t1\t30\tII", "rs103\t1\t31\t--"])
    return "\n".join(rows) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("detect", "annotate-stage24", "self-test"))
    ap.add_argument("--patient-root", type=Path)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--save", action="store_true")
    args = ap.parse_args()

    if args.mode == "self-test":
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inp = root / "input"
            inp.mkdir()
            (inp / "sample_normal.txt").write_text(_synthetic_text(), encoding="utf-8")
            profile = detect(root)
            assert profile["kind"] == KIND
            assert profile["assembly"] == "GRCh37"
            assert profile["source"] == "23andMe"
            assert profile["capabilities"]["whole_genome_negative_evidence"] is False
            assert profile["capabilities"]["sequence_resolved_indels"] is False
        print("SELF-TEST PASS | independent GRCh37 genotype-table TXT profile detection")
        return 0

    if args.patient_root is None:
        ap.error("--patient-root is required")
    if args.mode == "detect":
        profile = detect(args.patient_root)
        if args.save and profile.get("kind") == KIND:
            save_atomic(args.patient_root.resolve() / "input_profile.json", profile)
        if args.json:
            print(json.dumps(profile, ensure_ascii=False, sort_keys=True, separators=(",", ":")))
        else:
            print(profile["kind"])
        return 0

    state = annotate_stage24(args.patient_root)
    if args.json:
        print(json.dumps(state.get("input_profile", {}), ensure_ascii=False, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
