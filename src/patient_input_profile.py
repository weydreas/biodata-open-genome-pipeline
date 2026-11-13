#!/usr/bin/env python3
"""Detect the analytical input scope of one Dante patient workspace.

The detector is deliberately conservative.  It does not infer "whole genome"
from a generic VCF.  A single dense genotype VCF with many homozygous-reference
loci is treated as a genotyping-array VCF; the absence of a locus from such a
file is not evidence about the rest of the genome.
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
ARRAY_KIND = "GENOTYPING_ARRAY_VCF"
WGS_KIND = "DANTE_WGS_BUNDLE"
UNKNOWN_KIND = "UNKNOWN"
RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)

DANTE_PATTERNS = {
    "raw": re.compile(r"\.raw\.vcf(?:\.gz)?$", re.IGNORECASE),
    "snp": re.compile(r"\.filtered\.snp\.vcf(?:\.gz)?$", re.IGNORECASE),
    "indel": re.compile(r"\.filtered\.indel\.vcf(?:\.gz)?$", re.IGNORECASE),
    "cnv": re.compile(r"\.cnv\.vcf(?:\.gz)?$", re.IGNORECASE),
    "sv": re.compile(r"\.sv\.vcf(?:\.gz)?$", re.IGNORECASE),
}
VCF_RE = re.compile(r"\.vcf(?:\.gz)?$", re.IGNORECASE)
FASTQ_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def normalize_reference(value: str) -> str:
    text = (value or "").strip().lower()
    if any(token in text for token in ("grch37", "hg19", "b37")):
        return "GRCh37"
    if any(token in text for token in ("grch38", "hg38")):
        return "GRCh38"
    return value.strip() or "UNKNOWN"


def quick_identity(path: Path) -> dict[str, Any]:
    st = path.stat()
    return {
        "name": path.name,
        "length": int(st.st_size),
        "modified_utc_ticks": int(st.st_mtime_ns // 100),
    }


def inspect_vcf(path: Path, max_records: int = 10_000) -> dict[str, Any]:
    meta: dict[str, str] = {}
    counts: Counter[str] = Counter()
    sample_names: list[str] = []
    format_ids: set[str] = set()
    filter_ids: set[str] = set()

    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("##source="):
                meta["source"] = raw.split("=", 1)[1].strip()
                continue
            if raw.startswith("##reference="):
                meta["reference"] = raw.split("=", 1)[1].strip()
                continue
            if raw.startswith("##FORMAT=<ID="):
                m = re.match(r"##FORMAT=<ID=([^,>]+)", raw)
                if m:
                    format_ids.add(m.group(1))
                continue
            if raw.startswith("##FILTER=<ID="):
                m = re.match(r"##FILTER=<ID=([^,>]+)", raw)
                if m:
                    filter_ids.add(m.group(1))
                continue
            if raw.startswith("#CHROM"):
                cols = raw.rstrip("\r\n").split("\t")
                sample_names = cols[9:] if len(cols) > 9 else []
                continue
            if raw.startswith("#") or not raw.strip():
                continue

            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                counts["malformed"] += 1
                continue
            counts["records_sampled"] += 1
            _, _, vid, ref, alt, qual, filt, _ = fields[:8]
            if vid and vid != ".":
                if RSID_RE.fullmatch(vid):
                    counts["rsid_records"] += 1
                else:
                    counts["named_non_rsid_records"] += 1
            else:
                counts["missing_id_records"] += 1
            if qual == ".":
                counts["qual_missing"] += 1
            if filt == ".":
                counts["filter_missing"] += 1
            if filt == "PASS":
                counts["filter_pass"] += 1

            alts = alt.split(",")
            if len(ref) == 1 and alts and all(len(a) == 1 for a in alts):
                counts["snv_records"] += 1
            else:
                counts["non_snv_records"] += 1

            if len(fields) >= 10 and sample_names:
                fmt = fields[8].split(":")
                vals = fields[9].split(":")
                sample = {k: vals[i] if i < len(vals) else "." for i, k in enumerate(fmt)}
                gt = sample.get("GT", ".")
                if gt in {".", "./.", ".|."}:
                    counts["missing_gt"] += 1
                else:
                    tokens = re.split(r"[/|]", gt)
                    if tokens and all(t.isdigit() for t in tokens):
                        alleles = [int(t) for t in tokens]
                        if all(x == 0 for x in alleles):
                            counts["hom_ref_gt"] += 1
                        elif any(x > 0 for x in alleles):
                            counts["non_ref_gt"] += 1
                    else:
                        counts["other_gt"] += 1

            if counts["records_sampled"] >= max_records:
                break

    sampled = max(1, counts["records_sampled"])
    return {
        "source": meta.get("source", "UNKNOWN"),
        "reference_raw": meta.get("reference", "UNKNOWN"),
        "assembly": normalize_reference(meta.get("reference", "")),
        "samples": sample_names,
        "sample_count": len(sample_names),
        "format_ids": sorted(format_ids),
        "filter_ids": sorted(filter_ids),
        "sampled_records": int(counts["records_sampled"]),
        "rsid_fraction": counts["rsid_records"] / sampled,
        "hom_ref_fraction": counts["hom_ref_gt"] / sampled,
        "non_ref_fraction": counts["non_ref_gt"] / sampled,
        "missing_gt_fraction": counts["missing_gt"] / sampled,
        "qual_missing_fraction": counts["qual_missing"] / sampled,
        "filter_missing_fraction": counts["filter_missing"] / sampled,
        "counts": dict(counts),
    }


def classify_single_vcf(path: Path) -> tuple[str, dict[str, Any]]:
    detail = inspect_vcf(path)
    is_array_like = (
        detail["sample_count"] == 1
        and detail["assembly"] == "GRCh37"
        and detail["rsid_fraction"] >= 0.50
        and detail["hom_ref_fraction"] >= 0.10
    )
    if is_array_like:
        return ARRAY_KIND, detail
    return UNKNOWN_KIND, detail


def detect(patient_root: Path) -> dict[str, Any]:
    patient_root = patient_root.resolve()
    input_dir = patient_root / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    files = sorted((p for p in input_dir.iterdir() if p.is_file()), key=lambda p: p.name.casefold())
    vcfs = [p for p in files if VCF_RE.search(p.name)]
    fastqs = [p for p in files if FASTQ_RE.search(p.name)]

    role_files: dict[str, list[Path]] = {
        role: [p for p in files if pattern.search(p.name)]
        for role, pattern in DANTE_PATTERNS.items()
    }
    wgs_core = all(len(role_files[x]) == 1 for x in ("raw", "snp", "indel"))
    wgs_structural = all(len(role_files[x]) == 1 for x in ("cnv", "sv"))

    if wgs_core:
        kind = WGS_KIND
        source_file = None
        assembly = "GRCh37"
        source = "Dante WGS"
        detail: dict[str, Any] = {}
    elif len(vcfs) == 1:
        source_file = vcfs[0]
        kind, detail = classify_single_vcf(source_file)
        assembly = detail.get("assembly", "UNKNOWN")
        source = detail.get("source", "UNKNOWN")
    else:
        kind = UNKNOWN_KIND
        source_file = None
        assembly = "UNKNOWN"
        source = "UNKNOWN"
        detail = {}

    if kind == ARRAY_KIND:
        capabilities = {
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
        }
        scope = (
            "Interrogated genotyping-array loci only. Homozygous-reference calls are valid "
            "only at loci present in the source VCF; absence elsewhere is unknown."
        )
    elif kind == WGS_KIND:
        capabilities = {
            "direct_rsid_export": False,
            "small_variant_annotation": True,
            "structural_variant_annotation": wgs_structural,
            "fastq_alignment": bool(fastqs),
            "reference_genome_completion": bool(fastqs),
            "telomere_from_fastq": bool(fastqs),
            "local_bam_browser": bool(fastqs),
            "knowledge_enrichment": wgs_structural,
            "final_reports": wgs_structural,
            "whole_genome_negative_evidence": True,
        }
        scope = "Dante WGS patient bundle; whole-genome analytical contracts apply where required inputs are present."
    else:
        capabilities = {
            "direct_rsid_export": False,
            "small_variant_annotation": False,
            "structural_variant_annotation": False,
            "fastq_alignment": bool(fastqs),
            "reference_genome_completion": False,
            "telomere_from_fastq": bool(fastqs),
            "local_bam_browser": False,
            "knowledge_enrichment": False,
            "final_reports": False,
            "whole_genome_negative_evidence": False,
        }
        scope = "Input type is not yet supported automatically. No whole-genome assumptions are made."

    return {
        "schema": SCHEMA,
        "kind": kind,
        "assembly": assembly,
        "source": source,
        "scope": scope,
        "capabilities": capabilities,
        "source_file": quick_identity(source_file) if source_file else None,
        "vcf_files": [quick_identity(x) for x in vcfs],
        "fastq_files": [quick_identity(x) for x in fastqs],
        "dante_roles": {k: [quick_identity(x) for x in v] for k, v in role_files.items()},
        "inspection": detail,
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
    """Attach patient input scope to completed Stage-24 outputs without changing findings."""
    patient_root = patient_root.resolve()
    profile_path = patient_root / "input_profile.json"
    if not profile_path.is_file():
        save_atomic(profile_path, detect(patient_root))
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    output_dir = patient_root / "Dante_results" / "knowledge_enrichment"
    state_path = output_dir / "Dante_knowledge_enrichment.state.json"
    report_path = output_dir / "Dante_knowledge_enrichment_report.txt"
    if not state_path.is_file() or not report_path.is_file():
        raise FileNotFoundError("Stage 24 state/report is missing")

    state = json.loads(state_path.read_text(encoding="utf-8"))
    scope_block = [
        "INPUT SCOPE",
        f"Input profile: {profile.get('kind', 'UNKNOWN')}",
        f"Source: {profile.get('source', 'UNKNOWN')}",
        f"Assembly: {profile.get('assembly', 'UNKNOWN')}",
        f"Scope: {profile.get('scope', '')}",
        "Whole-genome negative evidence: "
        + ("YES" if profile.get("capabilities", {}).get("whole_genome_negative_evidence") else "NO"),
    ]
    if profile.get("kind") == ARRAY_KIND:
        scope_block.extend([
            "Structural variants: NOT TESTED by this source; CNV/SV Stage-24 inputs are explicit N/A sentinels only.",
            "Technical confidence: array calls have no read-level DP/GQ/QUAL/FILTER contract in the source VCF.",
        ])
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
        "schema": profile.get("schema"),
        "kind": profile.get("kind"),
        "source": profile.get("source"),
        "assembly": profile.get("assembly"),
        "scope": profile.get("scope"),
        "whole_genome_negative_evidence": bool(
            profile.get("capabilities", {}).get("whole_genome_negative_evidence")
        ),
        "profile_sha256": sha256_file(profile_path),
    }
    outputs = state.get("outputs", {})
    for _, identity in outputs.items():
        if identity.get("name") == report_path.name:
            identity["length"] = report_path.stat().st_size
            identity["sha256"] = sha256_file(report_path)
            break
    save_atomic(state_path, state)
    return state


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("detect", "annotate-stage24", "self-test"))
    ap.add_argument("--patient-root", type=Path)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    if args.mode == "self-test":
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inp = root / "input"
            inp.mkdir()
            (inp / "sample.vcf").write_text(
                "##fileformat=VCFv4.2\n##source=Genotek\n##reference=hg19\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"
                "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0|0\n"
                "chr1\t200\trs2\tC\tT\t.\t.\t.\tGT\t0|1\n",
                encoding="utf-8",
            )
            result = detect(root)
            assert result["kind"] == ARRAY_KIND, result
            assert result["assembly"] == "GRCh37"
            assert result["capabilities"]["direct_rsid_export"]
            assert not result["capabilities"]["structural_variant_annotation"]
        print("SELF-TEST PASS | patient input profile detection")
        return 0

    if args.patient_root is None:
        ap.error("--patient-root is required")
    if args.mode == "annotate-stage24":
        payload = annotate_stage24(args.patient_root)
        if args.json:
            print(json.dumps(payload.get("input_profile", {}), ensure_ascii=False, sort_keys=True))
        else:
            print("Stage 24 input scope annotated")
        return 0

    payload = detect(args.patient_root)
    save_atomic(args.patient_root / "input_profile.json", payload)
    if args.json:
        print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    else:
        print(f"{payload['kind']} | {payload['assembly']} | {payload['source']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
