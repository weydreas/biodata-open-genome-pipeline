#!/usr/bin/env python3
"""Install and audit open reference-knowledge snapshots for the Dante WGS pipeline.

No third-party Python packages are required. Runtime downloads are versioned and
recorded in a reproducible manifest together with source URLs, hashes, licensing
classes and attribution text for later insertion into the final PDF report.

The installer intentionally does NOT download the 26/28 GiB full human dbSNP VCFs.
It downloads their tabix indexes and assembly reports and records the remote VCF
URLs for a later target-sliced normalization stage. This avoids a ~54 GiB local
copy while preserving an indexed, reproducible route to the public dbSNP data.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, asdict
from typing import Any

PIPELINE_RELEASE = "9.0.0"
BUNDLE_SCHEMA = "open-reference-bundle-v2"
USER_AGENT = f"DanteWGS/{PIPELINE_RELEASE} open-reference-manager"
DOWNLOAD_TIMEOUT = 120
MAX_ATTEMPTS = 5
CHUNK = 4 * 1024 * 1024

LICENSE_CLASSES = {
    "OPEN_NO_ATTRIBUTION",
    "OPEN_ATTRIBUTION_REQUIRED",
    "CUSTOM_OPEN_ATTRIBUTION",
}


@dataclass(frozen=True)
class Resource:
    source: str
    key: str
    url: str
    relative_path: str
    license_name: str
    license_class: str
    license_url: str
    attribution: str
    expected_md5: str | None = None
    content_kind: str = "file"
    notes: str = ""
    medical_disclaimer: str = ""
    static: bool = False


HGNC_ATTR = (
    "HGNC (HUGO Gene Nomenclature Committee), genenames.org; "
    "HGNC recommends attribution although its data are CC0."
)
GENCC_ATTR = (
    "Gene Curation Coalition (GenCC) gene-disease relationship data; "
    "GenCC requests attribution to GenCC and contributing sources."
)
CLINGEN_ATTR = (
    "Curated content obtained from the Clinical Genome Resource (ClinGen); "
    "include the date accessed in the final report."
)
GIAB_ATTR = (
    "Genome in a Bottle / NIST genome stratifications; acknowledge NIST as the source "
    "and preserve the snapshot/version in derived work."
)
MONDO_ATTR = "Mondo Disease Ontology, Monarch Initiative, CC BY 4.0; cite the version/date used."
HPO_ATTR = (
    "This service/product uses the Human Phenotype Ontology (HPO) (version information). "
    "Find out more at http://www.human-phenotype-ontology.org. "
    "The HPO Consortium must be acknowledged; publicly displayed HPO files must include "
    "their version/date and must not be altered."
)
DBSNP_ATTR = "NCBI dbSNP / RefSNP human variation data; record assembly and access date."


GIAB_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "genome-stratifications/v3.6/GRCh37%40all"
)

RESOURCES: tuple[Resource, ...] = (
    Resource(
        "HGNC", "complete", 
        "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
        "HGNC/hgnc_complete_set.tsv", "CC0-1.0", "OPEN_NO_ATTRIBUTION",
        "https://creativecommons.org/publicdomain/zero/1.0/", HGNC_ATTR,
        content_kind="tsv",
    ),
    Resource(
        "HGNC", "withdrawn",
        "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/withdrawn.txt",
        "HGNC/withdrawn.tsv", "CC0-1.0", "OPEN_NO_ATTRIBUTION",
        "https://creativecommons.org/publicdomain/zero/1.0/", HGNC_ATTR,
        content_kind="tsv",
    ),
    Resource(
        "GenCC", "submissions_new_tsv",
        "https://thegencc.org/download/action/submissions-export-tsv?format=new",
        "GenCC/gencc_submissions.tsv", "CC0-1.0", "OPEN_NO_ATTRIBUTION",
        "https://creativecommons.org/publicdomain/zero/1.0/", GENCC_ATTR,
        content_kind="tsv",
        medical_disclaimer=(
            "GenCC states its information is not intended for direct diagnostic use or medical "
            "decision-making without genetics-professional review and is provided as-is."
        ),
    ),
    Resource(
        "ClinGen", "gene_disease_validity",
        "https://search.clinicalgenome.org/kb/gene-validity/download",
        "ClinGen/gene_disease_validity.csv", "CC0-1.0", "OPEN_NO_ATTRIBUTION",
        "https://creativecommons.org/publicdomain/zero/1.0/", CLINGEN_ATTR,
        content_kind="csv",
    ),
    Resource(
        "ClinGen", "dosage_genes_regions",
        "https://search.clinicalgenome.org/kb/gene-dosage/downloadall",
        "ClinGen/dosage_genes_regions.csv", "CC0-1.0", "OPEN_NO_ATTRIBUTION",
        "https://creativecommons.org/publicdomain/zero/1.0/", CLINGEN_ATTR,
        content_kind="csv",
    ),
    Resource(
        "ClinGen", "actionability_adult",
        "https://actionability.clinicalgenome.org/ac/Adult/api/summ",
        "ClinGen/actionability_adult.json", "CC0-1.0", "OPEN_NO_ATTRIBUTION",
        "https://creativecommons.org/publicdomain/zero/1.0/", CLINGEN_ATTR,
        content_kind="json",
    ),
    Resource(
        "ClinGen", "actionability_pediatric",
        "https://actionability.clinicalgenome.org/ac/Pediatric/api/summ",
        "ClinGen/actionability_pediatric.json", "CC0-1.0", "OPEN_NO_ATTRIBUTION",
        "https://creativecommons.org/publicdomain/zero/1.0/", CLINGEN_ATTR,
        content_kind="json",
    ),
    Resource(
        "GIAB", "all_difficult",
        f"{GIAB_BASE}/Union/GRCh37_alldifficultregions.bed.gz",
        "GIAB/GRCh37_alldifficultregions.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "lowmap_segdup_union",
        f"{GIAB_BASE}/Union/GRCh37_alllowmapandsegdupregions.bed.gz",
        "GIAB/GRCh37_alllowmapandsegdupregions.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "low_mappability",
        f"{GIAB_BASE}/Mappability/GRCh37_lowmappabilityall.bed.gz",
        "GIAB/GRCh37_lowmappabilityall.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "segmental_duplications",
        f"{GIAB_BASE}/SegmentalDuplications/GRCh37_segdups.bed.gz",
        "GIAB/GRCh37_segdups.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "gaps_15kb",
        f"{GIAB_BASE}/OtherDifficult/GRCh37_gaps_slop15kb.bed.gz",
        "GIAB/GRCh37_gaps_slop15kb.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "mhc",
        f"{GIAB_BASE}/OtherDifficult/GRCh37_MHC.bed.gz",
        "GIAB/GRCh37_MHC.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "kir",
        f"{GIAB_BASE}/OtherDifficult/GRCh37_KIR.bed.gz",
        "GIAB/GRCh37_KIR.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "vdj",
        f"{GIAB_BASE}/OtherDifficult/GRCh37_VDJ.bed.gz",
        "GIAB/GRCh37_VDJ.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "chrX_PAR",
        f"{GIAB_BASE}/XY/GRCh37_chrX_PAR.bed.gz",
        "GIAB/GRCh37_chrX_PAR.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "chrY_PAR",
        f"{GIAB_BASE}/XY/GRCh37_chrY_PAR.bed.gz",
        "GIAB/GRCh37_chrY_PAR.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "chrX_XTR",
        f"{GIAB_BASE}/XY/GRCh37_chrX_XTR.bed.gz",
        "GIAB/GRCh37_chrX_XTR.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "GIAB", "chrY_XTR",
        f"{GIAB_BASE}/XY/GRCh37_chrY_XTR.bed.gz",
        "GIAB/GRCh37_chrY_XTR.bed.gz", "NIST-OPEN-DATA", "OPEN_ATTRIBUTION_REQUIRED",
        "https://www.nist.gov/open/license", GIAB_ATTR, content_kind="gzip", static=True,
    ),
    Resource(
        "Mondo", "ontology_obo",
        "https://purl.obolibrary.org/obo/mondo.obo",
        "Mondo/mondo.obo", "CC-BY-4.0", "OPEN_ATTRIBUTION_REQUIRED",
        "https://creativecommons.org/licenses/by/4.0/", MONDO_ATTR,
        content_kind="obo",
    ),
    Resource(
        "HPO", "ontology_obo",
        "https://purl.obolibrary.org/obo/hp.obo",
        "HPO/hp.obo", "HPO-CUSTOM-OPEN", "CUSTOM_OPEN_ATTRIBUTION",
        "https://hpo.jax.org/app/license", HPO_ATTR,
        content_kind="obo",
        notes=(
            "Ontology only. HPO disease/gene annotation products are intentionally excluded here "
            "because they aggregate third-party sources such as OMIM/Orphanet with separate terms."
        ),
    ),
    Resource(
        "dbSNP", "grch37_index",
        "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi",
        "dbSNP/GCF_000001405.25.gz.tbi", "NCBI-MOLECULAR-DATA", "OPEN_NO_ATTRIBUTION",
        "https://www.ncbi.nlm.nih.gov/home/about/policies/", DBSNP_ATTR,
        content_kind="binary",
    ),
    Resource(
        "dbSNP", "grch37_index_md5",
        "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi.md5",
        "dbSNP/GCF_000001405.25.gz.tbi.md5", "NCBI-MOLECULAR-DATA", "OPEN_NO_ATTRIBUTION",
        "https://www.ncbi.nlm.nih.gov/home/about/policies/", DBSNP_ATTR,
        content_kind="text",
    ),
    Resource(
        "dbSNP", "grch38_index",
        "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi",
        "dbSNP/GCF_000001405.40.gz.tbi", "NCBI-MOLECULAR-DATA", "OPEN_NO_ATTRIBUTION",
        "https://www.ncbi.nlm.nih.gov/home/about/policies/", DBSNP_ATTR,
        content_kind="binary",
    ),
    Resource(
        "dbSNP", "grch38_index_md5",
        "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi.md5",
        "dbSNP/GCF_000001405.40.gz.tbi.md5", "NCBI-MOLECULAR-DATA", "OPEN_NO_ATTRIBUTION",
        "https://www.ncbi.nlm.nih.gov/home/about/policies/", DBSNP_ATTR,
        content_kind="text",
    ),
    Resource(
        "dbSNP", "grch37_assembly_report",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt",
        "dbSNP/GCF_000001405.25_GRCh37.p13_assembly_report.txt", "NCBI-MOLECULAR-DATA", "OPEN_NO_ATTRIBUTION",
        "https://www.ncbi.nlm.nih.gov/home/about/policies/", DBSNP_ATTR,
        content_kind="text", static=True,
    ),
    Resource(
        "dbSNP", "grch38_assembly_report",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt",
        "dbSNP/GCF_000001405.40_GRCh38.p14_assembly_report.txt", "NCBI-MOLECULAR-DATA", "OPEN_NO_ATTRIBUTION",
        "https://www.ncbi.nlm.nih.gov/home/about/policies/", DBSNP_ATTR,
        content_kind="text", static=True,
    ),
)

DBSNP_REMOTE_VCFS = {
    "GRCh37.p13": "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz",
    "GRCh38.p14": "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz",
}


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def fmt_bytes(value: int | float | None) -> str:
    if value is None:
        return "?"
    x = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if x < 1024 or unit == "TiB":
            return f"{x:.1f} {unit}"
        x /= 1024
    return f"{x:.1f} TiB"


def fmt_eta(seconds: float | None) -> str:
    if seconds is None or seconds < 0 or seconds == float("inf"):
        return "?"
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def file_hash(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def request(url: str, *, method: str = "GET", headers: dict[str, str] | None = None):
    merged = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if headers:
        merged.update(headers)
    return urllib.request.Request(url, method=method, headers=merged)


def head_metadata(url: str) -> dict[str, Any]:
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request(url, method="HEAD"), timeout=30) as response:
                length = response.headers.get("Content-Length")
                return {
                    "url": url,
                    "final_url": response.geturl(),
                    "content_length": int(length) if length and length.isdigit() else None,
                    "etag": response.headers.get("ETag"),
                    "last_modified": response.headers.get("Last-Modified"),
                    "http_status": getattr(response, "status", 200),
                }
        except Exception as exc:
            if attempt == 3:
                return {"url": url, "head_error": f"{type(exc).__name__}: {exc}"}
            time.sleep(attempt * 2)
    return {"url": url, "head_error": "unknown"}


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def save_json_atomic(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(part, path)


def download_resource(
    resource: Resource,
    output_dir: Path,
    previous: dict[str, Any],
    force: bool,
) -> dict[str, Any]:
    dest = output_dir / resource.relative_path
    dest.parent.mkdir(parents=True, exist_ok=True)
    previous_art = previous.get("artifacts", {}).get(resource.key, {}) if previous else {}

    if dest.is_file() and not force:
        expected_len = previous_art.get("length")
        expected_sha = previous_art.get("sha256")
        if expected_len == dest.stat().st_size and expected_sha:
            print(f"[CURRENT] {resource.source}/{resource.key} | {fmt_bytes(dest.stat().st_size)}")
            return dict(previous_art)
        if resource.static and resource.expected_md5:
            print(f"[VERIFY]  {resource.source}/{resource.key} | checking pinned MD5")
            if file_hash(dest, "md5").lower() == resource.expected_md5.lower():
                return artifact_record(resource, dest, {}, reused=True)

    part = dest.with_name(dest.name + ".part")
    resume_from = part.stat().st_size if part.is_file() and not force else 0
    if force and part.exists():
        part.unlink()
        resume_from = 0

    conditional: dict[str, str] = {}
    if dest.is_file() and not force and previous_art:
        if previous_art.get("etag"):
            conditional["If-None-Match"] = str(previous_art["etag"])
        elif previous_art.get("last_modified"):
            conditional["If-Modified-Since"] = str(previous_art["last_modified"])

    if resume_from:
        conditional["Range"] = f"bytes={resume_from}-"

    last_error: Exception | None = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            req = request(resource.url, headers=conditional)
            try:
                response = urllib.request.urlopen(req, timeout=DOWNLOAD_TIMEOUT)
            except urllib.error.HTTPError as exc:
                if exc.code == 304 and dest.is_file():
                    print(f"[CURRENT] {resource.source}/{resource.key} | remote unchanged")
                    return artifact_record(resource, dest, dict(exc.headers), reused=True)
                if exc.code == 416 and part.is_file():
                    part.unlink()
                    resume_from = 0
                    conditional.pop("Range", None)
                    continue
                raise

            with response:
                status = getattr(response, "status", 200)
                if resume_from and status != 206:
                    resume_from = 0
                mode = "ab" if resume_from and status == 206 else "wb"
                if mode == "wb" and part.exists():
                    part.unlink()
                content_length = response.headers.get("Content-Length")
                remaining = int(content_length) if content_length and content_length.isdigit() else None
                total = (resume_from + remaining) if remaining is not None else None
                downloaded = resume_from
                start = time.monotonic()
                last_report = start
                with part.open(mode) as out:
                    while True:
                        block = response.read(CHUNK)
                        if not block:
                            break
                        out.write(block)
                        downloaded += len(block)
                        now = time.monotonic()
                        if now - last_report >= 2.0:
                            elapsed = max(0.001, now - start)
                            speed = (downloaded - resume_from) / elapsed
                            eta = ((total - downloaded) / speed) if total and speed > 0 else None
                            pct = (100.0 * downloaded / total) if total else None
                            pct_text = f"{pct:6.2f}%" if pct is not None else "   ?  "
                            print(
                                f"[DOWNLOAD] {resource.source}/{resource.key} | {pct_text} | "
                                f"{fmt_bytes(downloaded)} / {fmt_bytes(total)} | "
                                f"{fmt_bytes(speed)}/s | ETA {fmt_eta(eta)}",
                                flush=True,
                            )
                            last_report = now
                headers = dict(response.headers.items())
            if part.stat().st_size <= 0:
                raise RuntimeError("download produced an empty file")
            os.replace(part, dest)
            validate_resource(resource, dest)
            if resource.expected_md5:
                actual_md5 = file_hash(dest, "md5")
                if actual_md5.lower() != resource.expected_md5.lower():
                    raise RuntimeError(
                        f"MD5 mismatch for {resource.source}/{resource.key}: "
                        f"expected {resource.expected_md5}, got {actual_md5}"
                    )
            record = artifact_record(resource, dest, headers, reused=False)
            print(
                f"[SAVED]   {resource.source}/{resource.key} | {fmt_bytes(dest.stat().st_size)} | "
                f"sha256={record['sha256'][:16]}...",
                flush=True,
            )
            return record
        except Exception as exc:
            last_error = exc
            print(
                f"[RETRY]   {resource.source}/{resource.key} | attempt {attempt}/{MAX_ATTEMPTS} | "
                f"{type(exc).__name__}: {exc}",
                file=sys.stderr,
                flush=True,
            )
            if attempt < MAX_ATTEMPTS:
                time.sleep(min(30, 2 ** attempt))
    raise RuntimeError(f"Failed to download {resource.source}/{resource.key}: {last_error}")


def artifact_record(resource: Resource, dest: Path, headers: dict[str, Any], reused: bool) -> dict[str, Any]:
    return {
        "path": resource.relative_path.replace("\\", "/"),
        "url": resource.url,
        "length": dest.stat().st_size,
        "sha256": file_hash(dest),
        "md5": file_hash(dest, "md5") if resource.expected_md5 else None,
        "etag": headers.get("ETag") or headers.get("Etag"),
        "last_modified": headers.get("Last-Modified"),
        "downloaded_or_verified_utc": utcnow(),
        "reused": bool(reused),
    }


def validate_resource(resource: Resource, path: Path) -> None:
    if path.stat().st_size <= 0:
        raise RuntimeError(f"Empty resource: {path}")
    kind = resource.content_kind
    if kind == "gzip":
        with gzip.open(path, "rb") as handle:
            if not handle.read(64):
                raise RuntimeError(f"Empty gzip payload: {path}")
    elif kind == "json":
        with path.open("rt", encoding="utf-8") as handle:
            json.load(handle)
    elif kind == "obo":
        text = path.read_text(encoding="utf-8", errors="strict")[:4096]
        if "format-version:" not in text and "ontology:" not in text:
            raise RuntimeError(f"File does not look like an OBO ontology: {path}")
    elif kind in {"tsv", "csv", "text"}:
        text = path.read_text(encoding="utf-8", errors="replace")[:4096]
        if not text.strip():
            raise RuntimeError(f"Empty text payload: {path}")
        if resource.source == "GenCC" and "sgc_id" not in text.lower():
            raise RuntimeError("GenCC download is not the expected new SGC-ID format")


def source_group(resource: Resource) -> dict[str, Any]:
    return {
        "license_name": resource.license_name,
        "license_class": resource.license_class,
        "license_url": resource.license_url,
        "attribution": resource.attribution,
        "notes": resource.notes,
        "medical_disclaimer": resource.medical_disclaimer,
    }


def write_attributions(output_dir: Path, manifest: dict[str, Any]) -> None:
    lines = [
        "# Open reference-data attribution ledger",
        "",
        "Generated automatically by the Dante WGS pipeline. The future PDF report generator should "
        "include the applicable entries for every source actually used in a published report.",
        "",
        f"Snapshot installed: `{manifest['installed_utc']}`",
        "",
    ]
    for source in sorted(manifest["sources"]):
        meta = manifest["sources"][source]
        lines += [
            f"## {source}",
            "",
            f"- License: `{meta['license_name']}`",
            f"- License class: `{meta['license_class']}`",
            f"- License URL: {meta['license_url']}",
            f"- Attribution: {meta['attribution']}",
        ]
        if meta.get("notes"):
            lines.append(f"- Notes: {meta['notes']}")
        if meta.get("medical_disclaimer"):
            lines.append(f"- Source disclaimer: {meta['medical_disclaimer']}")
        lines.append("")
    lines += [
        "## dbSNP storage policy",
        "",
        "The full human dbSNP VCFs are intentionally **not** copied into this project. Their current "
        "GRCh37 and GRCh38 remote URLs, content metadata and local tabix indexes are recorded in the "
        "bundle manifest. A later normalization stage can query only target regions/records.",
        "",
        "## HPO scope",
        "",
        "Only the HPO ontology itself (`hp.obo`) is installed here. HPO disease/gene annotation "
        "products are excluded from this bundle because they aggregate third-party sources with "
        "separate licensing/provenance requirements.",
        "",
    ]
    (output_dir / "ATTRIBUTIONS.md").write_text("\n".join(lines), encoding="utf-8")


def install(root: Path, output_dir: Path, bcftools_image: str, force: bool) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "open_reference_manifest.json"
    previous = load_json(manifest_path)

    free = shutil.disk_usage(output_dir).free
    # Open bundle excludes non-commercial AlphaMissense; 2 GiB is ample atomic headroom.
    if free < 2 * 1024**3:
        raise RuntimeError(
            f"Open-reference install requires at least 2 GiB free for safe atomic downloads; only {fmt_bytes(free)} available."
        )

    print("OPEN REFERENCE KNOWLEDGE BUNDLE")
    print(f"  Pipeline release:      v{PIPELINE_RELEASE}")
    print(f"  Output:                {output_dir}")
    print(f"  Free disk:             {fmt_bytes(free)}")
    print("  Full dbSNP VCFs:       remote-indexed only (26 GiB GRCh37 + 28 GiB GRCh38 are NOT downloaded)")
    print("  HPO scope:             ontology only; third-party-derived disease annotations excluded")
    print("  Python dependencies:   standard library only")
    print("")

    grouped: dict[str, dict[str, Any]] = {}
    for resource in RESOURCES:
        if resource.license_class not in LICENSE_CLASSES:
            raise RuntimeError(f"Blocked license class in installer catalog: {resource.source} {resource.license_class}")
        group = grouped.setdefault(resource.source, source_group(resource) | {"artifacts": {}})
        prev_source = previous.get("sources", {}).get(resource.source, {}) if previous else {}
        rec = download_resource(resource, output_dir, prev_source, force)
        group["artifacts"][resource.key] = rec

    remote_dbsnp = {assembly: head_metadata(url) for assembly, url in DBSNP_REMOTE_VCFS.items()}
    remote_manifest = output_dir / "dbSNP" / "remote_vcf_sources.json"
    save_json_atomic(remote_manifest, remote_dbsnp)
    grouped["dbSNP"]["artifacts"]["remote_vcf_sources"] = {
        "path": str(remote_manifest.relative_to(output_dir)).replace("\\", "/"),
        "length": remote_manifest.stat().st_size,
        "sha256": file_hash(remote_manifest),
        "created_utc": utcnow(),
    }

    vep_cache_candidates = [
        root / "annotation_data" / "vep_rsid",
        root / "annotation_data" / "vep_rsid" / "homo_sapiens",
    ]
    vep_present = any(p.exists() for p in vep_cache_candidates)

    manifest = {
        "schema": BUNDLE_SCHEMA,
        "status": "complete",
        "pipeline_release": PIPELINE_RELEASE,
        "installed_utc": utcnow(),
        "policy": {
            "allowed_license_classes": sorted(LICENSE_CLASSES),
            "blocked_by_design": [
                "non-commercial-only datasets",
                "research-only datasets",
                "PanelApp custom commercial/healthcare restrictions",
                "OMIM licensed content",
                "HPO annotation products carrying third-party provenance",
            ],
            "pdf_attribution_ledger": "ATTRIBUTIONS.md",
        },
        "sources": grouped,
        "existing_vep_cache": {
            "present": bool(vep_present),
            "expected_embedded_annotations": [
                "gnomAD allele frequencies",
                "1000 Genomes frequencies",
                "SIFT",
                "PolyPhen",
            ],
            "note": "No duplicate population/predictor bundle is downloaded by this stage; Stage 9 VEP is the canonical source.",
        },
        "dbsnp_strategy": {
            "full_vcfs_downloaded": False,
            "reason": "Avoid ~54 GiB duplicate local storage; keep local tabix indexes + assembly reports and remote indexed VCF URLs.",
            "future_target_map": "A later deterministic enrichment stage will build a compact target-only dual-build map from these resources.",
        },
    }
    write_attributions(output_dir, manifest)
    manifest["attributions"] = {
        "path": "ATTRIBUTIONS.md",
        "length": (output_dir / "ATTRIBUTIONS.md").stat().st_size,
        "sha256": file_hash(output_dir / "ATTRIBUTIONS.md"),
    }
    save_json_atomic(manifest_path, manifest)

    print("")
    print("[DONE] Open reference knowledge bundle installed")
    print(f"       Manifest:      {manifest_path}")
    print(f"       Attributions:  {output_dir / 'ATTRIBUTIONS.md'}")
    print(f"       Sources:       {', '.join(sorted(grouped))}")
    if not vep_present:
        print("[NOTE] VEP cache was not detected by the lightweight path probe; item 1/9 remains the canonical VEP setup.")
    return 0


def status(output_dir: Path) -> int:
    manifest_path = output_dir / "open_reference_manifest.json"
    manifest = load_json(manifest_path)
    if not manifest or manifest.get("schema") != BUNDLE_SCHEMA or manifest.get("status") != "complete":
        print("NOT RUN" if not manifest else "INCOMPLETE")
        return 1
    missing = []
    for source, meta in manifest.get("sources", {}).items():
        for key, art in meta.get("artifacts", {}).items():
            rel = art.get("path")
            if not rel:
                continue
            p = output_dir / rel
            if not p.is_file() or p.stat().st_size != art.get("length"):
                missing.append(f"{source}/{key}")
    attr = manifest.get("attributions", {}).get("path")
    if attr and not (output_dir / attr).is_file():
        missing.append("ATTRIBUTIONS.md")
    if missing:
        print("INCOMPLETE: " + ", ".join(missing[:8]))
        return 1
    print("READY")
    return 0


def self_test() -> int:
    errors: list[str] = []
    if len({(r.source, r.key) for r in RESOURCES}) != len(RESOURCES):
        errors.append("duplicate source/key in resource catalog")
    for r in RESOURCES:
        if r.license_class not in LICENSE_CLASSES:
            errors.append(f"blocked license class: {r.source}/{r.key} {r.license_class}")
        policy_text = " ".join([r.license_name, r.license_class, r.notes]).upper()
        for blocked in ("CC-BY-NC", "NON-COMMERCIAL", "RESEARCH-ONLY", "PANELAPP", "OMIM LICENSE"):
            if blocked in policy_text:
                errors.append(f"blocked policy leaked into catalog: {r.source}/{r.key}: {blocked}")
        if r.source == "HPO" and r.key != "ontology_obo":
            errors.append("HPO annotations must not be bundled")
        if r.source == "dbSNP" and r.relative_path.endswith(".gz") and not r.relative_path.endswith(".tbi.gz"):
            # Full dbSNP VCF names end .gz; indexes are .gz.tbi and therefore do not match.
            errors.append(f"full dbSNP VCF must not be bundled: {r.relative_path}")
    if any(r.source == "AlphaMissense" for r in RESOURCES):
        errors.append("AlphaMissense is CC BY-NC-SA 4.0 and must not appear in the open bundle")
    required_sources = {"HGNC", "dbSNP", "GIAB", "GenCC", "ClinGen", "Mondo", "HPO"}
    found_sources = {r.source for r in RESOURCES}
    missing = required_sources - found_sources
    if missing:
        errors.append("missing required sources: " + ", ".join(sorted(missing)))
    if errors:
        for error in errors:
            print("SELF-TEST ERROR:", error)
        return 1
    print(f"SELF-TEST PASS | {len(RESOURCES)} download artifacts | sources={','.join(sorted(found_sources))}")
    print("Policy: commercial-safe/open attribution sources only; AlphaMissense excluded; HPO ontology-only; full dbSNP VCFs remote-only")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("install", "status", "self-test"))
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bcftools-image", default="quay.io/biocontainers/bcftools:1.24--h487d631_1")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    output = (args.output or (root / "reference_data" / "open_reference")).resolve()
    if args.mode == "self-test":
        return self_test()
    if args.mode == "status":
        return status(output)
    return install(root, output, args.bcftools_image, args.force)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted; completed downloads and .part files are preserved for resume.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
