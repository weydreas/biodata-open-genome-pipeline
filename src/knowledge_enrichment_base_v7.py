#!/usr/bin/env python3
"""Full-genome deterministic knowledge enrichment for Dante WGS v9.0.0.

Stage 24 scans the canonical GRCh37 genomic outputs and applies the installed
open and research-only reference bundles without using OpenAI.  It deliberately
keeps commercial/open findings physically separate from research-only overlays.

Design boundaries:
- every row of the unified small-variant catalog is scanned;
- only knowledge-relevant patient events are materialized in the patient SQLite
  database (the full catalog remains the canonical lossless allele catalog);
- open/reference evidence never becomes research-only merely because a CADD or
  PanelApp overlay also exists;
- AlphaMissense is supporting prediction evidence, never a pathogenicity label;
- CADD and PanelApp are research-only and are never written into the open TSV;
- GIAB difficult-region annotations change technical confidence, not biological
  classification;
- inheritance-aware flags are hypothesis labels only; unphased variants are
  never called compound heterozygous.

Only the Python standard library is required.  Indexed lookups use the pipeline's
already-approved bcftools/htslib Docker image via its `tabix` / `bcftools`
entrypoints.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
import traceback
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Iterable, Iterator, TextIO

PIPELINE_RELEASE = "9.0.0"
ENRICHMENT_SCHEMA = "dante-knowledge-enrichment-v1"
INDEX_SCHEMA = "dante-knowledge-reference-index-v1"
RULES_SCHEMA = "dante-knowledge-rules-v7"
USAGE_SCHEMA = "knowledge-usage-v1"
ASSEMBLY = "GRCh37"
DEFAULT_MAX_AF = 0.01

OPEN_FINDING_COLUMNS = [
    "finding_id", "priority_rank", "evidence_rank", "patient_relevance", "reason_code", "evidence_class",
    "event_type", "event_key", "source_layers", "chrom", "position", "end",
    "ref", "alt", "genotype", "zygosity", "gene", "hgnc_id",
    "consequence", "impact", "max_af", "rsid", "resolved_rsid",
    "clinvar_significance", "clinvar_review_stars", "clinvar_conditions",
    "disease_id", "disease_name", "gene_disease_classification", "inheritance",
    "actionability", "dosage_evidence", "possible_biallelic_unphased", "giab_flags",
    "technical_confidence", "evidence_sources", "interpretation_boundary",
    "provenance_json",
]

RESEARCH_COLUMNS = [
    "overlay_id", "base_open_finding_id", "priority_rank", "reason_code",
    "event_type", "event_key", "chrom", "position", "end", "ref", "alt",
    "gene", "consequence", "max_af", "panelapp_panels", "panelapp_confidence",
    "panelapp_inheritance", "alpha_missense_pathogenicity", "alpha_missense_class",
    "cadd_raw", "cadd_phred", "giab_flags",
    "technical_confidence", "evidence_sources", "research_banner",
    "interpretation_boundary", "provenance_json",
]

QC_COLUMNS = ["event_key", "event_type", "chrom", "position", "end", "flags", "technical_confidence"]

LOF_TERMS = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "start_lost",
}
MISSENSE_TERMS = {"missense_variant"}
FUNCTIONAL_TERMS = LOF_TERMS | MISSENSE_TERMS | {
    "stop_lost", "inframe_insertion", "inframe_deletion", "protein_altering_variant",
    "splice_region_variant",
}
STRONG_VALIDITY = {"definitive", "strong"}


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def norm_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (value or "").lower())


def normalize_chrom(value: str) -> str:
    value = (value or "").strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value == "M":
        value = "MT"
    return value


def canonical_chrom(value: str) -> bool:
    return normalize_chrom(value) in {*(str(i) for i in range(1, 23)), "X", "Y", "MT"}


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {"name": path.name, "length": stat.st_size, "sha256": sha256_file(path)}


def quick_identity(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "name": path.name,
        "length": stat.st_size,
        "modified_utc_ticks": int(stat.st_mtime_ns // 100),
    }


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except Exception:
        return {}


def save_json_atomic(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(part, path)


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("rt", encoding="utf-8", errors="replace", newline="")


@contextmanager
def atomic_text(path: Path, gzip_output: bool = False) -> Iterator[TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    raw = None
    gz = None
    text = None
    try:
        if gzip_output or path.name.lower().endswith(".gz"):
            raw = part.open("wb")
            gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
            text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
        else:
            text = part.open("wt", encoding="utf-8", newline="")
        yield text
        text.close()
        text = None
        if gz is not None and not gz.closed:
            gz.close()
        if raw is not None and not raw.closed:
            raw.close()
        os.replace(part, path)
    except Exception:
        try:
            if text is not None and not text.closed:
                text.close()
            if gz is not None and not gz.closed:
                gz.close()
            if raw is not None and not raw.closed:
                raw.close()
        finally:
            part.unlink(missing_ok=True)
        raise


def split_multi(value: str) -> list[str]:
    value = (value or "").strip().strip('"')
    if not value or value in {".", "-"}:
        return []
    values = re.split(r"\s*[|;,]\s*", value)
    return [v.strip().strip('"') for v in values if v.strip().strip('"')]


def split_consequences(value: str) -> set[str]:
    if not value or value in {".", "-"}:
        return set()
    return {v.strip().lower() for v in re.split(r"[&,;]", value) if v.strip()}


def first_value(row: dict[str, str], names: Iterable[str], contains: Iterable[str] = ()) -> str:
    normalized = {norm_header(k): (v or "") for k, v in row.items()}
    for name in names:
        key = norm_header(name)
        if key in normalized and normalized[key].strip():
            return normalized[key].strip()
    contains_norm = [norm_header(x) for x in contains]
    for key, value in normalized.items():
        if value.strip() and any(token in key for token in contains_norm):
            return value.strip()
    return ""


def number(value: str | None) -> float | None:
    if value is None:
        return None
    values: list[float] = []
    for token in re.split(r"[&,;|]", str(value)):
        token = token.strip()
        if not token or token in {".", "-", "NA", "na"}:
            continue
        try:
            values.append(float(token))
        except ValueError:
            pass
    return max(values) if values else None


def rarity_status(value: str, threshold: float) -> str:
    af = number(value)
    if af is None:
        return "UNKNOWN"
    return "RARE" if af <= threshold else "COMMON"


def is_rare(value: str, threshold: float) -> bool:
    """True only when rarity is positively supported by a numeric population AF.

    Missing AF is uncertainty, not evidence of rarity.  Stage 24 preserves such
    variants in explicit AF_UNKNOWN candidate classes instead of promoting them
    under a reason code whose name claims they are rare.
    """
    return rarity_status(value, threshold) == "RARE"


def review_stars(status: str) -> int:
    normalized = (status or "").lower().replace("_", " ")
    if "practice guideline" in normalized:
        return 4
    if "reviewed by expert panel" in normalized:
        return 3
    if "multiple submitters" in normalized and "no conflicts" in normalized:
        return 2
    if "criteria provided" in normalized:
        return 1
    return 0


def clinvar_bucket(significance: str) -> str:
    """Conservatively normalize aggregate ClinVar CLNSIG text.

    Aggregate CLNSIG can combine interpretations for different conditions.
    Therefore pathogenic + protective/risk/drug-response/association is not a
    clean P/LP assertion and is kept as MIXED rather than promoted as P/LP.
    """
    text = (significance or "").strip().lower().replace("_", " ")
    if not text.strip(" ."):
        return "NONE"
    if "conflicting" in text:
        return "CONFLICT"
    if "association not found" in text:
        return "OTHER"

    tokens = [re.sub(r"\s+", " ", t.strip()) for t in re.split(r"[|/,;]+", text) if t.strip()]
    path_tokens = {"pathogenic", "likely pathogenic"}
    benign_tokens = {"benign", "likely benign"}
    association_tokens = {"drug response", "risk factor", "protective", "association", "affects"}
    pathogenic = any(t in path_tokens for t in tokens)
    benign = any(t in benign_tokens for t in tokens)
    association = any(t in association_tokens for t in tokens)
    if pathogenic and benign:
        return "CONFLICT"
    if pathogenic and association:
        return "MIXED"
    if pathogenic:
        return "PATHOGENIC"
    if association:
        return "ASSOCIATION"
    if any(("uncertain" in t) or t == "vus" for t in tokens):
        return "VUS"
    if benign:
        return "BENIGN"
    return "OTHER"


def zygosity(gt: str, ref: str, alt: str) -> str:
    gt = (gt or "").replace("|", "/")
    if gt in {"", ".", "./."}:
        return "UNKNOWN"
    if re.fullmatch(r"\d+(?:/\d+)?", gt):
        nums = [int(x) for x in gt.split("/")]
        if len(nums) == 1:
            return "HEMIZYGOUS_ALT" if nums[0] > 0 else "HEMIZYGOUS_REF"
        if all(x == 0 for x in nums):
            return "HOM_REF"
        if len(set(nums)) == 1 and nums[0] > 0:
            return "HOM_ALT"
        if any(x > 0 for x in nums):
            return "HET_ALT"
    # Stage-17 genotype may be bases rather than numeric GT.
    bases = re.sub(r"[/|]", "", gt.upper())
    if bases and set(bases) <= {"A", "C", "G", "T"} and len(ref) == len(alt) == 1:
        alt_count = bases.count(alt.upper())
        ref_count = bases.count(ref.upper())
        if alt_count == len(bases) and alt_count > 0:
            return "HOM_ALT" if len(bases) > 1 else "HEMIZYGOUS_ALT"
        if alt_count > 0 and ref_count > 0:
            return "HET_ALT"
    return "UNKNOWN"


def canonical_event_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"GRCh37:{normalize_chrom(chrom)}:{int(pos)}:{ref}:{alt}"


def deterministic_id(prefix: str, *parts: str) -> str:
    digest = hashlib.sha256("\x1f".join(str(p) for p in parts).encode("utf-8")).hexdigest()[:20]
    return f"{prefix}_{digest}"


def classification_rank(value: str) -> int:
    text = (value or "").lower()
    if "definitive" in text:
        return 6
    if re.search(r"\bstrong\b", text):
        return 5
    if "moderate" in text:
        return 4
    if "limited" in text:
        return 3
    if "disputed" in text:
        return 1
    if "refuted" in text:
        return 0
    return 2


def normalized_classification(value: str) -> str:
    value = (value or "").strip()
    if not value:
        return "UNKNOWN"
    for label in ("Definitive", "Strong", "Moderate", "Limited", "Disputed", "Refuted"):
        if label.lower() in value.lower():
            return label.upper()
    return value.upper().replace(" ", "_")


def inheritance_bucket(value: str) -> str:
    """Normalize mode-of-inheritance text to a small idempotent vocabulary.

    The function is intentionally idempotent because collapsed GeneDisease rows
    store the canonical bucket (for example ``X_LINKED``) and may be normalized
    again later while findings are built.  Hotfix1 did not recognize its own
    canonical ``X_LINKED`` token, which could turn a specific X-linked assertion
    back into ``UNKNOWN`` and collide with a genuinely unknown assertion.
    """
    raw = (value or "").strip()
    canonical = re.sub(r"[\s-]+", "_", raw.upper())
    if canonical in {"AR", "AD", "X_LINKED", "MITOCHONDRIAL", "UNKNOWN"}:
        return canonical

    text = raw.lower()
    # Check X-linked before dominant/recessive so phrases such as
    # "X-linked recessive" stay in the X-linked bucket rather than AR.
    if "x-linked" in text or "x linked" in text or "x_linked" in text:
        return "X_LINKED"
    if "mitochond" in text or "maternal inheritance" in text:
        return "MITOCHONDRIAL"
    if "recessive" in text or re.search(r"\bar\b", text):
        return "AR"
    if "dominant" in text or re.search(r"\bad\b", text):
        return "AD"
    return "UNKNOWN"


@dataclass
class GeneDisease:
    gene: str
    hgnc_id: str
    disease_id: str
    disease_name: str
    classification: str
    inheritance: str
    source: str


@dataclass(frozen=True)
class DosageRegion:
    region_id: str
    name: str
    chrom: str
    start: int
    end: int
    haploinsufficiency: str
    triplosensitivity: str

    @property
    def evidence(self) -> str:
        return (
            f"{self.region_id}|{self.name}|GRCh37:{self.chrom}:{self.start}-{self.end}|"
            f"HI={self.haploinsufficiency or '.'};TS={self.triplosensitivity or '.'}"
        )


@dataclass
class GeneProfile:
    gene: str
    hgnc_id: str = ""
    relations: list[GeneDisease] = field(default_factory=list)
    actionable: bool = False
    actionability_context: list[str] = field(default_factory=list)
    dosage: list[str] = field(default_factory=list)
    panelapp: list[dict[str, str]] = field(default_factory=list)

    @property
    def strong_relations(self) -> list[GeneDisease]:
        return [r for r in self.relations if normalized_classification(r.classification).lower() in STRONG_VALIDITY]


class GeneKnowledge:
    def __init__(self, db: Path) -> None:
        self.db = db
        self.alias_to_gene: dict[str, tuple[str, str]] = {}
        self.profiles: dict[str, GeneProfile] = {}
        self.dosage_regions: dict[str, list[DosageRegion]] = defaultdict(list)
        con = sqlite3.connect(db)
        try:
            for alias, symbol, hgnc_id in con.execute("SELECT alias,symbol,hgnc_id FROM gene_alias"):
                self.alias_to_gene[str(alias).upper()] = (str(symbol), str(hgnc_id or ""))
            for symbol, hgnc_id in con.execute("SELECT symbol,hgnc_id FROM genes"):
                self.profiles[str(symbol)] = GeneProfile(str(symbol), str(hgnc_id or ""))
            for row in con.execute(
                "SELECT gene,hgnc_id,disease_id,disease_name,classification,inheritance,source FROM gene_disease"
            ):
                rel = GeneDisease(*(str(v or "") for v in row))
                self.profiles.setdefault(rel.gene, GeneProfile(rel.gene, rel.hgnc_id)).relations.append(rel)
            for gene, context in con.execute("SELECT gene,context FROM actionability"):
                p = self.profiles.setdefault(str(gene), GeneProfile(str(gene)))
                p.actionable = True
                if context:
                    p.actionability_context.append(str(context))
            for gene, evidence in con.execute("SELECT gene,evidence FROM dosage"):
                self.profiles.setdefault(str(gene), GeneProfile(str(gene))).dosage.append(str(evidence or ""))
            try:
                for region_id, name, chrom, start, end, hi, ts in con.execute(
                    "SELECT region_id,name,chrom,start,end,haploinsufficiency,triplosensitivity FROM dosage_region"
                ):
                    region = DosageRegion(
                        str(region_id or ""), str(name or ""), normalize_chrom(str(chrom or "")),
                        int(start), int(end), str(hi or ""), str(ts or ""),
                    )
                    self.dosage_regions[region.chrom].append(region)
                for chrom in self.dosage_regions:
                    self.dosage_regions[chrom].sort(key=lambda r: (r.start, r.end, r.region_id))
            except sqlite3.OperationalError:
                # Backward-safe while rebuilding a pre-v9 shared index; the state contract
                # normally causes an immediate rebuild before this class is instantiated.
                pass
            for gene, panel_id, panel_name, confidence, inheritance, phenotype in con.execute(
                "SELECT gene,panel_id,panel_name,confidence,inheritance,phenotype FROM panelapp_gene"
            ):
                self.profiles.setdefault(str(gene), GeneProfile(str(gene))).panelapp.append({
                    "panel_id": str(panel_id), "panel_name": str(panel_name or ""),
                    "confidence": str(confidence or ""), "inheritance": str(inheritance or ""),
                    "phenotype": str(phenotype or ""),
                })
        finally:
            con.close()

    def canonicalize(self, symbol: str) -> tuple[str, str]:
        token = (symbol or "").strip()
        if not token or token in {".", "-"}:
            return "", ""
        return self.alias_to_gene.get(token.upper(), (token, ""))

    def profile(self, symbol: str) -> GeneProfile | None:
        canonical, hgnc_id = self.canonicalize(symbol)
        if not canonical:
            return None
        p = self.profiles.get(canonical)
        if p is None:
            return GeneProfile(canonical, hgnc_id)
        return p

    def overlapping_dosage_regions(self, chrom: str, start: int, end: int) -> list[DosageRegion]:
        # The ClinGen dosage region set is small enough that a chromosome-local scan
        # is cheap relative to AnnotSV I/O and avoids a second dependency/index format.
        lo, hi = min(start, end), max(start, end)
        return [r for r in self.dosage_regions.get(normalize_chrom(chrom), []) if r.start <= hi and r.end >= lo]


# ----------------------------- reference-index build -----------------------------

def source_artifact(manifest: dict[str, Any], root: Path, source: str, key: str) -> Path | None:
    try:
        rel = manifest["sources"][source]["artifacts"][key]["path"]
    except Exception:
        return None
    path = root / str(rel)
    return path if path.is_file() else None


def parse_obo(path: Path) -> Iterator[tuple[str, str]]:
    current_id = ""
    current_name = ""
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                if current_id and current_name:
                    yield current_id, current_name
                current_id = current_name = ""
            elif line.startswith("id: "):
                current_id = line[4:].strip()
            elif line.startswith("name: "):
                current_name = line[6:].strip()
        if current_id and current_name:
            yield current_id, current_name


def _read_dict_rows(path: Path, delimiter: str | None = None) -> Iterator[dict[str, str]]:
    with open_text(path) as handle:
        # Consume comment lines while retaining the actual header.
        first = ""
        while True:
            first = handle.readline()
            if not first:
                return
            if first.strip() and not first.startswith("#"):
                break
        if delimiter is None:
            delimiter = "\t" if first.count("\t") >= first.count(",") else ","
        remainder = first + handle.read()
        reader = csv.DictReader(io.StringIO(remainder), delimiter=delimiter)
        for row in reader:
            yield {str(k or ""): str(v or "") for k, v in row.items()}


def _read_clingen_rows(path: Path, required_header_tokens: Iterable[str]) -> Iterator[dict[str, str]]:
    """Read ClinGen CSV downloads that contain a human-readable preamble.

    Current ClinGen exports start with title/date/webpage/separator rows before the
    actual CSV header. A normal DictReader would silently treat the title as the
    header and therefore index zero useful records.
    """
    with open_text(path) as handle:
        lines = handle.readlines()
    required = [norm_header(x) for x in required_header_tokens]
    header_index = None
    for i, line in enumerate(lines):
        parsed = next(csv.reader([line]))
        normalized = [norm_header(x) for x in parsed]
        joined = " ".join(normalized)
        if all(token in joined for token in required):
            header_index = i
            break
    if header_index is None:
        raise RuntimeError(
            f"ClinGen table header not found in {path.name}; expected tokens: {', '.join(required_header_tokens)}"
        )
    reader = csv.DictReader(io.StringIO("".join(lines[header_index:])))
    for row in reader:
        values = [str(v or "").strip() for v in row.values()]
        if not any(values):
            continue
        if values and values[0] and set(values[0]) <= {"+"}:
            continue
        yield {str(k or ""): str(v or "") for k, v in row.items()}


def parse_grch37_interval(value: str) -> tuple[str, int, int] | None:
    match = re.fullmatch(r"chr([0-9XYM]+):(\d+)-(\d+)", (value or "").strip(), flags=re.I)
    if not match:
        return None
    chrom = normalize_chrom(match.group(1))
    if not canonical_chrom(chrom):
        return None
    start, end = int(match.group(2)), int(match.group(3))
    return chrom, min(start, end), max(start, end)


def dosage_is_sufficient(value: str) -> bool:
    text = (value or "").strip().lower()
    return "sufficient evidence" in text or text in {"3", "score 3", "3.0"}


def annotsv_class_number(value: str) -> int | None:
    """Extract ACMG/AnnotSV class 1..5 from values such as '5', 'full=5', 'NA'."""
    match = re.search(r"(?<!\d)([1-5])(?!\d)", str(value or ""))
    return int(match.group(1)) if match else None


def dosage_evidence_supports(evidence: str, sv_type: str) -> bool:
    """Return whether ClinGen dosage evidence supports the event copy-number direction."""
    sv = (sv_type or "").upper()
    text = evidence or ""
    hi_match = re.search(r"HI=([^;|]+)", text, flags=re.I)
    ts_match = re.search(r"TS=([^;|]+)", text, flags=re.I)
    hi = hi_match.group(1).strip() if hi_match else ""
    ts = ts_match.group(1).strip() if ts_match else ""
    if sv == "DEL":
        return dosage_is_sufficient(hi)
    if sv == "DUP":
        return dosage_is_sufficient(ts)
    return False


def dosage_region_supports(region: DosageRegion, sv_type: str) -> bool:
    sv = (sv_type or "").upper()
    if sv == "DEL":
        return dosage_is_sufficient(region.haploinsufficiency)
    if sv == "DUP":
        return dosage_is_sufficient(region.triplosensitivity)
    return False


def region_coverage_fraction(start: int, end: int, region: DosageRegion) -> float:
    lo, hi = min(start, end), max(start, end)
    overlap = max(0, min(hi, region.end) - max(lo, region.start) + 1)
    denom = max(1, region.end - region.start + 1)
    return overlap / denom


def build_reference_index(
    open_root: Path,
    research_root: Path,
    index_dir: Path,
    force: bool = False,
) -> Path:
    open_manifest_path = open_root / "open_reference_manifest.json"
    research_manifest_path = research_root / "research_reference_manifest.json"
    open_manifest = load_json(open_manifest_path)
    research_manifest = load_json(research_manifest_path)
    if open_manifest.get("schema") != "open-reference-bundle-v2" or open_manifest.get("status") != "complete":
        raise RuntimeError("Open reference bundle is not READY. Run item 22 first.")
    if research_manifest.get("schema") != "research-reference-bundle-v2" or research_manifest.get("status") != "complete":
        raise RuntimeError("Research reference bundle is not READY. Run item 23 first.")

    index_dir.mkdir(parents=True, exist_ok=True)
    db = index_dir / "knowledge_index.sqlite3"
    state_path = index_dir / "knowledge_index.state.json"
    signature = {
        "schema": INDEX_SCHEMA,
        "open_manifest_sha256": sha256_file(open_manifest_path),
        "research_manifest_sha256": sha256_file(research_manifest_path),
        "rules_schema": RULES_SCHEMA,
    }
    state = load_json(state_path)
    if not force and db.is_file() and state.get("status") == "complete" and all(
        state.get(k) == v for k, v in signature.items()
    ):
        try:
            con = sqlite3.connect(db)
            ok = con.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
            con.close()
            if ok:
                print(f"[CURRENT] Knowledge reference index: {db}")
                return db
        except Exception:
            pass

    tmp = db.with_name(db.name + ".part")
    tmp.unlink(missing_ok=True)
    con = sqlite3.connect(tmp)
    counts: Counter[str] = Counter()
    try:
        con.executescript("""
            PRAGMA journal_mode=OFF;
            PRAGMA synchronous=OFF;
            CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE genes(symbol TEXT PRIMARY KEY, hgnc_id TEXT, name TEXT);
            CREATE TABLE gene_alias(alias TEXT PRIMARY KEY, symbol TEXT NOT NULL, hgnc_id TEXT);
            CREATE TABLE ontology_terms(term_id TEXT PRIMARY KEY, name TEXT, ontology TEXT);
            CREATE TABLE gene_disease(
                id INTEGER PRIMARY KEY, gene TEXT NOT NULL, hgnc_id TEXT, disease_id TEXT,
                disease_name TEXT, classification TEXT, inheritance TEXT, source TEXT NOT NULL,
                UNIQUE(gene,disease_id,disease_name,classification,inheritance,source)
            );
            CREATE TABLE actionability(id INTEGER PRIMARY KEY, gene TEXT NOT NULL, context TEXT, source TEXT NOT NULL,
                UNIQUE(gene,context,source));
            CREATE TABLE dosage(id INTEGER PRIMARY KEY, gene TEXT NOT NULL, evidence TEXT, source TEXT NOT NULL,
                UNIQUE(gene,evidence,source));
            CREATE TABLE dosage_region(
                id INTEGER PRIMARY KEY, region_id TEXT NOT NULL, name TEXT, chrom TEXT NOT NULL,
                start INTEGER NOT NULL, end INTEGER NOT NULL, haploinsufficiency TEXT,
                triplosensitivity TEXT, source TEXT NOT NULL,
                UNIQUE(region_id,chrom,start,end,haploinsufficiency,triplosensitivity,source)
            );
            CREATE TABLE panelapp_gene(
                id INTEGER PRIMARY KEY, gene TEXT NOT NULL, panel_id TEXT, panel_name TEXT,
                confidence TEXT, inheritance TEXT, phenotype TEXT,
                UNIQUE(gene,panel_id,confidence,inheritance,phenotype)
            );
            CREATE TABLE source_snapshot(source TEXT PRIMARY KEY, manifest_sha256 TEXT, installed_utc TEXT, scope TEXT);
        """)

        # HGNC canonical symbols and aliases.
        hgnc = source_artifact(open_manifest, open_root, "HGNC", "complete")
        if hgnc is None:
            raise RuntimeError("HGNC complete set is missing from open reference bundle")
        for row in _read_dict_rows(hgnc, "\t"):
            symbol = first_value(row, ["symbol"])
            if not symbol:
                continue
            hgnc_id = first_value(row, ["hgnc_id", "hgncid"])
            name = first_value(row, ["name"])
            con.execute("INSERT OR REPLACE INTO genes(symbol,hgnc_id,name) VALUES(?,?,?)", (symbol, hgnc_id, name))
            aliases = {symbol}
            aliases.update(split_multi(first_value(row, ["alias_symbol", "aliassymbol"])))
            aliases.update(split_multi(first_value(row, ["prev_symbol", "prevsymbol"])))
            for alias in aliases:
                con.execute(
                    "INSERT OR IGNORE INTO gene_alias(alias,symbol,hgnc_id) VALUES(?,?,?)",
                    (alias.upper(), symbol, hgnc_id),
                )
            counts["genes"] += 1
        con.commit()

        alias_map = {str(a).upper(): (str(s), str(h or "")) for a, s, h in con.execute(
            "SELECT alias,symbol,hgnc_id FROM gene_alias"
        )}
        known_symbols = {symbol.upper() for symbol, _ in alias_map.values()}

        def canonical_gene(value: str) -> tuple[str, str]:
            for token in split_multi(value) or [value.strip()]:
                if token.upper() in alias_map:
                    return alias_map[token.upper()]
            return (value.strip(), "") if value.strip() else ("", "")

        # Mondo and HPO ontology names (HPO ontology itself is not phenotype-gene evidence).
        for source, key, ontology in (("Mondo", "ontology_obo", "MONDO"), ("HPO", "ontology_obo", "HPO")):
            p = source_artifact(open_manifest, open_root, source, key)
            if p:
                for term_id, name in parse_obo(p):
                    con.execute("INSERT OR REPLACE INTO ontology_terms(term_id,name,ontology) VALUES(?,?,?)", (term_id, name, ontology))
                    counts[f"ontology_{ontology}"] += 1
        con.commit()
        ontology_names = {str(i): str(n) for i, n in con.execute("SELECT term_id,name FROM ontology_terms")}

        # GenCC submissions.
        gencc = source_artifact(open_manifest, open_root, "GenCC", "submissions_new_tsv")
        if gencc:
            for row in _read_dict_rows(gencc, "\t"):
                gene_raw = first_value(row, ["gene_symbol", "genesymbol", "submitted_as_gene_symbol"], ["genesymbol"])
                gene, hgnc_id = canonical_gene(gene_raw)
                if not gene:
                    continue
                disease_id = first_value(row, ["disease_curie", "disease_id", "diseaseid"], ["diseasecurie"])
                disease_name = first_value(row, ["disease_title", "disease_name", "diseasename"], ["diseasetitle"])
                disease_name = ontology_names.get(disease_id, disease_name)
                classification = first_value(row, ["classification_title", "classification"], ["classification"])
                inheritance = first_value(row, ["moi_title", "mode_of_inheritance", "inheritance"], ["moi"])
                con.execute(
                    "INSERT OR IGNORE INTO gene_disease(gene,hgnc_id,disease_id,disease_name,classification,inheritance,source) VALUES(?,?,?,?,?,?,?)",
                    (gene, hgnc_id, disease_id, disease_name, classification, inheritance, "GenCC"),
                )
                counts["gencc_relations"] += 1

        # ClinGen gene-disease validity. Current downloads contain a preamble
        # before the true CSV header, so use the dedicated preamble-aware reader.
        validity = source_artifact(open_manifest, open_root, "ClinGen", "gene_disease_validity")
        if validity:
            for row in _read_clingen_rows(validity, ("gene", "classification", "disease")):
                gene_raw = first_value(row, ["gene symbol", "gene_symbol", "genesymbol"], ["genesymbol"])
                gene, hgnc_id = canonical_gene(gene_raw)
                if not gene:
                    continue
                disease_id = first_value(row, ["disease id", "disease_id", "mondo id", "mondoid"], ["diseaseid", "mondo"])
                disease_name = first_value(row, ["disease label", "disease name", "disease"], ["diseaselabel"])
                disease_name = ontology_names.get(disease_id, disease_name)
                classification = first_value(row, ["classification", "gene-disease validity classification"], ["classification"])
                inheritance = first_value(row, ["moi", "mode of inheritance", "inheritance"], ["inheritance", "moi"])
                con.execute(
                    "INSERT OR IGNORE INTO gene_disease(gene,hgnc_id,disease_id,disease_name,classification,inheritance,source) VALUES(?,?,?,?,?,?,?)",
                    (gene, hgnc_id, disease_id, disease_name, classification, inheritance, "ClinGen"),
                )
                counts["clingen_validity_relations"] += 1

        # ClinGen dosage sensitivity. Only sufficient HI/TS evidence may trigger K5.
        # ISCA region rows are retained with their native GRCh37 intervals so CNV/SV
        # events can be intersected even when no single HGNC gene explains the region.
        dosage = source_artifact(open_manifest, open_root, "ClinGen", "dosage_genes_regions")
        if dosage:
            for row in _read_clingen_rows(dosage, ("gene", "haploinsufficiency", "triplosensitivity")):
                entity = first_value(row, ["gene/region", "gene region", "gene symbol", "symbol"], ["generegion", "genesymbol"])
                entity_id = first_value(row, ["hgnc/isca", "hgnc isca", "hgnc_id"], ["hgncisca"])
                coord37 = first_value(row, ["grch37"], ["grch37"])
                hi = first_value(row, ["haploinsufficiency", "haploinsufficiency score"], ["haploinsufficiency"])
                ts = first_value(row, ["triplosensitivity", "triplosensitivity score"], ["triplosensitivity"])
                if not (dosage_is_sufficient(hi) or dosage_is_sufficient(ts)):
                    counts["clingen_dosage_non_sufficient_skipped"] += 1
                    continue
                evidence = f"HI={hi or '.'};TS={ts or '.'}"
                if entity_id.upper().startswith("ISCA-"):
                    interval = parse_grch37_interval(coord37)
                    if interval:
                        chrom, start37, end37 = interval
                        con.execute(
                            "INSERT OR IGNORE INTO dosage_region(region_id,name,chrom,start,end,haploinsufficiency,triplosensitivity,source) VALUES(?,?,?,?,?,?,?,?)",
                            (entity_id, entity, chrom, start37, end37, hi, ts, "ClinGen"),
                        )
                        counts["clingen_dosage_regions_sensitive"] += 1
                else:
                    # Older/test ClinGen dosage exports may omit the HGNC/ISCA
                    # identifier and GRCh37 columns. A canonical HGNC symbol is
                    # still sufficient to index a gene-level dosage assertion.
                    gene, _ = canonical_gene(entity)
                    if gene and gene.upper() in known_symbols:
                        con.execute("INSERT OR IGNORE INTO dosage(gene,evidence,source) VALUES(?,?,?)", (gene, evidence, "ClinGen"))
                        counts["clingen_dosage_genes_sensitive"] += 1

        # ClinGen actionability JSON: source schemas can evolve, so use key-aware recursive extraction.
        for key, label in (("actionability_adult", "adult"), ("actionability_pediatric", "pediatric")):
            p = source_artifact(open_manifest, open_root, "ClinGen", key)
            if not p:
                continue
            data = json.loads(p.read_text(encoding="utf-8"))

            def walk(obj: Any, path: tuple[str, ...] = (), context: str = "") -> None:
                gene_context = any("gene" in token for token in path)
                if isinstance(obj, dict):
                    local = " | ".join(
                        str(obj.get(k, "")) for k in ("condition", "disease", "name", "title") if obj.get(k)
                    )[:1000]
                    next_context = local or context
                    for k, v in obj.items():
                        walk(v, path + (str(k).lower(),), next_context)
                elif isinstance(obj, list):
                    for item in obj:
                        walk(item, path, context)
                elif isinstance(obj, str) and gene_context:
                    for token in re.split(r"[^A-Za-z0-9._-]+", obj):
                        gene, _ = canonical_gene(token)
                        if gene and token.upper() in alias_map:
                            con.execute(
                                "INSERT OR IGNORE INTO actionability(gene,context,source) VALUES(?,?,?)",
                                (gene, f"{label}:{context}".strip(":"), "ClinGen"),
                            )
                            counts["actionability"] += 1
            walk(data)

        # PanelApp research snapshot.
        panel_db = source_artifact(research_manifest, research_root, "PanelApp", "snapshot_sqlite")
        if panel_db:
            pc = sqlite3.connect(panel_db)
            try:
                for panel_id, panel_name, payload_json in pc.execute("SELECT panel_id,name,payload_json FROM panels"):
                    try:
                        payload = json.loads(payload_json)
                    except Exception:
                        continue
                    genes = payload.get("genes", []) if isinstance(payload, dict) else []
                    if not isinstance(genes, list):
                        continue
                    for entry in genes:
                        if not isinstance(entry, dict):
                            continue
                        gene_raw = str(entry.get("entity_name") or entry.get("gene_data", {}).get("gene_symbol") or entry.get("gene") or "")
                        gene, _ = canonical_gene(gene_raw)
                        if not gene or gene_raw.upper() not in alias_map:
                            continue
                        confidence = str(entry.get("confidence_level") or entry.get("confidence") or "")
                        inheritance = str(entry.get("mode_of_inheritance") or entry.get("inheritance") or "")
                        phen = entry.get("phenotypes") or entry.get("phenotype") or ""
                        if isinstance(phen, list):
                            phenotype = " | ".join(str(x) for x in phen)
                        else:
                            phenotype = str(phen)
                        con.execute(
                            "INSERT OR IGNORE INTO panelapp_gene(gene,panel_id,panel_name,confidence,inheritance,phenotype) VALUES(?,?,?,?,?,?)",
                            (gene, str(panel_id), str(panel_name or ""), confidence, inheritance, phenotype[:2000]),
                        )
                        counts["panelapp_gene_rows"] += 1
            finally:
                pc.close()

        con.execute("INSERT OR REPLACE INTO source_snapshot VALUES(?,?,?,?)", (
            "OPEN_BUNDLE", signature["open_manifest_sha256"], str(open_manifest.get("installed_utc", "")), "COMMERCIAL_ALLOWED"
        ))
        con.execute("INSERT OR REPLACE INTO source_snapshot VALUES(?,?,?,?)", (
            "RESEARCH_BUNDLE", signature["research_manifest_sha256"], str(research_manifest.get("installed_utc", "")), "RESEARCH_ONLY"
        ))
        con.execute("INSERT OR REPLACE INTO meta VALUES('schema',?)", (INDEX_SCHEMA,))
        con.execute("INSERT OR REPLACE INTO meta VALUES('rules_schema',?)", (RULES_SCHEMA,))
        con.execute("INSERT OR REPLACE INTO meta VALUES('pipeline_release',?)", (PIPELINE_RELEASE,))
        con.executescript("""
            CREATE INDEX idx_gene_disease_gene ON gene_disease(gene);
            CREATE INDEX idx_actionability_gene ON actionability(gene);
            CREATE INDEX idx_dosage_gene ON dosage(gene);
            CREATE INDEX idx_dosage_region_coord ON dosage_region(chrom,start,end);
            CREATE INDEX idx_panelapp_gene ON panelapp_gene(gene);
        """)
        con.commit()
        if con.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise RuntimeError("Knowledge reference SQLite integrity check failed")
    finally:
        con.close()

    os.replace(tmp, db)
    state_payload = {
        **signature,
        "status": "complete",
        "pipeline_release": PIPELINE_RELEASE,
        "created_utc": utcnow(),
        "database": file_identity(db),
        "counts": dict(counts),
        "note": "HPO ontology terms are indexed only as ontology vocabulary; no HPO gene/disease annotations are inferred.",
    }
    save_json_atomic(state_path, state_payload)
    print(f"[BUILT] Knowledge reference index: {db} | genes={counts['genes']:,}")
    return db


# ----------------------------- patient events -----------------------------

@dataclass
class PatientEvent:
    event_type: str
    event_key: str
    source_layers: set[str]
    chrom: str
    position: int
    end: int
    ref: str = "."
    alt: str = "."
    genotype: str = "."
    zygosity: str = "UNKNOWN"
    genes: set[str] = field(default_factory=set)
    primary_gene: str = ""
    clinvar_genes: set[str] = field(default_factory=set)
    consequence: str = "."
    impact: str = "."
    max_af: str = "."
    rsid: str = "."
    resolved_rsid: str = "."
    clinvar_significance: str = "."
    clinvar_review_status: str = "."
    clinvar_review_stars: int = 0
    clinvar_conditions: str = "."
    alpha_score: str = "."
    alpha_class: str = "."
    cadd_raw: str = "."
    cadd_phred: str = "."
    giab_flags: set[str] = field(default_factory=set)
    dosage_region_evidence: list[str] = field(default_factory=list)
    dosage_gene_evidence: dict[str, list[str]] = field(default_factory=dict)
    open_reasons: set[str] = field(default_factory=set)
    research_reasons: set[str] = field(default_factory=set)
    possible_biallelic_unphased: bool = False
    biallelic_genes: set[str] = field(default_factory=set)
    raw_context: dict[str, str] = field(default_factory=dict)


def canonical_row_gene(value: str, knowledge: GeneKnowledge) -> str:
    for token in re.split(r"[&,;|]", value or ""):
        token = token.strip()
        if not token or token in {".", "-"}:
            continue
        canonical, _ = knowledge.canonicalize(token)
        if canonical:
            return canonical
    return ""


def clinvar_gene_symbols(row: dict[str, str], knowledge: GeneKnowledge) -> set[str]:
    result: set[str] = set()
    for token in re.split(r"[|;,]", row.get("clinvar_current_geneinfo", "") or ""):
        symbol = token.split(":", 1)[0].strip()
        canonical, _ = knowledge.canonicalize(symbol)
        if canonical:
            result.add(canonical)
    return result


def row_symbols(row: dict[str, str], knowledge: GeneKnowledge) -> set[str]:
    result: set[str] = set()
    for field_name in ("symbol", "all_symbols"):
        for token in re.split(r"[&,;|]", row.get(field_name, "") or ""):
            token = token.strip()
            if not token or token in {".", "-"}:
                continue
            canonical, _ = knowledge.canonicalize(token)
            if canonical:
                result.add(canonical)
    # ClinVar GENEINFO is preserved separately for exact ClinVar scoping, but
    # remains part of the event-level gene union for provenance/search.
    result.update(clinvar_gene_symbols(row, knowledge))
    return result


def candidate_functional_class(event: PatientEvent) -> str:
    terms = split_consequences(event.consequence)
    if terms & LOF_TERMS:
        return "LOF"
    if terms & MISSENSE_TERMS:
        return "MISSENSE"
    if terms & FUNCTIONAL_TERMS:
        return "FUNCTIONAL"
    if event.impact.upper() == "HIGH":
        return "LOF"
    if event.impact.upper() == "MODERATE":
        return "FUNCTIONAL"
    return "OTHER"


def scan_small_catalog(path: Path, knowledge: GeneKnowledge, max_af: float, progress_every: int = 500_000) -> tuple[dict[str, PatientEvent], Counter[str]]:
    events: dict[str, PatientEvent] = {}
    counts: Counter[str] = Counter()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"Small-variant catalog has no TSV header: {path}")
        required = {"variant_key", "chrom", "position", "ref", "alt", "alt_carried", "vcf_filter"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise RuntimeError("Small-variant catalog missing required fields: " + ", ".join(sorted(missing)))
        for row in reader:
            counts["catalog_rows"] += 1
            if progress_every and counts["catalog_rows"] % progress_every == 0:
                print(f"[SCAN] small variants {counts['catalog_rows']:,} rows | candidates {len(events):,}", flush=True)
            if row.get("alt_carried", "").lower() != "yes" or row.get("vcf_filter", "") != "PASS":
                continue
            counts["carried_pass"] += 1
            try:
                pos = int(row.get("position", "0"))
            except ValueError:
                continue
            chrom = normalize_chrom(row.get("chrom", ""))
            if not canonical_chrom(chrom):
                continue
            ref, alt = row.get("ref", "."), row.get("alt", ".")
            primary_gene = canonical_row_gene(row.get("symbol", "") or row.get("all_symbols", ""), knowledge)
            cv_genes = clinvar_gene_symbols(row, knowledge) if row.get("clinvar_exact_match", "") == "yes" else set()
            event = PatientEvent(
                "SNV_INDEL", row.get("variant_key") or canonical_event_key(chrom, pos, ref, alt),
                {row.get("source_layer", "CATALOG")}, chrom, pos, pos + max(0, len(ref) - 1),
                ref=ref, alt=alt, genotype=row.get("GT", "."),
                zygosity=zygosity(row.get("GT", "."), ref, alt),
                genes=row_symbols(row, knowledge), primary_gene=primary_gene, clinvar_genes=cv_genes,
                consequence=row.get("consequence") or row.get("all_consequences", "."),
                impact=row.get("impact", "."), max_af=row.get("max_af", "."), rsid=row.get("rsid", "."),
                resolved_rsid=row.get("clinvar_current_rsid", ".") if row.get("clinvar_exact_match", "") == "yes" else ".",
                clinvar_significance=row.get("clinvar_current_significance", ".") if row.get("clinvar_exact_match", "") == "yes" else ".",
                clinvar_review_status=row.get("clinvar_current_review_status", ".") if row.get("clinvar_exact_match", "") == "yes" else ".",
                clinvar_review_stars=int(row.get("clinvar_review_stars", "0") or 0) if str(row.get("clinvar_review_stars", "0")).isdigit() else 0,
                clinvar_conditions=row.get("clinvar_current_conditions", ".") if row.get("clinvar_exact_match", "") == "yes" else ".",
                raw_context={
                    "hgvsc": row.get("hgvsc", "."), "hgvsp": row.get("hgvsp", "."),
                    "sift": row.get("sift", "."), "polyphen": row.get("polyphen", "."),
                    "DP": row.get("DP", "."), "GQ": row.get("GQ", "."), "AD": row.get("AD", "."),
                    "primary_gene": primary_gene or ".",
                    "clinvar_genes": ";".join(sorted(cv_genes)) or ".",
                    "all_consequences": row.get("all_consequences", ".") or ".",
                    "selected_consequence": row.get("consequence", ".") or ".",
                },
            )
            # Gene-specific rules must use the consequence paired with the selected
            # VEP gene/transcript, not the union of consequences from neighbouring
            # genes/transcripts.  Keep the union only in provenance.
            event.consequence = row.get("consequence", "") or row.get("all_consequences", "") or "."

            bucket = clinvar_bucket(event.clinvar_significance)
            observed_af = number(event.max_af)
            if bucket == "PATHOGENIC":
                # Very common one-star P/LP ClinVar records are often phenotype/trait
                # assertions (e.g. blood groups) rather than rare Mendelian disease
                # alleles. Preserve them, but do not let them outrank genuinely rare
                # or independently reviewed P/LP records.
                if observed_af is not None and observed_af > 0.05 and event.clinvar_review_stars < 2:
                    event.open_reasons.add("K1_EXACT_CLINVAR_PLP_COMMON")
                else:
                    event.open_reasons.add("K1_EXACT_CLINVAR_PLP")
            elif bucket == "MIXED":
                event.open_reasons.add("K1_EXACT_CLINVAR_MIXED")
            elif bucket == "CONFLICT":
                event.open_reasons.add("K1_EXACT_CLINVAR_CONFLICT")
            elif bucket == "ASSOCIATION":
                event.open_reasons.add("K1_EXACT_CLINVAR_ASSOCIATION")

            func = candidate_functional_class(event)
            rarity = rarity_status(event.max_af, max_af)
            rare = rarity == "RARE"
            af_unknown = rarity == "UNKNOWN"
            benign_exact = bucket == "BENIGN"
            # Gene-level functional rules are scoped to the selected VEP gene and
            # require either positively established rarity or an explicit AF_UNKNOWN
            # candidate label. An absent population frequency is not evidence of rarity.
            primary_profile = knowledge.profile(event.primary_gene) if event.primary_gene else None
            if primary_profile is not None and not benign_exact:
                if func == "LOF" and primary_profile.strong_relations:
                    if rare:
                        event.open_reasons.add("K2_RARE_LOF_VALIDATED_GENE")
                    elif af_unknown:
                        event.open_reasons.add("K2_LOF_VALIDATED_GENE_AF_UNKNOWN")
                if func == "MISSENSE" and primary_profile.strong_relations:
                    if rare:
                        event.open_reasons.add("K3_RARE_MISSENSE_VALIDATED_GENE")
                    elif af_unknown:
                        event.open_reasons.add("K3_MISSENSE_VALIDATED_GENE_AF_UNKNOWN")
                if func in {"LOF", "MISSENSE", "FUNCTIONAL"} and primary_profile.actionable:
                    if rare:
                        event.open_reasons.add("K4_CLINGEN_ACTIONABLE")
                    elif af_unknown:
                        event.open_reasons.add("K4_CLINGEN_ACTIONABLE_AF_UNKNOWN")
                # PanelApp is explicitly research-only; AF-unknown functional variants
                # remain useful research hypotheses, but common or ClinVar-benign alleles
                # are not promoted merely because the gene occurs on a panel.
                if (rare or af_unknown) and func in {"LOF", "MISSENSE", "FUNCTIONAL"} and primary_profile.panelapp:
                    event.research_reasons.add("K7_PANELAPP_RESEARCH")
            # Predictor targets are research-only: AlphaMissense and CADD are
            # research-only. Include AF-unknown missense candidates, but not exact
            # ClinVar-benign alleles.
            if func == "MISSENSE" and (rare or af_unknown) and not benign_exact:
                if event.open_reasons or (primary_profile is not None and primary_profile.strong_relations):
                    event.raw_context["alpha_missense_target"] = "yes"
                if event.open_reasons or event.research_reasons:
                    event.research_reasons.add("K8_CADD_SUPPORT_TARGET")

            if event.open_reasons or event.research_reasons:
                events[event.event_key] = event
                counts["candidate_rows"] += 1
                for reason in event.open_reasons:
                    counts[f"open:{reason}"] += 1
                for reason in event.research_reasons:
                    counts[f"research:{reason}"] += 1
    return events, counts


def parse_reference_targets(path: Path, wanted_rsids: set[str]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    if not path.is_file() or not wanted_rsids:
        return result
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rsid = row.get("rsid", "")
            if rsid in wanted_rsids:
                result[rsid] = {str(k): str(v or "") for k, v in row.items()}
    return result


def add_stage17_fastq_events(
    events: dict[str, PatientEvent], evidence_path: Path | None, targets_path: Path | None,
    knowledge: GeneKnowledge, clinvar_path: Path | None,
) -> Counter[str]:
    counts: Counter[str] = Counter()
    if evidence_path is None or targets_path is None or not evidence_path.is_file() or not targets_path.is_file():
        return counts
    evidence_rows: list[dict[str, str]] = []
    wanted: set[str] = set()
    with open_text(evidence_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("tier") not in {"A", "B"} or row.get("selected_in_export") != "yes":
                continue
            genotype = row.get("genotype", "").replace("/", "").replace("|", "")
            if not genotype or genotype in {"--", "."}:
                continue
            evidence_rows.append({str(k): str(v or "") for k, v in row.items()})
            wanted.add(row.get("rsid", ""))
    targets = parse_reference_targets(targets_path, wanted)
    candidate_keys: set[str] = set()
    stage17_events: dict[str, PatientEvent] = {}
    for row in evidence_rows:
        target = targets.get(row.get("rsid", ""))
        if not target:
            continue
        ref = target.get("ref", ".").upper()
        genotype = row.get("genotype", "").upper().replace("/", "").replace("|", "")
        try:
            pos = int(row.get("position", "0"))
        except ValueError:
            continue
        chrom = normalize_chrom(row.get("chromosome", ""))
        if ref in {"", "."} or not canonical_chrom(chrom):
            continue
        for alt in split_multi(target.get("alts", "")):
            alt = alt.upper()
            if not alt or alt == ref or alt not in genotype:
                continue
            key = canonical_event_key(chrom, pos, ref, alt)
            if key in events:
                events[key].source_layers.add("FASTQ_STAGE17")
                continue
            ev = PatientEvent(
                "SNV_INDEL", key, {"FASTQ_STAGE17"}, chrom, pos, pos + max(0, len(ref) - 1),
                ref=ref, alt=alt, genotype=row.get("genotype", "."), zygosity=zygosity(row.get("genotype", "."), ref, alt),
                rsid=row.get("rsid", "."), resolved_rsid=row.get("rsid", "."),
                raw_context={"fastq_tier": row.get("tier", "."), "fastq_dp": row.get("dp", "."), "fastq_quality": row.get("quality", ".")},
            )
            stage17_events[key] = ev
            candidate_keys.add(key)
            counts["stage17_alt_events"] += 1

    # Exact current ClinVar lookup for the small set of FASTQ-discovered ALT alleles.
    if clinvar_path and clinvar_path.is_file() and candidate_keys:
        with open_text(clinvar_path) as handle:
            for raw in handle:
                if not raw or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 8:
                    continue
                chrom, pos_s, ident, ref, alts = fields[:5]
                try:
                    pos = int(pos_s)
                except ValueError:
                    continue
                info = {}
                for token in fields[7].split(";"):
                    k, sep, v = token.partition("=")
                    info[k] = v if sep else "1"
                for alt in alts.split(","):
                    key = canonical_event_key(chrom, pos, ref, alt)
                    ev = stage17_events.get(key)
                    if ev is None:
                        continue
                    ev.clinvar_significance = info.get("CLNSIG", ".")
                    ev.clinvar_review_status = info.get("CLNREVSTAT", ".")
                    ev.clinvar_review_stars = review_stars(ev.clinvar_review_status)
                    ev.clinvar_conditions = info.get("CLNDN", ".")
                    if info.get("RS"):
                        ev.resolved_rsid = "rs" + re.sub(r"^rs", "", info["RS"], flags=re.I)
                    genes = info.get("GENEINFO", "")
                    cv_genes: set[str] = set()
                    for token in genes.split("|"):
                        sym = token.split(":", 1)[0]
                        canonical, _ = knowledge.canonicalize(sym)
                        if canonical:
                            ev.genes.add(canonical)
                            cv_genes.add(canonical)
                    ev.clinvar_genes.update(cv_genes)
                    if not ev.primary_gene and len(cv_genes) == 1:
                        ev.primary_gene = next(iter(cv_genes))
                    ev.raw_context["clinvar_genes"] = ";".join(sorted(cv_genes)) or "."
                    bucket = clinvar_bucket(ev.clinvar_significance)
                    if bucket == "PATHOGENIC":
                        # Stage-17 FASTQ evidence has no population-frequency column,
                        # so an exact pure P/LP ClinVar assertion is retained without
                        # attempting the common-allele demotion used by the catalog.
                        ev.open_reasons.add("K1_EXACT_CLINVAR_PLP")
                    elif bucket == "MIXED":
                        ev.open_reasons.add("K1_EXACT_CLINVAR_MIXED")
                    elif bucket == "CONFLICT":
                        ev.open_reasons.add("K1_EXACT_CLINVAR_CONFLICT")
                    elif bucket == "ASSOCIATION":
                        ev.open_reasons.add("K1_EXACT_CLINVAR_ASSOCIATION")
        for key, ev in stage17_events.items():
            if ev.open_reasons:
                events[key] = ev
                counts["stage17_promoted_by_clinvar"] += 1
            else:
                counts["stage17_not_promoted"] += 1
    return counts


def annotsv_gene_fields(row: dict[str, str]) -> list[str]:
    candidates = []
    for key, value in row.items():
        nk = norm_header(key)
        if any(x in nk for x in ("genename", "genesymbol", "geneoverlap", "gene1name", "gene2name")):
            candidates.extend(split_multi(value))
    return candidates


def parse_info_kv(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for token in (value or "").split(";"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, val = token.split("=", 1)
            result[key] = val
        else:
            result[token] = "1"
    return result


def parse_format_sample(format_value: str, sample_value: str) -> dict[str, str]:
    keys = (format_value or "").split(":")
    vals = (sample_value or "").split(":")
    return {k: vals[i] if i < len(vals) else "" for i, k in enumerate(keys) if k}


def second_int(value: str) -> int:
    try:
        parts = (value or "").split(",")
        return int(float(parts[1])) if len(parts) > 1 and parts[1] not in {"", "."} else 0
    except Exception:
        return 0


def structural_cluster_key(row: dict[str, str], kind: str, chrom: str, start: int, end: int, sv_type: str) -> str:
    caller_id = first_value(row, ["ID", "id"], ["id"])
    info = parse_info_kv(first_value(row, ["INFO", "info"], ["info"]))
    mate_id = info.get("MATEID", "")
    if caller_id.startswith("MantaBND") and mate_id.startswith("MantaBND"):
        return "MANTA_PAIR:" + "|".join(sorted((caller_id, mate_id)))
    if caller_id:
        return f"{kind.upper()}_CALL:{caller_id}"
    annotsv_id = first_value(row, ["AnnotSV_ID"], ["annotsvid"])
    if annotsv_id:
        return f"{kind.upper()}_ANNOTSV:{annotsv_id}"
    return f"{kind.upper()}:{chrom}:{start}:{end}:{sv_type}"


def structural_member_context(row: dict[str, str]) -> dict[str, Any]:
    info_text = first_value(row, ["INFO", "info"], ["info"])
    info = parse_info_kv(info_text)
    sample = parse_format_sample(
        first_value(row, ["FORMAT", "format"], ["format"]),
        first_value(row, ["GFX0242962"], []),
    )
    gq = 0
    try:
        gq = int(float(sample.get("GQ", "0") or 0))
    except Exception:
        pass
    qual = 0.0
    try:
        qual = float(first_value(row, ["QUAL", "qual"], ["qual"]) or 0)
    except Exception:
        pass
    caller_id = first_value(row, ["ID", "id"], ["id"])
    caller = "MANTA" if caller_id.startswith("Manta") else ("DRAGEN" if caller_id.startswith("DRAGEN") else "OTHER")
    return {
        "caller_id": caller_id,
        "caller": caller,
        "filter": first_value(row, ["FILTER", "filter"], ["filter"]) or ".",
        "qual": qual,
        "gq": gq,
        "pr_alt": second_int(sample.get("PR", "")),
        "sr_alt": second_int(sample.get("SR", "")),
        "imprecise": "IMPRECISE" in info,
        "info": info_text,
        "sample": first_value(row, ["GFX0242962"], []) or ".",
        "annotsv_id": first_value(row, ["AnnotSV_ID"], ["annotsvid"]) or ".",
    }


def _merge_event_payload(dst: PatientEvent, src: PatientEvent) -> PatientEvent:
    dst.source_layers.update(src.source_layers)
    dst.genes.update(src.genes)
    dst.clinvar_genes.update(src.clinvar_genes)
    dst.open_reasons.update(src.open_reasons)
    dst.research_reasons.update(src.research_reasons)
    dst.giab_flags.update(src.giab_flags)
    dst.dosage_region_evidence = sorted(set(dst.dosage_region_evidence) | set(src.dosage_region_evidence))
    for gene, values in src.dosage_gene_evidence.items():
        dst.dosage_gene_evidence.setdefault(gene, [])
        dst.dosage_gene_evidence[gene] = sorted(set(dst.dosage_gene_evidence[gene]) | set(values))
    if src.position < dst.position:
        dst.position = src.position
    if src.end > dst.end:
        dst.end = src.end
    members = []
    for payload in (dst.raw_context.get("structural_members"), src.raw_context.get("structural_members")):
        if isinstance(payload, list):
            members.extend(payload)
    if members:
        dst.raw_context["structural_members"] = members
    return dst


def scan_structural(path: Path | None, kind: str, knowledge: GeneKnowledge) -> tuple[dict[str, PatientEvent], Counter[str]]:
    """Scan AnnotSV once per *biological call*, not once per split gene row.

    AnnotSV writes one `full` event row followed by multiple `split` rows.  Older
    Stage-24 rules accidentally materialized split rows as independent variants,
    multiplying structural findings.  Manta reciprocal BND mates can also describe
    the same event twice with slightly different homology-adjusted breakpoints.
    This scanner collapses those representations before any clinical knowledge is
    attached.
    """
    result: dict[str, PatientEvent] = {}
    counts: Counter[str] = Counter()
    if path is None or not path.is_file():
        return result, counts

    clusters: dict[str, dict[str, Any]] = {}
    with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return result, counts
        for row in reader:
            counts["rows"] += 1
            annotation_mode = (first_value(row, ["Annotation_mode"], ["annotationmode"]) or "").strip().lower()
            if annotation_mode and annotation_mode != "full":
                counts["split_rows_skipped"] += 1
                continue
            counts["full_rows"] += 1
            chrom = normalize_chrom(first_value(row, ["SV_chrom", "chrom", "chr"], ["svchrom"]))
            try:
                start = int(float(first_value(row, ["SV_start", "start"], ["svstart"]) or 0))
                end = int(float(first_value(row, ["SV_end", "end"], ["svend"]) or start))
            except ValueError:
                continue
            if not chrom or start <= 0 or not canonical_chrom(chrom):
                continue
            sv_type = (first_value(row, ["SV_type", "type"], ["svtype"]) or kind.upper()).upper()
            cluster_id = structural_cluster_key(row, kind, chrom, start, end, sv_type)
            genes: set[str] = set()
            for raw_gene in annotsv_gene_fields(row):
                canonical, _ = knowledge.canonicalize(raw_gene)
                if canonical:
                    genes.add(canonical)
            annotsv_class_raw = first_value(row, ["ACMG_class", "AnnotSV_class"], ["acmgclass", "annotsvclass"])
            annotsv_class = annotsv_class_number(annotsv_class_raw)
            member = structural_member_context(row)
            bucket = clusters.setdefault(cluster_id, {
                "kind": kind.upper(), "chrom": chrom, "start": start, "end": end, "sv_type": sv_type,
                "genes": set(), "annotsv_class": None, "annotsv_class_raw": ".", "members": [],
            })
            bucket["start"] = min(bucket["start"], start)
            bucket["end"] = max(bucket["end"], end)
            bucket["genes"].update(genes)
            bucket["members"].append(member)
            if annotsv_class is not None and (bucket["annotsv_class"] is None or annotsv_class > bucket["annotsv_class"]):
                bucket["annotsv_class"] = annotsv_class
                bucket["annotsv_class_raw"] = annotsv_class_raw

    counts["reciprocal_or_duplicate_full_rows_collapsed"] = sum(max(0, len(x["members"]) - 1) for x in clusters.values())

    for cluster_id, bucket in clusters.items():
        chrom = bucket["chrom"]
        start = int(bucket["start"])
        end = int(bucket["end"])
        sv_type = bucket["sv_type"]
        annotsv_class = bucket["annotsv_class"]
        annotsv_class_raw = bucket["annotsv_class_raw"]
        genes = set(bucket["genes"])
        key = f"GRCh37:{kind.upper()}:{chrom}:{start}:{end}:{sv_type}:{deterministic_id('SC', cluster_id)}"
        ev = PatientEvent(kind.upper(), key, {f"ANNOTSV_{kind.upper()}"}, chrom, start, end, genes=genes)
        members = bucket["members"]
        callers = sorted({str(x.get("caller", "")) for x in members if x.get("caller")})
        ev.raw_context = {
            "sv_type": sv_type,
            "structural_cluster_id": cluster_id,
            "annotsv_class": annotsv_class_raw,
            "structural_members": members,
            "callers": callers,
            "sv_length_abs": str(abs(end - start)),
        }

        # Biological evidence is attached here; technical call confidence is evaluated
        # later after Manta/DRAGEN cross-caller reconciliation.
        for gene in genes:
            profile = knowledge.profile(gene)
            if not profile:
                continue
            if sv_type in {"DEL", "DUP"} and (annotsv_class is not None and annotsv_class >= 3):
                hits = [e for e in profile.dosage if dosage_evidence_supports(e, sv_type)]
                if hits:
                    ev.dosage_gene_evidence[gene] = sorted(set(hits))
                    ev.open_reasons.add("K5_DOSAGE_SENSITIVE_CNV")
            if profile.strong_relations:
                if annotsv_class is not None and annotsv_class >= 4:
                    ev.open_reasons.add("K6_HIGH_EVIDENCE_SV_GENE")
                elif annotsv_class is None or annotsv_class == 3:
                    ev.open_reasons.add("K6_STRUCTURAL_GENE_CANDIDATE")
            if profile.panelapp and (annotsv_class is None or annotsv_class >= 3):
                ev.research_reasons.add("K7_PANELAPP_RESEARCH")

        if sv_type in {"DEL", "DUP"} and (annotsv_class is not None and annotsv_class >= 3):
            region_hits = [
                r for r in knowledge.overlapping_dosage_regions(chrom, start, end)
                if dosage_region_supports(r, sv_type) and region_coverage_fraction(start, end, r) >= 0.50
            ]
            if region_hits:
                ev.dosage_region_evidence = [r.evidence for r in region_hits]
                ev.open_reasons.add("K5_DOSAGE_SENSITIVE_CNV")

        result[key] = ev
        if ev.open_reasons or ev.research_reasons:
            counts["candidates"] += 1
    counts["normalized_events"] = len(result)
    return result, counts


def structural_length(event: PatientEvent) -> int:
    return max(0, int(event.end) - int(event.position))


def reciprocal_overlap(a: PatientEvent, b: PatientEvent) -> tuple[float, float]:
    if a.chrom != b.chrom:
        return 0.0, 0.0
    overlap = max(0, min(a.end, b.end) - max(a.position, b.position) + 1)
    la = max(1, a.end - a.position + 1)
    lb = max(1, b.end - b.position + 1)
    return overlap / la, overlap / lb


def _member_support(event: PatientEvent) -> dict[str, Any]:
    members = event.raw_context.get("structural_members", [])
    if not isinstance(members, list):
        members = []
    return {
        "callers": sorted({str(x.get("caller", "")) for x in members if isinstance(x, dict) and x.get("caller")}),
        "filter_pass": bool(members) and all(str(x.get("filter", "")).upper() == "PASS" for x in members if isinstance(x, dict)),
        "gq": max([int(x.get("gq", 0) or 0) for x in members if isinstance(x, dict)] or [0]),
        "qual": max([float(x.get("qual", 0) or 0) for x in members if isinstance(x, dict)] or [0.0]),
        "pr_alt": max([int(x.get("pr_alt", 0) or 0) for x in members if isinstance(x, dict)] or [0]),
        "sr_alt": max([int(x.get("sr_alt", 0) or 0) for x in members if isinstance(x, dict)] or [0]),
        "imprecise": any(bool(x.get("imprecise")) for x in members if isinstance(x, dict)),
    }


def single_caller_structural_confidence(event: PatientEvent) -> str:
    support = _member_support(event)
    sv_type = str(event.raw_context.get("sv_type", "")).upper()
    length = structural_length(event)
    callers = set(support["callers"])

    if "DRAGEN" in callers:
        # DRAGEN CNV PASS calls already incorporate normalized depth/copy-number
        # evidence. They are still caller-supported rather than independently confirmed.
        return "CALLER_SUPPORTED_DRAGEN_CNV"

    if "MANTA" in callers:
        if not support["filter_pass"]:
            return "CAUTION_FILTER"
        if support["imprecise"]:
            return "CAUTION_IMPRECISE"
        strong_reads = (support["sr_alt"] >= 3) or (support["pr_alt"] >= 5)
        if support["gq"] < 20 or not strong_reads:
            return "CAUTION_LOW_SUPPORT"
        if sv_type in {"DEL", "DUP"} and length >= 1_000_000:
            # A multi-megabase copy-number event in a viable adult is a high-impact
            # claim. One caller, even with read support, is not technical confirmation.
            return "CAUTION_LARGE_UNCONFIRMED"
        return "SINGLE_CALLER_STRONG"

    return "CAUTION_UNKNOWN_CALLER"


def collapse_near_duplicate_structural(events: dict[str, PatientEvent]) -> tuple[dict[str, PatientEvent], int]:
    """Collapse near-identical long DEL/DUP/INV representations from the same caller.

    This is intentionally conservative: same chromosome/type, both >=1 kb,
    starts and ends within 1 kb, and reciprocal overlap >=95%.
    TRA/BND are excluded because their AnnotSV interval representation is not a
    simple genomic interval and must not be merged by span geometry.
    """
    keys = list(events)
    by_group: dict[tuple[str, str], list[str]] = defaultdict(list)
    for key, ev in events.items():
        sv_type = str(ev.raw_context.get("sv_type", "")).upper()
        if sv_type in {"DEL", "DUP", "INV"} and structural_length(ev) >= 1000:
            by_group[(ev.chrom, sv_type)].append(key)

    removed: set[str] = set()
    collapsed = 0
    for _, group_keys in by_group.items():
        ordered = sorted(group_keys, key=lambda k: (events[k].position, events[k].end, k))
        for i, key in enumerate(ordered):
            if key in removed:
                continue
            base = events[key]
            for other_key in ordered[i + 1:]:
                if other_key in removed:
                    continue
                other = events[other_key]
                if other.position - base.position > 1000:
                    break
                if abs(other.end - base.end) > 1000:
                    continue
                ro1, ro2 = reciprocal_overlap(base, other)
                if min(ro1, ro2) < 0.95:
                    continue
                _merge_event_payload(base, other)
                members = base.raw_context.setdefault("near_duplicate_event_keys", [])
                if isinstance(members, list):
                    members.append(other_key)
                removed.add(other_key)
                collapsed += 1

    return {k: v for k, v in events.items() if k not in removed}, collapsed


def merge_cross_caller_structural(cnv_events: dict[str, PatientEvent], sv_events: dict[str, PatientEvent]) -> tuple[dict[str, PatientEvent], Counter[str]]:
    """Merge same-direction Manta + DRAGEN DEL/DUP calls with high reciprocal overlap."""
    counts: Counter[str] = Counter()
    cnv_events, cnv_near = collapse_near_duplicate_structural(cnv_events)
    sv_events, sv_near = collapse_near_duplicate_structural(sv_events)
    counts["cnv_near_duplicates_collapsed"] = cnv_near
    counts["sv_near_duplicates_collapsed"] = sv_near
    cnv_available = set(cnv_events)
    matches: list[tuple[float, str, str, float, float]] = []
    for skey, sev in sv_events.items():
        sv_type = str(sev.raw_context.get("sv_type", "")).upper()
        if sv_type not in {"DEL", "DUP"}:
            continue
        for ckey in list(cnv_available):
            cev = cnv_events[ckey]
            if str(cev.raw_context.get("sv_type", "")).upper() != sv_type or cev.chrom != sev.chrom:
                continue
            ro_s, ro_c = reciprocal_overlap(sev, cev)
            score = min(ro_s, ro_c)
            if score >= 0.70:
                matches.append((score, skey, ckey, ro_s, ro_c))

    used_sv: set[str] = set()
    used_cnv: set[str] = set()
    for score, skey, ckey, ro_s, ro_c in sorted(matches, reverse=True):
        if skey in used_sv or ckey in used_cnv:
            continue
        sev = sv_events[skey]
        cev = cnv_events[ckey]
        _merge_event_payload(sev, cev)
        sev.event_type = "STRUCTURAL"
        sev.raw_context["technical_confidence"] = "CROSS_CALLER_CONFIRMED"
        sev.raw_context["cross_caller_match"] = {
            "cnv_event_key": ckey,
            "reciprocal_overlap_sv": round(ro_s, 6),
            "reciprocal_overlap_cnv": round(ro_c, 6),
        }
        used_sv.add(skey)
        used_cnv.add(ckey)
        counts["cross_caller_merged"] += 1

    merged: dict[str, PatientEvent] = {}
    for key, ev in sv_events.items():
        if "technical_confidence" not in ev.raw_context:
            ev.raw_context["technical_confidence"] = single_caller_structural_confidence(ev)
        merged[key] = ev
    for key, ev in cnv_events.items():
        if key in used_cnv:
            continue
        ev.raw_context["technical_confidence"] = single_caller_structural_confidence(ev)
        merged[key] = ev

    high_allowed = {"CROSS_CALLER_CONFIRMED", "SINGLE_CALLER_STRONG", "CALLER_SUPPORTED_DRAGEN_CNV"}
    for ev in merged.values():
        tech = str(ev.raw_context.get("technical_confidence", "CAUTION_UNKNOWN"))
        if "K6_HIGH_EVIDENCE_SV_GENE" in ev.open_reasons and tech not in high_allowed:
            ev.open_reasons.discard("K6_HIGH_EVIDENCE_SV_GENE")
            ev.open_reasons.add("K6_STRUCTURAL_GENE_CANDIDATE")
            counts["k6_high_demoted_technical"] += 1
    counts["merged_structural_events"] = len(merged)
    return merged, counts


def event_technical_confidence(event: PatientEvent) -> str:
    base = str(event.raw_context.get("technical_confidence", "") or "")
    if not base:
        base = "STANDARD"
    if event.giab_flags:
        if base.startswith("CAUTION"):
            return base + "+GIAB_CAUTION"
        return base + "+GIAB_CAUTION"
    return base


# ----------------------------- indexed overlays -----------------------------

def docker_tabix(image: str, data: str, regions: Path, output: Path, mounts: list[tuple[Path, str]], include_header: bool = True, retries: int = 3) -> None:
    """Run tabix against a local or remote indexed resource.

    Remote BGZF/range queries occasionally fail transiently in libcurl/HTTP2.
    Retry the whole deterministic target slice; partial output is always discarded.
    Local files keep the same behavior but also benefit from one clean retry if Docker
    itself had a transient failure.
    """
    args = ["docker", "run", "--rm"]
    for host, container in mounts:
        args += ["--mount", f"type=bind,source={host},target={container}"]
    args += ["--entrypoint", "tabix", image]
    if include_header:
        args += ["-h"]
    args += ["-R", "/work/" + regions.name, data]
    attempts = max(1, int(retries))
    last_error = ""
    for attempt in range(1, attempts + 1):
        output.unlink(missing_ok=True)
        with output.open("wb") as handle:
            proc = subprocess.run(args, stdout=handle, stderr=subprocess.PIPE)
        if proc.returncode == 0:
            return
        last_error = proc.stderr.decode("utf-8", errors="replace")[-4000:]
        output.unlink(missing_ok=True)
        if attempt < attempts:
            time.sleep(min(2 ** attempt, 8))
    raise RuntimeError(f"tabix query failed after {attempts} attempt(s) ({proc.returncode}): {last_error}")


def write_regions(events: Iterable[PatientEvent], path: Path, only_snv: bool = False, chrom_prefix: str = "") -> int:
    rows = sorted({(e.chrom, e.position, e.position) for e in events if (not only_snv or len(e.ref) == len(e.alt) == 1)})
    with path.open("wt", encoding="ascii", newline="\n") as handle:
        for chrom, start, end in rows:
            handle.write(f"{chrom_prefix}{chrom}\t{start}\t{end}\n")
    return len(rows)


def apply_alphamissense(events: dict[str, PatientEvent], research_root: Path, research_manifest: dict[str, Any], image: str, work: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    targets = [e for e in events.values() if e.raw_context.get("alpha_missense_target") == "yes" and len(e.ref) == len(e.alt) == 1]
    if not targets:
        return counts
    data = source_artifact(research_manifest, research_root, "AlphaMissense", "hg19_predictions")
    if data is None:
        raise RuntimeError("AlphaMissense data missing from research bundle")
    region_file = work / "AlphaMissense_targets.regions.tsv"
    raw_output = work / "AlphaMissense_target_slice.tsv"
    counts["targets"] = len(targets)
    # Official AlphaMissense hg19 files use chr1..chrX contigs and 1-based POS.
    write_regions(targets, region_file, only_snv=True, chrom_prefix="chr")
    docker_tabix(image, "/alpha/" + data.name, region_file, raw_output, [(work, "/work"), (data.parent, "/alpha")])
    target_map = {(e.chrom, e.position, e.ref.upper(), e.alt.upper()): e for e in targets}
    matched_keys: set[str] = set()
    header: list[str] | None = None
    with raw_output.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("#"):
                candidate = raw.lstrip("#").rstrip("\r\n").split("\t")
                if len(candidate) >= 6:
                    header = candidate
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 6:
                continue
            counts["records_returned"] += 1
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            ev = target_map.get((normalize_chrom(fields[0]), pos, fields[2].upper(), fields[3].upper()))
            if ev is None:
                continue
            if header and len(header) == len(fields):
                row = dict(zip(header, fields))
                ev.alpha_score = first_value(row, ["am_pathogenicity", "pathogenicity"], ["pathogenicity"]) or "."
                ev.alpha_class = first_value(row, ["am_class", "class"], ["amclass"]) or "."
            else:
                # Official file: first 4 variant columns + transcript/protein columns + score/class at the end.
                ev.alpha_score = fields[-2] if len(fields) >= 2 else "."
                ev.alpha_class = fields[-1] if len(fields) >= 1 else "."
            matched_keys.add(ev.event_key)
    counts["matched"] = len(matched_keys)
    counts["not_in_canonical_database"] = max(0, counts["targets"] - counts["matched"])
    if counts["targets"] and counts["matched"] == 0:
        raise RuntimeError(
            "AlphaMissense target slice returned zero matches for "
            f"{counts['targets']} missense SNV target(s); contig/build/index mismatch suspected"
        )
    return counts


def cadd_remote_url(research_root: Path, research_manifest: dict[str, Any]) -> str:
    p = source_artifact(research_manifest, research_root, "CADD", "remote_score_source")
    if p:
        data = load_json(p)
        if data.get("score_file_url"):
            return str(data["score_file_url"])
    return "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz"


def _cadd_api_lookup(targets: list[PatientEvent], work: Path, attempts: int = 3) -> tuple[dict[str, tuple[str, str]], Counter[str]]:
    """Retrieve a modest Stage-24 SNV target set through the documented CADD API.

    CADD explicitly documents its API for small SNV batches and recommends tabix
    for larger jobs.  Stage 24 normally has only O(10^2) targets, so the API is a
    suitable fallback when htslib/libcurl cannot sustain remote BGZF HTTP range
    requests (notably HTTP/2 framing errors observed on Windows/Docker).
    """
    results: dict[str, tuple[str, str]] = {}
    counts: Counter[str] = Counter()
    out = work / "CADD_target_api.tsv"
    lines = ["chrom\tposition\tref\talt\traw\tphred\tstatus"]
    failures: list[str] = []
    consecutive_transport_failures = 0
    for ev in targets:
        endpoint = (
            "https://cadd.gs.washington.edu/api/v1.0/GRCh37-v1.7/"
            f"{ev.chrom}:{ev.position}_{ev.ref.upper()}_{ev.alt.upper()}"
        )
        payload = None
        last_error = ""
        for attempt in range(1, max(1, attempts) + 1):
            try:
                req = urllib.request.Request(endpoint, headers={"User-Agent": "Dante-WGS/9.0 Stage24"})
                with urllib.request.urlopen(req, timeout=30) as response:
                    payload = json.loads(response.read().decode("utf-8", errors="replace"))
                counts["api_requests_ok"] += 1
                break
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                if attempt < attempts:
                    time.sleep(min(2 ** (attempt - 1), 4))
        if payload is None:
            failures.append(f"{ev.event_key}: {last_error}")
            counts["api_requests_failed"] += 1
            consecutive_transport_failures += 1
            lines.append(f"{ev.chrom}\t{ev.position}\t{ev.ref}\t{ev.alt}\t.\t.\tERROR")
            # If the API host itself is unreachable, do not turn a 100-target
            # fallback into hours of repeated 30-second timeouts. Three
            # consecutive fully retried transport failures are enough to mark
            # the optional lookup deferred; no partial scores are applied.
            if consecutive_transport_failures >= 3:
                out.write_text("\n".join(lines) + "\n", encoding="utf-8")
                raise RuntimeError(
                    "CADD API fallback aborted after 3 consecutive transport failures; first errors: "
                    + " | ".join(failures[:3])
                )
            continue
        consecutive_transport_failures = 0
        record = None
        if isinstance(payload, list):
            for item in payload:
                if not isinstance(item, dict):
                    continue
                try:
                    pos = int(item.get("Pos", item.get("pos", 0)))
                except (TypeError, ValueError):
                    continue
                if (
                    normalize_chrom(str(item.get("Chrom", item.get("chrom", "")))) == ev.chrom
                    and pos == ev.position
                    and str(item.get("Ref", item.get("ref", ""))).upper() == ev.ref.upper()
                    and str(item.get("Alt", item.get("alt", ""))).upper() == ev.alt.upper()
                ):
                    record = item
                    break
        if record is None:
            counts["api_empty_or_mismatch"] += 1
            lines.append(f"{ev.chrom}\t{ev.position}\t{ev.ref}\t{ev.alt}\t.\t.\tNO_EXACT_MATCH")
            continue
        raw_score = str(record.get("RawScore", record.get("Raw", record.get("raw", "."))))
        phred = str(record.get("PHRED", record.get("Phred", record.get("phred", "."))))
        if raw_score in {"", "None"}:
            raw_score = "."
        if phred in {"", "None"}:
            phred = "."
        results[ev.event_key] = (raw_score, phred)
        counts["matched"] += 1
        lines.append(f"{ev.chrom}\t{ev.position}\t{ev.ref}\t{ev.alt}\t{raw_score}\t{phred}\tMATCH")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if failures:
        raise RuntimeError(
            f"CADD API fallback failed for {len(failures)}/{len(targets)} target(s); first errors: "
            + " | ".join(failures[:3])
        )
    return results, counts


def apply_cadd(events: dict[str, PatientEvent], research_root: Path, research_manifest: dict[str, Any], image: str, work: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    targets = [e for e in events.values() if "K8_CADD_SUPPORT_TARGET" in e.research_reasons and len(e.ref) == len(e.alt) == 1]
    if not targets:
        return counts
    counts["targets"] = len(targets)
    region_file = work / "CADD_targets.regions.tsv"
    raw_output = work / "CADD_target_slice.tsv"
    write_regions(targets, region_file, only_snv=True)
    url = cadd_remote_url(research_root, research_manifest)
    target_map = {(e.chrom, e.position, e.ref.upper(), e.alt.upper()): e for e in targets}
    matched_keys: set[str] = set()
    try:
        docker_tabix(image, url, region_file, raw_output, [(work, "/work")], include_header=True)
        with raw_output.open("rt", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                if raw.startswith("#") or not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 6:
                    continue
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                ev = target_map.get((normalize_chrom(fields[0]), pos, fields[2].upper(), fields[3].upper()))
                if ev is None:
                    continue
                ev.cadd_raw = fields[4]
                ev.cadd_phred = fields[5]
                matched_keys.add(ev.event_key)
        counts["matched"] = len(matched_keys)
        counts["transport_tabix"] = 1
        if counts["targets"] and counts["matched"] == 0:
            raise RuntimeError("CADD tabix target slice returned zero exact matches")
        return counts
    except RuntimeError as tabix_exc:
        print(f"[WARNING] CADD remote tabix failed; trying documented small-batch API fallback: {tabix_exc}", file=sys.stderr)
        # Do not expose partially parsed tabix results.  The fallback is applied
        # atomically only after every HTTP request has completed successfully.
        for ev in targets:
            ev.cadd_raw = "."
            ev.cadd_phred = "."
        api_results, api_counts = _cadd_api_lookup(targets, work)
        for ev in targets:
            values = api_results.get(ev.event_key)
            if values:
                ev.cadd_raw, ev.cadd_phred = values
        counts.update(api_counts)
        counts["matched"] = len(api_results)
        counts["transport_api_fallback"] = 1
        if counts["targets"] and counts["matched"] == 0:
            raise RuntimeError(f"CADD API fallback returned zero exact matches after tabix failure: {tabix_exc}")
        return counts


def parse_assembly_report(path: Path) -> dict[str, str]:
    """Map canonical GRCh37 chromosomes to their primary RefSeq accessions.

    NCBI assembly reports also contain many alt-scaffolds assigned to chromosome
    1..22/X/Y.  A naive last-row-wins map silently replaces NC_000001.10 etc.
    with NT_/NW_ alternate-locus accessions.  Remote dbSNP VCF regions must use
    the assembled-molecule RefSeq accession, so assembled-molecule rows are
    explicitly preferred and alt/patch scaffolds are never allowed to overwrite
    them.
    """
    result: dict[str, str] = {}
    priority: dict[str, int] = {}
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 7:
                continue
            sequence_name = fields[0]
            sequence_role = (fields[1] if len(fields) > 1 else "").strip().lower()
            assigned = fields[2] if len(fields) > 2 else ""
            refseq = fields[6] if len(fields) > 6 else ""
            chrom = assigned if assigned not in {"na", ""} else sequence_name
            if not refseq or refseq == "na" or not canonical_chrom(chrom):
                continue
            canonical = normalize_chrom(chrom)
            # Primary chromosome rows are the only safe default for coordinate
            # queries.  Keep a low-priority fallback for unusual assembly-report
            # variants, but never let an alt/patch scaffold replace a primary row.
            row_priority = 2 if sequence_role == "assembled-molecule" else 0
            if row_priority >= priority.get(canonical, -1):
                if row_priority > priority.get(canonical, -1) or canonical not in result:
                    result[canonical] = refseq
                    priority[canonical] = row_priority
    return result


def apply_dbsnp(events: dict[str, PatientEvent], open_root: Path, open_manifest: dict[str, Any], image: str, work: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    targets = [e for e in events.values() if e.event_type == "SNV_INDEL" and (not e.resolved_rsid or e.resolved_rsid == ".")]
    if not targets:
        return counts
    counts["targets"] = len(targets)
    counts["known_rsid_targets"] = sum(1 for e in targets if e.rsid not in {"", "."})
    remote_p = source_artifact(open_manifest, open_root, "dbSNP", "remote_vcf_sources")
    assembly_p = source_artifact(open_manifest, open_root, "dbSNP", "grch37_assembly_report")
    if remote_p is None or assembly_p is None:
        return counts
    remote = load_json(remote_p)
    grch37 = {}
    if isinstance(remote, dict):
        # open_reference_manager currently keys the dual-build sources as GRCh37.p13 / GRCh38.p14.
        # Accept the shorter legacy key too so Stage 24 survives bundle-schema-compatible renames.
        candidate = remote.get("GRCh37")
        if isinstance(candidate, dict):
            grch37 = candidate
        else:
            for key, value in remote.items():
                if str(key).upper().startswith("GRCH37") and isinstance(value, dict):
                    grch37 = value
                    break
    url = grch37.get("url") or grch37.get("final_url")
    if not url:
        return counts
    contigs = parse_assembly_report(assembly_p)
    region_file = work / "dbSNP_targets.regions.tsv"
    target_map: dict[tuple[str, int, str, str], PatientEvent] = {}
    lines: list[tuple[str, int, int]] = []
    for e in targets:
        accession = contigs.get(e.chrom)
        if not accession:
            continue
        lines.append((accession, e.position, e.position))
        target_map[(accession, e.position, e.ref.upper(), e.alt.upper())] = e
        # Defensive alternate key in case a remote reader rewrites RefSeq contigs
        # back to assigned chromosome names in query output.
        target_map[(e.chrom, e.position, e.ref.upper(), e.alt.upper())] = e
    if not lines:
        return counts
    unique_lines = sorted(set(lines))
    counts["query_regions"] = len(unique_lines)
    with region_file.open("wt", encoding="ascii", newline="\n") as handle:
        for chrom, pos, end in unique_lines:
            handle.write(f"{chrom}\t{pos}\t{end}\n")
    out = work / "dbSNP_target_slice.tsv"
    args = [
        "docker", "run", "--rm", "--mount", f"type=bind,source={work},target=/work",
        "--entrypoint", "bcftools", image,
        "query", "-R", "/work/" + region_file.name,
        "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\n", str(url),
    ]
    with out.open("wb") as handle:
        proc = subprocess.run(args, stdout=handle, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        out.unlink(missing_ok=True)
        raise RuntimeError(f"dbSNP target query failed ({proc.returncode}): {proc.stderr.decode('utf-8', errors='replace')[-4000:]}")
    matched_keys: set[str] = set()
    queried_positions = {(chrom, pos) for chrom, pos, _ in unique_lines}
    with out.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 5:
                continue
            counts["records_returned"] += 1
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            out_chrom = fields[0]
            if (out_chrom, pos) in queried_positions:
                counts["records_at_requested_positions"] += 1
            for alt in fields[3].split(","):
                ev = target_map.get((out_chrom, pos, fields[2].upper(), alt.upper()))
                if ev is None:
                    ev = target_map.get((normalize_chrom(out_chrom), pos, fields[2].upper(), alt.upper()))
                if ev is None:
                    continue
                counts["exact_allele_records"] += 1
                rsids = re.findall(r"rs\d+", fields[4], flags=re.I)
                if rsids:
                    ev.resolved_rsid = ";".join(sorted(set(x.lower() for x in rsids), key=lambda x: int(x[2:])))
                    matched_keys.add(ev.event_key)
    counts["matched"] = len(matched_keys)
    if counts["known_rsid_targets"] and counts["matched"] == 0:
        raise RuntimeError(
            "dbSNP target slice returned zero exact matches despite "
            f"{counts['known_rsid_targets']} target(s) already carrying rsIDs; "
            f"regions={counts['query_regions']}, records={counts['records_returned']}, "
            f"requested_position_records={counts['records_at_requested_positions']}. "
            "Primary-contig/allele normalization mismatch suspected"
        )
    return counts


# ----------------------------- GIAB QC -----------------------------

def giab_beds(open_manifest: dict[str, Any], open_root: Path) -> list[tuple[str, Path]]:
    result = []
    source = open_manifest.get("sources", {}).get("GIAB", {})
    for key, artifact in source.get("artifacts", {}).items():
        rel = artifact.get("path") if isinstance(artifact, dict) else None
        if rel:
            p = open_root / str(rel)
            if p.is_file() and p.name.endswith(".bed.gz"):
                result.append((str(key).upper(), p))
    return result


def apply_giab(events: dict[str, PatientEvent], beds: list[tuple[str, Path]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    positions: dict[str, list[int]] = defaultdict(list)
    pos_events: dict[tuple[str, int], list[PatientEvent]] = defaultdict(list)
    for ev in events.values():
        # Small variants have one relevant coordinate. For CNV/SV use both breakpoints;
        # a long event must not be labelled difficult merely because its span crosses a GIAB BED.
        event_positions = {ev.position}
        if ev.event_type in {"CNV", "SV"} and ev.end > 0:
            event_positions.add(ev.end)
        for pos in event_positions:
            positions[ev.chrom].append(pos)
            pos_events[(ev.chrom, pos)].append(ev)
    for chrom in positions:
        positions[chrom] = sorted(set(positions[chrom]))
    for label, path in beds:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                if raw.startswith("#") or not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom = normalize_chrom(fields[0])
                vals = positions.get(chrom)
                if not vals:
                    continue
                try:
                    # BED is 0-based half-open; patient positions are 1-based.
                    start1 = int(fields[1]) + 1
                    end1 = int(fields[2])
                except ValueError:
                    continue
                lo = bisect.bisect_left(vals, start1)
                hi = bisect.bisect_right(vals, end1)
                for pos in vals[lo:hi]:
                    for ev in pos_events[(chrom, pos)]:
                        ev.giab_flags.add(label)
                        counts[label] += 1
    return counts


# ----------------------------- findings -----------------------------

def _relation_disease_key(rel: GeneDisease) -> str:
    disease_id = (rel.disease_id or "").strip()
    if disease_id and disease_id not in {".", "-"}:
        return disease_id.upper()
    disease_name = re.sub(r"\s+", " ", (rel.disease_name or "").strip()).casefold()
    return f"NAME:{disease_name}" if disease_name else "UNKNOWN_DISEASE"


def relation_rows(event: PatientEvent, gene: str, profile: GeneProfile | None) -> list[GeneDisease | None]:
    """Collapse duplicate GenCC/ClinGen assertions into semantic gene-disease rows.

    A finding is a patient event + disease relationship, not a source assertion.
    Assertions for the same disease are first grouped by disease identity and then
    by normalized mode of inheritance.  If one source supplies a specific mode
    (for example X-linked) while another source leaves the mode unknown, the
    unknown assertion is absorbed into the specific bucket and contributes only
    provenance.  This avoids duplicate patient findings while preserving genuine
    AD-vs-AR disagreements as distinct semantic rows.
    """
    if not profile or not profile.strong_relations:
        return [None]

    by_disease: dict[str, list[GeneDisease]] = defaultdict(list)
    for rel in profile.strong_relations:
        by_disease[_relation_disease_key(rel)].append(rel)

    collapsed: list[GeneDisease] = []
    for disease_key, disease_rels in by_disease.items():
        by_inheritance: dict[str, list[GeneDisease]] = defaultdict(list)
        for rel in disease_rels:
            by_inheritance[inheritance_bucket(rel.inheritance)].append(rel)

        specific = [bucket for bucket in by_inheritance if bucket != "UNKNOWN"]
        # A missing MOI is not an independent biological assertion when every
        # source that *does* specify MOI agrees on one bucket.  Merge its source
        # provenance into that bucket instead of emitting a second finding.
        if "UNKNOWN" in by_inheritance and len(specific) == 1:
            by_inheritance[specific[0]].extend(by_inheritance.pop("UNKNOWN"))

        for inh_bucket, rels in by_inheritance.items():
            best = max(rels, key=lambda r: (classification_rank(r.classification), r.source, r.disease_name))
            disease_id = next((r.disease_id for r in rels if r.disease_id and r.disease_id not in {".", "-"}), best.disease_id)
            disease_name = next((r.disease_name for r in rels if r.disease_name and r.disease_name not in {".", "-"}), best.disease_name)
            hgnc_id = next((r.hgnc_id for r in rels if r.hgnc_id), best.hgnc_id)
            sources = ";".join(sorted({r.source for r in rels if r.source}))
            collapsed.append(GeneDisease(
                gene=best.gene, hgnc_id=hgnc_id, disease_id=disease_id, disease_name=disease_name,
                classification=best.classification, inheritance=inh_bucket, source=sources,
            ))

    collapsed.sort(
        key=lambda r: (_relation_disease_key(r), inheritance_bucket(r.inheritance),
                       -classification_rank(r.classification), r.source)
    )
    return collapsed


def evidence_class_for_reason(reason: str, event: PatientEvent) -> tuple[int, str]:
    if reason == "K1_EXACT_CLINVAR_PLP":
        if event.clinvar_review_stars >= 2:
            return 1, "OPEN_CLINICAL_HIGH"
        return 2, "OPEN_CLINICAL"
    if reason == "K1_EXACT_CLINVAR_PLP_COMMON":
        return 4, "OPEN_CLINICAL_COMMON_PLP_REVIEW"
    if reason == "K1_EXACT_CLINVAR_MIXED":
        return 4, "OPEN_CLINICAL_MIXED"
    if reason == "K1_EXACT_CLINVAR_CONFLICT":
        return 4, "OPEN_CLINICAL_CONFLICT"
    if reason == "K5_DOSAGE_SENSITIVE_CNV":
        acmg = annotsv_class_number(event.raw_context.get("annotsv_class", ""))
        if acmg is not None and acmg >= 4:
            return 2, "OPEN_DOSAGE_HIGH"
        return 4, "OPEN_DOSAGE_CANDIDATE"
    if reason in {"K2_RARE_LOF_VALIDATED_GENE", "K6_HIGH_EVIDENCE_SV_GENE"}:
        return 3, "OPEN_GENE_DISEASE_SUPPORT"
    if reason in {
        "K2_LOF_VALIDATED_GENE_AF_UNKNOWN",
        "K3_MISSENSE_VALIDATED_GENE_AF_UNKNOWN",
        "K4_CLINGEN_ACTIONABLE_AF_UNKNOWN",
    }:
        return 5, "OPEN_GENE_CANDIDATE_AF_UNKNOWN"
    if reason == "K6_STRUCTURAL_GENE_CANDIDATE":
        return 5, "OPEN_STRUCTURAL_CANDIDATE"
    if reason in {"K3_RARE_MISSENSE_VALIDATED_GENE", "K4_CLINGEN_ACTIONABLE"}:
        return 4, "OPEN_HYPOTHESIS_SUPPORTED"
    if reason.startswith("K8_"):
        return 5, "OPEN_PREDICTOR_SUPPORT"
    return 5, "OPEN_CONTEXT"


def inheritance_context(event: PatientEvent, relation: GeneDisease | None) -> str:
    if relation is None:
        return "UNKNOWN"
    bucket = inheritance_bucket(relation.inheritance)
    if bucket == "AR" and event.zygosity in {"HET_ALT", "HEMIZYGOUS_ALT"}:
        return "AR_SINGLE_HET_CARRIER_OR_CANDIDATE"
    if bucket == "AR" and event.zygosity == "HOM_ALT":
        return "AR_HOMOZYGOUS_ALT"
    return bucket


def patient_priority(reason: str, event: PatientEvent, gene: str, profile: GeneProfile | None, inheritance: str, evidence_rank: int) -> tuple[int, str]:
    """Keep evidence strength separate from relevance to this person's genotype/call confidence."""
    rank = evidence_rank
    flags: list[str] = []

    if inheritance == "AR_SINGLE_HET_CARRIER_OR_CANDIDATE":
        rank = max(rank, 4)
        flags.append("AR_SINGLE_HET_CARRIER_OR_CANDIDATE")
    elif reason.startswith("K1_EXACT_CLINVAR_PLP") and gene and profile and event.zygosity == "HET_ALT":
        modes = {
            inheritance_bucket(r.inheritance)
            for r in profile.strong_relations
            if inheritance_bucket(r.inheritance) != "UNKNOWN"
        }
        if modes == {"AR"}:
            rank = max(rank, 4)
            flags.append("LIKELY_AR_CARRIER_CONTEXT")

    tech = event_technical_confidence(event)
    if tech.startswith("CAUTION"):
        if reason.startswith(("K5_", "K6_")):
            rank = max(rank, 5)
        flags.append("TECHNICAL_CALL_CAUTION")
    if "GIAB_CAUTION" in tech:
        flags.append("GIAB_DIFFICULT_REGION")

    return rank, (";".join(flags) if flags else "DIRECT_EVIDENCE_CONTEXT")


def mark_possible_biallelic(events: dict[str, PatientEvent], knowledge: GeneKnowledge, max_af: float) -> int:
    by_gene: dict[str, list[PatientEvent]] = defaultdict(list)
    for ev in events.values():
        if ev.event_type != "SNV_INDEL" or ev.zygosity != "HET_ALT":
            continue
        if not is_rare(ev.max_af, max_af) or clinvar_bucket(ev.clinvar_significance) == "BENIGN":
            continue
        if candidate_functional_class(ev) not in {"LOF", "MISSENSE", "FUNCTIONAL"}:
            continue
        gene = ev.primary_gene
        if not gene:
            continue
        p = knowledge.profile(gene)
        if p and any(inheritance_bucket(r.inheritance) == "AR" for r in p.strong_relations):
            by_gene[gene].append(ev)
    marked = 0
    for gene, variants in by_gene.items():
        unique = {v.event_key: v for v in variants}
        if len(unique) >= 2:
            for ev in unique.values():
                ev.possible_biallelic_unphased = True
                ev.biallelic_genes.add(gene)
                marked += 1
    return marked


def open_findings(events: dict[str, PatientEvent], knowledge: GeneKnowledge) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    gene_validity_reasons = {
        "K2_RARE_LOF_VALIDATED_GENE",
        "K2_LOF_VALIDATED_GENE_AF_UNKNOWN",
        "K3_RARE_MISSENSE_VALIDATED_GENE",
        "K3_MISSENSE_VALIDATED_GENE_AF_UNKNOWN",
        "K6_HIGH_EVIDENCE_SV_GENE",
        "K6_STRUCTURAL_GENE_CANDIDATE",
    }

    for ev in events.values():
        for reason in sorted(ev.open_reasons):
            # Scope each evidence type to the entity that actually generated it.
            # This prevents event-level evidence from leaking to every nearby VEP symbol.
            entries: list[tuple[str, str]] = []  # (gene, scope)
            if reason.startswith("K1_"):
                entries = [(g, "clinvar") for g in sorted(ev.clinvar_genes)] or [("", "clinvar")]
            elif reason in {
                "K2_RARE_LOF_VALIDATED_GENE", "K2_LOF_VALIDATED_GENE_AF_UNKNOWN",
                "K3_RARE_MISSENSE_VALIDATED_GENE", "K3_MISSENSE_VALIDATED_GENE_AF_UNKNOWN",
                "K4_CLINGEN_ACTIONABLE", "K4_CLINGEN_ACTIONABLE_AF_UNKNOWN",
            }:
                if ev.primary_gene:
                    entries = [(ev.primary_gene, "primary")]
            elif reason in {"K6_HIGH_EVIDENCE_SV_GENE", "K6_STRUCTURAL_GENE_CANDIDATE"}:
                entries = [
                    (g, "structural_gene") for g in sorted(ev.genes)
                    if (knowledge.profile(g) and knowledge.profile(g).strong_relations)
                ]
            elif reason == "K5_DOSAGE_SENSITIVE_CNV":
                entries.extend((g, "dosage_gene") for g in sorted(ev.dosage_gene_evidence))
                if ev.dosage_region_evidence:
                    entries.append(("", "dosage_region"))
            else:
                entries = [(g, "default") for g in (sorted(ev.genes) or [""])]

            for gene, scope in entries:
                profile = knowledge.profile(gene) if gene else None
                if reason in {"K4_CLINGEN_ACTIONABLE", "K4_CLINGEN_ACTIONABLE_AF_UNKNOWN"} and not (profile and profile.actionable):
                    continue
                if reason in gene_validity_reasons and not (profile and profile.strong_relations):
                    continue

                relations = relation_rows(ev, gene, profile) if reason in gene_validity_reasons else [None]
                for rel in relations:
                    evidence_rank, evidence_class = evidence_class_for_reason(reason, ev)
                    actionability = "yes" if profile and profile.actionable else "no"
                    if scope == "dosage_gene":
                        dosage_items = ev.dosage_gene_evidence.get(gene, [])
                    elif scope == "dosage_region":
                        dosage_items = ev.dosage_region_evidence
                    else:
                        dosage_items = []
                    dosage = " | ".join(sorted(set(dosage_items))) if dosage_items else "."
                    hgnc_id = profile.hgnc_id if profile else ""
                    disease_id = rel.disease_id if rel else "."
                    if rel is not None:
                        disease_name = rel.disease_name
                    elif reason.startswith("K1_"):
                        disease_name = ev.clinvar_conditions
                    else:
                        disease_name = "."
                    classification = normalized_classification(rel.classification) if rel else "."
                    inheritance = inheritance_context(ev, rel)
                    giab = ";".join(sorted(ev.giab_flags)) or "."
                    tech = event_technical_confidence(ev)
                    rank, patient_relevance = patient_priority(reason, ev, gene, profile, inheritance, evidence_rank)
                    sources = {"VEP/CATALOG"}
                    if reason.startswith("K1_"):
                        sources.add("ClinVar")
                    if reason in gene_validity_reasons and rel is not None:
                        sources.update(x for x in rel.source.split(";") if x)
                    if reason in {"K4_CLINGEN_ACTIONABLE", "K4_CLINGEN_ACTIONABLE_AF_UNKNOWN"}:
                        sources.add("ClinGen Actionability")
                    if reason == "K5_DOSAGE_SENSITIVE_CNV":
                        sources.add("ClinGen Dosage")
                    if ev.giab_flags:
                        sources.add("GIAB")
                    if ev.resolved_rsid not in {"", "."} and ev.resolved_rsid != ev.rsid:
                        sources.add("dbSNP")
                    provenance = {
                        "source_layers": sorted(ev.source_layers),
                        "raw_context": ev.raw_context,
                        "evidence_scope": scope,
                        "primary_gene": ev.primary_gene or ".",
                        "clinvar_genes": sorted(ev.clinvar_genes),
                        "gene_profile": {
                            "strong_relation_sources": sorted({r.source for r in (profile.strong_relations if profile else [])}),
                            "actionability_context": profile.actionability_context if profile else [],
                        },
                        "dosage_region_evidence": ev.dosage_region_evidence,
                        "dosage_gene_evidence": ev.dosage_gene_evidence.get(gene, []) if gene else [],
                    }
                    disease_identity = disease_id if disease_id not in {"", ".", "-"} else f"NAME:{(disease_name or '.').strip().casefold()}"
                    finding_id = deterministic_id(
                        "OF", reason, ev.event_key, gene, disease_identity, classification, inheritance, scope
                    )
                    rows.append({
                        "finding_id": finding_id,
                        "priority_rank": str(rank), "evidence_rank": str(evidence_rank),
                        "patient_relevance": patient_relevance,
                        "reason_code": reason, "evidence_class": evidence_class,
                        "event_type": ev.event_type, "event_key": ev.event_key,
                        "source_layers": ";".join(sorted(ev.source_layers)), "chrom": ev.chrom,
                        "position": str(ev.position), "end": str(ev.end), "ref": ev.ref, "alt": ev.alt,
                        "genotype": ev.genotype, "zygosity": ev.zygosity, "gene": gene or ".",
                        "hgnc_id": hgnc_id or ".", "consequence": ev.consequence, "impact": ev.impact,
                        "max_af": ev.max_af, "rsid": ev.rsid, "resolved_rsid": ev.resolved_rsid,
                        "clinvar_significance": ev.clinvar_significance,
                        "clinvar_review_stars": str(ev.clinvar_review_stars), "clinvar_conditions": ev.clinvar_conditions,
                        "disease_id": disease_id, "disease_name": disease_name or ".",
                        "gene_disease_classification": classification, "inheritance": inheritance,
                        "actionability": actionability, "dosage_evidence": dosage,
                        "possible_biallelic_unphased": "yes" if gene and gene in ev.biallelic_genes else "no",
                        "giab_flags": giab, "technical_confidence": tech,
                        "evidence_sources": ";".join(sorted(sources)),
                        "interpretation_boundary": (
                            "Candidate/evidence finding, not a diagnosis. Predictor scores are supporting evidence only; "
                            "AR multiple heterozygotes remain unphased unless independently phased."
                        ),
                        "provenance_json": json.dumps(provenance, ensure_ascii=False, sort_keys=True),
                    })

    # The finding_id is the semantic primary key used in the patient SQLite.
    identity_fields = (
        "reason_code", "event_key", "gene", "disease_id", "gene_disease_classification", "inheritance"
    )
    merged: dict[str, dict[str, str]] = {}
    merged_counts: Counter[str] = Counter()
    for row in rows:
        fid = row["finding_id"]
        if fid not in merged:
            merged[fid] = row
            continue
        current = merged[fid]
        if any(current.get(k, "") != row.get(k, "") for k in identity_fields):
            details = [
                {k: x.get(k, "") for k in ("finding_id", *identity_fields, "disease_name")}
                for x in (current, row)
            ]
            raise RuntimeError(
                "Deterministic finding-id collision across different semantic identities: "
                + json.dumps(details, ensure_ascii=False)
            )
        merged_counts[fid] += 1
        for field_name in ("evidence_sources", "source_layers", "giab_flags"):
            tokens = set()
            for value in (current.get(field_name, ""), row.get(field_name, "")):
                tokens.update(t for t in value.split(";") if t and t != ".")
            current[field_name] = ";".join(sorted(tokens)) or "."
        if current.get("dosage_evidence") in {"", "."} and row.get("dosage_evidence") not in {"", "."}:
            current["dosage_evidence"] = row["dosage_evidence"]
        if row.get("actionability") == "yes":
            current["actionability"] = "yes"
        if row.get("technical_confidence", "").startswith("CAUTION"):
            current["technical_confidence"] = row["technical_confidence"]
        try:
            p1 = json.loads(current.get("provenance_json", "{}"))
            p2 = json.loads(row.get("provenance_json", "{}"))
            records = p1.get("merged_provenance") if isinstance(p1, dict) else None
            if not isinstance(records, list):
                records = [p1]
            records.append(p2)
            current["provenance_json"] = json.dumps(
                {"merged_provenance": records, "semantic_duplicate_rows_collapsed": len(records)},
                ensure_ascii=False, sort_keys=True
            )
        except Exception:
            pass

    rows = list(merged.values())
    if merged_counts:
        print(
            f"[DEDUP] collapsed {sum(merged_counts.values()):,} duplicate semantic open-finding row(s) "
            f"across {len(merged_counts):,} finding(s)",
            flush=True,
        )
    rows.sort(key=lambda r: (int(r["priority_rank"]), r["gene"], r["chrom"], int(r["position"]), r["reason_code"]))
    return rows


def research_overlays(events: dict[str, PatientEvent], open_rows: list[dict[str, str]], knowledge: GeneKnowledge) -> list[dict[str, str]]:
    base_by_event: dict[str, list[str]] = defaultdict(list)
    for row in open_rows:
        base_by_event[row["event_key"]].append(row["finding_id"])
    rows: list[dict[str, str]] = []
    for ev in events.values():
        if not ev.research_reasons:
            continue
        for reason in sorted(ev.research_reasons):
            if reason == "K8_CADD_SUPPORT_TARGET":
                if ev.cadd_phred in {"", "."}:
                    continue
                genes = [ev.primary_gene] if ev.primary_gene else [""]
            elif reason == "K8_ALPHAMISSENSE_SUPPORT":
                if ev.alpha_score in {"", "."}:
                    continue
                genes = [ev.primary_gene] if ev.primary_gene else [""]
            elif reason == "K7_PANELAPP_RESEARCH":
                if ev.event_type == "SNV_INDEL" and ev.primary_gene:
                    genes = [ev.primary_gene]
                else:
                    genes = sorted(ev.genes) or [""]
            else:
                genes = sorted(ev.genes) or [""]

            for gene in genes:
                profile = knowledge.profile(gene) if gene else None
                panels = profile.panelapp if profile else []
                if reason == "K7_PANELAPP_RESEARCH" and not panels:
                    continue
                panel_names = ";".join(sorted({p["panel_name"] for p in panels if p.get("panel_name")})) or "."
                panel_conf = ";".join(sorted({p["confidence"] for p in panels if p.get("confidence")})) or "."
                panel_inh = ";".join(sorted({p["inheritance"] for p in panels if p.get("inheritance")})) or "."
                rank = 5 if reason == "K7_PANELAPP_RESEARCH" else 6
                sources = []
                if reason == "K7_PANELAPP_RESEARCH":
                    sources.append("PanelApp")
                if reason == "K8_CADD_SUPPORT_TARGET":
                    sources.append("CADD")
                if reason == "K8_ALPHAMISSENSE_SUPPORT":
                    sources.append("AlphaMissense")
                if ev.giab_flags:
                    sources.append("GIAB")
                base_ids = ";".join(sorted(base_by_event.get(ev.event_key, []))) or "."
                giab = ";".join(sorted(ev.giab_flags)) or "."
                tech = event_technical_confidence(ev)
                overlay_id = deterministic_id("RO", reason, ev.event_key, gene)
                row_panel_names = panel_names if reason == "K7_PANELAPP_RESEARCH" else "."
                row_panel_conf = panel_conf if reason == "K7_PANELAPP_RESEARCH" else "."
                row_panel_inh = panel_inh if reason == "K7_PANELAPP_RESEARCH" else "."
                row_alpha_score = ev.alpha_score if reason == "K8_ALPHAMISSENSE_SUPPORT" else "."
                row_alpha_class = ev.alpha_class if reason == "K8_ALPHAMISSENSE_SUPPORT" else "."
                row_cadd_raw = ev.cadd_raw if reason == "K8_CADD_SUPPORT_TARGET" else "."
                row_cadd_phred = ev.cadd_phred if reason == "K8_CADD_SUPPORT_TARGET" else "."
                rows.append({
                    "overlay_id": overlay_id, "base_open_finding_id": base_ids,
                    "priority_rank": str(rank), "reason_code": reason,
                    "event_type": ev.event_type, "event_key": ev.event_key, "chrom": ev.chrom,
                    "position": str(ev.position), "end": str(ev.end), "ref": ev.ref, "alt": ev.alt,
                    "gene": gene or ".", "consequence": ev.consequence, "max_af": ev.max_af,
                    "panelapp_panels": row_panel_names, "panelapp_confidence": row_panel_conf,
                    "panelapp_inheritance": row_panel_inh,
                    "alpha_missense_pathogenicity": row_alpha_score, "alpha_missense_class": row_alpha_class,
                    "cadd_raw": row_cadd_raw, "cadd_phred": row_cadd_phred,
                    "giab_flags": giab, "technical_confidence": tech,
                    "evidence_sources": ";".join(sorted(set(sources))),
                    "research_banner": "RESEARCH ONLY — NOT FOR COMMERCIAL USE",
                    "interpretation_boundary": (
                        "Research/hypothesis overlay only. PanelApp/AlphaMissense/CADD do not independently establish pathogenicity, "
                        "diagnosis, treatment, or medical decision-making."
                    ),
                    "provenance_json": json.dumps({"source_layers": sorted(ev.source_layers), "primary_gene": ev.primary_gene or "."}, ensure_ascii=False, sort_keys=True),
                })
    rows.sort(key=lambda r: (int(r["priority_rank"]), r["gene"], r["chrom"], int(r["position"]), r["reason_code"]))
    return rows


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with atomic_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_patient_sqlite(path: Path, events: dict[str, PatientEvent], open_rows: list[dict[str, str]], research_rows: list[dict[str, str]], scan_stats: dict[str, Any]) -> None:
    part = path.with_name(path.name + ".part")
    part.unlink(missing_ok=True)
    con = sqlite3.connect(part)
    try:
        con.executescript("""
            PRAGMA journal_mode=OFF;
            PRAGMA synchronous=OFF;
            CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE event_candidates(
                event_key TEXT PRIMARY KEY, event_type TEXT, chrom TEXT, position INTEGER, end INTEGER,
                ref TEXT, alt TEXT, genotype TEXT, zygosity TEXT, genes TEXT, source_layers TEXT,
                consequence TEXT, impact TEXT, max_af TEXT, rsid TEXT, resolved_rsid TEXT,
                clinvar_significance TEXT, alpha_score TEXT, alpha_class TEXT, cadd_raw TEXT, cadd_phred TEXT,
                giab_flags TEXT, possible_biallelic_unphased INTEGER
            );
            CREATE TABLE open_findings(finding_id TEXT PRIMARY KEY, event_key TEXT, priority_rank INTEGER, reason_code TEXT, gene TEXT, payload_json TEXT NOT NULL);
            CREATE TABLE research_overlays(overlay_id TEXT PRIMARY KEY, event_key TEXT, priority_rank INTEGER, reason_code TEXT, gene TEXT, payload_json TEXT NOT NULL);
            CREATE TABLE scan_stats(key TEXT PRIMARY KEY, value TEXT NOT NULL);
        """)
        con.executemany(
            "INSERT INTO event_candidates VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [(
                e.event_key, e.event_type, e.chrom, e.position, e.end, e.ref, e.alt, e.genotype, e.zygosity,
                ";".join(sorted(e.genes)), ";".join(sorted(e.source_layers)), e.consequence, e.impact, e.max_af,
                e.rsid, e.resolved_rsid, e.clinvar_significance, e.alpha_score, e.alpha_class, e.cadd_raw, e.cadd_phred,
                ";".join(sorted(e.giab_flags)), 1 if e.possible_biallelic_unphased else 0,
            ) for e in events.values()]
        )
        con.executemany(
            "INSERT INTO open_findings VALUES(?,?,?,?,?,?)",
            [(r["finding_id"], r["event_key"], int(r["priority_rank"]), r["reason_code"], r["gene"], json.dumps(r, ensure_ascii=False, sort_keys=True)) for r in open_rows],
        )
        con.executemany(
            "INSERT INTO research_overlays VALUES(?,?,?,?,?,?)",
            [(r["overlay_id"], r["event_key"], int(r["priority_rank"]), r["reason_code"], r["gene"], json.dumps(r, ensure_ascii=False, sort_keys=True)) for r in research_rows],
        )
        for key, value in scan_stats.items():
            con.execute("INSERT OR REPLACE INTO scan_stats VALUES(?,?)", (str(key), json.dumps(value, ensure_ascii=False, sort_keys=True)))
        for key, value in {
            "schema": ENRICHMENT_SCHEMA, "rules_schema": RULES_SCHEMA, "pipeline_release": PIPELINE_RELEASE,
            "assembly": ASSEMBLY, "created_utc": utcnow(),
        }.items():
            con.execute("INSERT OR REPLACE INTO meta VALUES(?,?)", (key, value))
        con.executescript("""
            CREATE INDEX idx_events_gene ON event_candidates(genes);
            CREATE INDEX idx_open_event ON open_findings(event_key);
            CREATE INDEX idx_research_event ON research_overlays(event_key);
        """)
        con.commit()
        if con.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise RuntimeError("Patient knowledge SQLite integrity check failed")
    finally:
        con.close()
    os.replace(part, path)


def source_usage_from_rows(open_rows: list[dict[str, str]], research_rows: list[dict[str, str]]) -> tuple[set[str], set[str]]:
    open_sources: set[str] = set()
    research_sources: set[str] = set()
    for row in open_rows:
        open_sources.update(x for x in row["evidence_sources"].split(";") if x and x not in {"VEP/CATALOG"})
    for row in research_rows:
        research_sources.update(x for x in row["evidence_sources"].split(";") if x)
    # HGNC canonicalization is applied to every materialized gene event. VEP/catalog is
    # the canonical patient annotation. Mondo is used when a disease relation carries
    # a MONDO identifier whose label can be normalized from the ontology snapshot.
    open_sources.update({"HGNC", "Ensembl VEP"})
    if any(str(row.get("disease_id", "")).upper().startswith("MONDO:") for row in open_rows):
        open_sources.add("Mondo")
    return open_sources, research_sources


def usage_record(source: str, research: bool = False) -> dict[str, Any]:
    if source == "PanelApp":
        return {
            "source": source, "publication_scope": "RESEARCH_ONLY",
            "prohibited_uses": ["commercial", "diagnostic use", "medical decision-making", "providing a healthcare service"],
        }
    if source == "AlphaMissense":
        return {
            "source": source, "publication_scope": "RESEARCH_ONLY",
            "prohibited_uses": ["commercial"],
            "license": "CC BY-NC-SA 4.0",
        }
    if research or source == "CADD":
        return {"source": source, "publication_scope": "RESEARCH_ONLY", "prohibited_uses": ["commercial"]}
    return {"source": source, "publication_scope": "COMMERCIAL_ALLOWED", "prohibited_uses": []}


def write_usage_manifests(output_dir: Path, open_rows: list[dict[str, str]], research_rows: list[dict[str, str]]) -> tuple[Path, Path]:
    open_sources, research_sources = source_usage_from_rows(open_rows, research_rows)
    open_manifest = output_dir / "Dante_knowledge_usage_manifest.json"
    research_manifest = output_dir / "Dante_knowledge_usage_research_manifest.json"
    save_json_atomic(open_manifest, {
        "schema": USAGE_SCHEMA,
        "scope": "OPEN_FINDINGS_ONLY",
        "sources_used": [usage_record(x, False) for x in sorted(open_sources)],
        "note": "Use this manifest with commercial publication checks. Research-only overlays are physically excluded.",
    })
    research_only_sources = research_sources - open_sources
    save_json_atomic(research_manifest, {
        "schema": USAGE_SCHEMA,
        "scope": "OPEN_PLUS_RESEARCH_OVERLAYS",
        "sources_used": [usage_record(x, False) for x in sorted(open_sources)] + [usage_record(x, True) for x in sorted(research_only_sources)],
        "note": "Research-only publication manifest. PanelApp/CADD restrictions apply where those sources are present.",
    })
    return open_manifest, research_manifest


def validate_output(path: Path, state: dict[str, Any]) -> bool:
    if state.get("schema") != ENRICHMENT_SCHEMA or state.get("status") != "complete" or state.get("rules_schema") != RULES_SCHEMA:
        return False
    for key, identity in state.get("outputs", {}).items():
        p = path / str(identity.get("name", ""))
        if not p.is_file() or p.stat().st_size != int(identity.get("length", -1)):
            return False
        expected_sha = str(identity.get("sha256", ""))
        if expected_sha and sha256_file(p) != expected_sha:
            return False
    return True


def state_current(state: dict[str, Any], signatures: dict[str, str], output_dir: Path) -> bool:
    if not validate_output(output_dir, state):
        return False
    saved = state.get("input_signatures", {})
    return all(saved.get(k) == v for k, v in signatures.items())


def input_signatures(paths: dict[str, Path | None], open_manifest: Path, research_manifest: Path, index_state: Path, max_af: float) -> dict[str, str]:
    result = {
        "open_manifest_sha256": sha256_file(open_manifest),
        "research_manifest_sha256": sha256_file(research_manifest),
        "knowledge_index_state_sha256": sha256_file(index_state),
        "rules_schema": RULES_SCHEMA,
        "max_af": f"{max_af:.12g}",
    }
    for key, path in paths.items():
        if path and path.is_file():
            # Large genomic inputs are already immutable pipeline outputs; content hashes are required for reproducibility.
            result[f"{key}_sha256"] = sha256_file(path)
        else:
            result[f"{key}_sha256"] = "MISSING"
    return result


def run(args: argparse.Namespace) -> int:
    root = args.root.resolve()
    patient_root = args.patient_root.resolve()
    open_root = args.open_root.resolve()
    research_root = args.research_root.resolve()
    output_dir = args.output_dir.resolve()
    work = output_dir / "work"
    output_dir.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    open_manifest_path = open_root / "open_reference_manifest.json"
    research_manifest_path = research_root / "research_reference_manifest.json"
    open_manifest = load_json(open_manifest_path)
    research_manifest = load_json(research_manifest_path)
    if open_manifest.get("status") != "complete":
        raise RuntimeError("Open reference bundle is not complete")
    if research_manifest.get("status") != "complete":
        raise RuntimeError("Research reference bundle is not complete")

    index_dir = root / "reference_data" / "knowledge_index"
    index_db = build_reference_index(open_root, research_root, index_dir, args.force_index)
    knowledge = GeneKnowledge(index_db)

    paths = {
        "catalog": args.catalog.resolve(),
        "cnv": args.cnv.resolve() if args.cnv else None,
        "sv": args.sv.resolve() if args.sv else None,
        "reference_targets": args.reference_targets.resolve() if args.reference_targets else None,
        "reference_evidence": args.reference_evidence.resolve() if args.reference_evidence else None,
        "clinvar": args.clinvar.resolve() if args.clinvar else None,
    }
    for required_name in ("catalog", "cnv", "sv"):
        required_path = paths[required_name]
        if required_path is None or not required_path.is_file():
            raise FileNotFoundError(f"Required Stage-24 input is missing ({required_name}): {required_path}")

    index_state = index_dir / "knowledge_index.state.json"
    signatures = input_signatures(paths, open_manifest_path, research_manifest_path, index_state, args.max_af)
    state_path = output_dir / "Dante_knowledge_enrichment.state.json"
    previous = load_json(state_path)
    if not args.force and state_current(previous, signatures, output_dir):
        print("[CURRENT] Stage 24 knowledge enrichment is current for identical genomic/reference inputs.")
        return 0

    started = time.monotonic()
    print("STAGE 24 - FULL GENOME KNOWLEDGE ENRICHMENT")
    print(f"  Pipeline:          v{PIPELINE_RELEASE}")
    print(f"  Patient workspace: {patient_root}")
    print(f"  Catalog:           {paths['catalog']}")
    print(f"  Open bundle:       {open_manifest.get('installed_utc', 'unknown')}")
    print(f"  Research bundle:   {research_manifest.get('installed_utc', 'unknown')}")
    print("  OpenAI:            NOT USED")
    print("")

    events, small_counts = scan_small_catalog(paths["catalog"], knowledge, args.max_af, args.progress_every)
    stage17_counts = add_stage17_fastq_events(
        events, paths["reference_evidence"], paths["reference_targets"], knowledge, paths["clinvar"]
    )
    cnv_events, cnv_counts = scan_structural(paths["cnv"], "CNV", knowledge)
    sv_events, sv_counts = scan_structural(paths["sv"], "SV", knowledge)
    structural_events, structural_counts = merge_cross_caller_structural(cnv_events, sv_events)
    structural_events = {
        key: ev for key, ev in structural_events.items()
        if ev.open_reasons or ev.research_reasons
    }
    structural_counts["knowledge_relevant_after_reconciliation"] = len(structural_events)
    events.update(structural_events)
    print(f"[SCAN COMPLETE] candidate events materialized: {len(events):,}")

    biallelic = mark_possible_biallelic(events, knowledge, args.max_af)
    giab_counts = apply_giab(events, giab_beds(open_manifest, open_root))

    remote_status: dict[str, Any] = {}
    try:
        alpha_counts = apply_alphamissense(events, research_root, research_manifest, args.bcftools_image, work)
        remote_status["AlphaMissense"] = {"status": "complete", "targets": alpha_counts.get("targets", 0), "matched": alpha_counts.get("matched", 0)}
    except Exception as exc:
        if args.strict_remote:
            raise
        remote_status["AlphaMissense"] = {"status": "deferred", "error": f"{type(exc).__name__}: {exc}"}
        print(f"[WARNING] AlphaMissense target query deferred: {exc}", file=sys.stderr)

    try:
        dbsnp_counts = apply_dbsnp(events, open_root, open_manifest, args.bcftools_image, work)
        remote_status["dbSNP"] = {"status": "complete", "targets": dbsnp_counts.get("targets", 0), "known_rsid_targets": dbsnp_counts.get("known_rsid_targets", 0), "matched": dbsnp_counts.get("matched", 0)}
    except Exception as exc:
        if args.strict_remote:
            raise
        remote_status["dbSNP"] = {"status": "deferred", "error": f"{type(exc).__name__}: {exc}"}
        print(f"[WARNING] dbSNP target query deferred: {exc}", file=sys.stderr)

    try:
        cadd_counts = apply_cadd(events, research_root, research_manifest, args.bcftools_image, work)
        remote_status["CADD"] = {"status": "complete", "targets": cadd_counts.get("targets", 0), "matched": cadd_counts.get("matched", 0)}
    except Exception as exc:
        if args.strict_remote:
            raise
        remote_status["CADD"] = {"status": "deferred", "error": f"{type(exc).__name__}: {exc}"}
        print(f"[WARNING] CADD remote target query deferred: {exc}", file=sys.stderr)

    open_rows = open_findings(events, knowledge)
    research_rows = research_overlays(events, open_rows, knowledge)

    open_tsv = output_dir / "Dante_knowledge_findings_open.tsv"
    research_tsv = output_dir / "Dante_knowledge_findings_research.tsv"
    qc_tsv = output_dir / "Dante_knowledge_qc.tsv"
    report = output_dir / "Dante_knowledge_enrichment_report.txt"
    patient_db = output_dir / "Dante_knowledge.sqlite3"

    write_tsv(open_tsv, OPEN_FINDING_COLUMNS, open_rows)
    write_tsv(research_tsv, RESEARCH_COLUMNS, research_rows)
    qc_rows = [
        {"event_key": e.event_key, "event_type": e.event_type, "chrom": e.chrom, "position": str(e.position), "end": str(e.end),
         "flags": ";".join(sorted(e.giab_flags)), "technical_confidence": event_technical_confidence(e)}
        for e in events.values() if e.giab_flags or event_technical_confidence(e).startswith("CAUTION")
    ]
    qc_rows.sort(key=lambda r: (r["chrom"], int(r["position"]), r["event_key"]))
    write_tsv(qc_tsv, QC_COLUMNS, qc_rows)

    scan_stats = {
        "small": dict(small_counts), "stage17": dict(stage17_counts), "cnv": dict(cnv_counts), "sv": dict(sv_counts),
        "structural_reconciliation": dict(structural_counts),
        "candidate_events": len(events), "open_findings": len(open_rows), "research_overlays": len(research_rows),
        "possible_biallelic_unphased_variant_rows": biallelic, "giab": dict(giab_counts), "remote_queries": remote_status,
    }
    build_patient_sqlite(patient_db, events, open_rows, research_rows, scan_stats)
    usage_open, usage_research = write_usage_manifests(output_dir, open_rows, research_rows)

    reason_counts = Counter(r["reason_code"] for r in open_rows)
    research_reason_counts = Counter(r["reason_code"] for r in research_rows)
    deferred = [k for k, v in remote_status.items() if v.get("status") != "complete"]
    report_lines = [
        "DANTE WGS v9.0.0 - STAGE 24 FULL GENOME KNOWLEDGE ENRICHMENT",
        "=" * 78,
        f"Schema: {ENRICHMENT_SCHEMA}",
        f"Rules: {RULES_SCHEMA}",
        f"Assembly: {ASSEMBLY}",
        f"Catalog rows scanned: {small_counts.get('catalog_rows', 0):,}",
        f"Carried PASS catalog rows: {small_counts.get('carried_pass', 0):,}",
        f"Knowledge-relevant candidate events materialized: {len(events):,}",
        f"Open/commercial-safe findings: {len(open_rows):,}",
        f"Research-only overlays: {len(research_rows):,}",
        f"GIAB-caution candidate events: {len(qc_rows):,}",
        f"Possible biallelic-unphased AR variant rows: {biallelic:,}",
        "",
        "OPEN FINDING REASON CODES",
        *[f"  {k}: {v:,}" for k, v in sorted(reason_counts.items())],
        "",
        "RESEARCH-ONLY OVERLAY REASON CODES",
        *[f"  {k}: {v:,}" for k, v in sorted(research_reason_counts.items())],
        "",
        "REMOTE/TABIX TARGET QUERIES",
        *[f"  {k}: {json.dumps(v, ensure_ascii=False, sort_keys=True)}" for k, v in sorted(remote_status.items())],
        "",
        "CONTRACT BOUNDARIES",
        "- The full small-variant catalog was scanned, but only knowledge-relevant events are duplicated into Dante_knowledge.sqlite3.",
        "- Open findings and research-only overlays are physically separate outputs.",
        "- AlphaMissense and CADD are supporting predictors, not clinical classifications.",
        "- PanelApp, AlphaMissense and CADD are RESEARCH ONLY and are excluded from the open findings TSV.",
        "- Structural biological evidence and technical call confidence are separate axes.",
        "- AnnotSV pathogenicity class never proves that a structural call itself is real.",
        "- Reciprocal Manta BND mates and AnnotSV split rows are collapsed before findings are generated.",
        "- DEL/DUP concordance between Manta and DRAGEN is recorded as CROSS_CALLER_CONFIRMED.",
        "- GIAB difficult-region overlap adds GIAB_CAUTION; it never deletes a finding.",
        "- Multiple heterozygous variants in an AR gene are labelled POSSIBLE_BIALLELIC_UNPHASED, never compound heterozygous.",
        "- HPO ontology is vocabulary only in the installed open bundle; no phenotype→gene inference is performed.",
        "- No OpenAI/API evidence research is performed in Stage 24; that begins only in Stage 25.",
        "",
    ]
    if deferred:
        report_lines += [
            "WARNING",
            "The following optional indexed target lookups were deferred: " + ", ".join(deferred) + ".",
            "Local full-genome scanning and source separation completed, but rerun Stage 24 with network/Docker available for the missing overlays.",
            "",
        ]
    report.write_text("\n".join(report_lines), encoding="utf-8")

    outputs = {p.stem: file_identity(p) for p in (open_tsv, research_tsv, qc_tsv, report, patient_db, usage_open, usage_research)}
    state = {
        "schema": ENRICHMENT_SCHEMA,
        "status": "complete" if not deferred or not args.strict_remote else "incomplete",
        "pipeline_release": PIPELINE_RELEASE,
        "rules_schema": RULES_SCHEMA,
        "assembly": ASSEMBLY,
        "completed_utc": utcnow(),
        "duration_seconds": round(time.monotonic() - started, 3),
        "input_signatures": signatures,
        "scan_stats": scan_stats,
        "remote_deferred": deferred,
        "outputs": outputs,
    }
    save_json_atomic(state_path, state)
    print(f"[DONE] Stage 24 | open findings {len(open_rows):,} | research overlays {len(research_rows):,} | {state['duration_seconds']:.1f}s")
    if deferred:
        print("[READY WITH WARNINGS] Deferred target lookups: " + ", ".join(deferred))
    return 0


def status(output_dir: Path) -> int:
    state = load_json(output_dir / "Dante_knowledge_enrichment.state.json")
    if not state:
        print("NOT RUN")
        return 1
    if not validate_output(output_dir, state):
        print("INCOMPLETE")
        return 1
    deferred = state.get("remote_deferred", [])
    if deferred:
        print("READY WITH WARNINGS: " + ", ".join(str(x) for x in deferred))
    else:
        print("READY")
    return 0


def self_test() -> int:
    assert clinvar_bucket("Pathogenic") == "PATHOGENIC"
    assert clinvar_bucket("Pathogenic|Benign") == "CONFLICT"
    assert clinvar_bucket("Likely_pathogenic|protective") == "MIXED"
    assert is_rare("0.001", 0.01)
    assert not is_rare("0.1", 0.01)
    assert rarity_status(".", 0.01) == "UNKNOWN"
    assert not is_rare(".", 0.01)
    assert zygosity("0/1", "A", "G") == "HET_ALT"
    assert inheritance_bucket("Autosomal recessive") == "AR"
    assert deterministic_id("X", "a", "b") == deterministic_id("X", "a", "b")
    print("SELF-TEST PASS | v9.0.0 hotfix5 knowledge rules-v7 + research-only AlphaMissense + structural caller reconciliation")
    return 0


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Dante WGS v9 full-genome knowledge enrichment")
    sub = ap.add_subparsers(dest="mode", required=True)
    runp = sub.add_parser("run")
    runp.add_argument("--root", type=Path, required=True)
    runp.add_argument("--patient-root", type=Path, required=True)
    runp.add_argument("--open-root", type=Path, required=True)
    runp.add_argument("--research-root", type=Path, required=True)
    runp.add_argument("--catalog", type=Path, required=True)
    runp.add_argument("--cnv", type=Path, required=True)
    runp.add_argument("--sv", type=Path, required=True)
    runp.add_argument("--reference-targets", type=Path)
    runp.add_argument("--reference-evidence", type=Path)
    runp.add_argument("--clinvar", type=Path)
    runp.add_argument("--output-dir", type=Path, required=True)
    runp.add_argument("--bcftools-image", required=True)
    runp.add_argument("--max-af", type=float, default=DEFAULT_MAX_AF)
    runp.add_argument("--progress-every", type=int, default=500_000)
    runp.add_argument("--strict-remote", action="store_true", help="Fail Stage 24 if AlphaMissense/dbSNP/CADD target queries cannot complete")
    runp.add_argument("--force", action="store_true")
    runp.add_argument("--force-index", action="store_true")
    statusp = sub.add_parser("status")
    statusp.add_argument("--output-dir", type=Path, required=True)
    sub.add_parser("self-test")
    return ap


def main() -> int:
    args = parser().parse_args()
    if args.mode == "self-test":
        return self_test()
    if args.mode == "status":
        return status(args.output_dir.resolve())
    return run(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nERROR: cancelled by user. Published outputs are atomic; completed reference index is reusable.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        traceback.print_exc()
        raise SystemExit(1)
