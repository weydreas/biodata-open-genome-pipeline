#!/usr/bin/env python3
"""Install non-commercial / research-only reference knowledge for human WGS research."""
from __future__ import annotations
import argparse, datetime as dt, gzip, hashlib, json, os, random, shutil, sqlite3, subprocess, sys, time
from pathlib import Path
from dataclasses import dataclass
from typing import Any
import urllib.error, urllib.request

PIPELINE_RELEASE = "8.5.2"
BUNDLE_SCHEMA = "research-reference-bundle-v1"
USER_AGENT = f"DanteWGS/{PIPELINE_RELEASE} research-reference-manager"
DOWNLOAD_TIMEOUT = 120
MAX_ATTEMPTS = 5
CHUNK = 4 * 1024 * 1024
RESEARCH_LICENSE_CLASS = "RESEARCH_ONLY_NONCOMMERCIAL"

# PanelApp snapshotting is intentionally conservative: this is background reference
# installation, so avoiding server throttling is more important than burst speed.
PANELAPP_MIN_INTERVAL_SECONDS = 1.0
PANELAPP_429_DEFAULT_WAIT_SECONDS = 60.0
PANELAPP_429_MAX_WAIT_SECONDS = 300.0
PANELAPP_MAX_ATTEMPTS = 12
_panelapp_last_request_monotonic = 0.0

PANELAPP_BASE = "https://panelapp.genomicsengland.co.uk/api/v1"
PANELAPP_TERMS = (
    "https://prod-media-panelapp.genomicsengland.co.uk/media/files/"
    "GEL_-_PanelApp_Terms_of_Use_December_2019.pdf"
)
PANELAPP_ATTR = (
    "Genomics England PanelApp. PanelApp and identified contributors must be acknowledged. "
    "Use is restricted by the PanelApp Terms of Use."
)

CADD_BASE = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37"
CADD_SNV_URL = f"{CADD_BASE}/whole_genome_SNVs.tsv.gz"
CADD_TBI_URL = f"{CADD_SNV_URL}.tbi"
CADD_TBI_MD5 = "6bc06b08fd9a521d05d970a06939e52d"
CADD_LICENSE_URL = "https://cadd.bihealth.org/impressum"
CADD_ATTR = (
    "CADD v1.7 (Combined Annotation Dependent Depletion). CADD scores are free for "
    "non-commercial applications; commercial use requires a separate licence."
)

@dataclass(frozen=True)
class Resource:
    source: str
    key: str
    url: str
    relative_path: str
    license_name: str
    license_url: str
    attribution: str
    content_kind: str = "file"
    expected_md5: str | None = None
    notes: str = ""
    prohibited_uses: tuple[str, ...] = ("commercial",)

RESOURCES: tuple[Resource, ...] = (
    Resource(
        "CADD", "grch37_v1_7_snv_tabix",
        CADD_TBI_URL,
        "CADD/GRCh37_v1.7_whole_genome_SNVs.tsv.gz.tbi",
        "CADD-NONCOMMERCIAL", CADD_LICENSE_URL, CADD_ATTR,
        content_kind="binary", expected_md5=CADD_TBI_MD5,
        notes=(
            "Only the small tabix index is installed. The ~79 GiB score file remains remote; "
            "a later enrichment stage can query only requested GRCh37 ranges/variants."
        )
    ),
)

def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()

def fmt_bytes(value: int | float | None) -> str:
    if value is None: return "?"
    x = float(value)
    for unit in ("B","KiB","MiB","GiB","TiB"):
        if x < 1024 or unit == "TiB":
            return f"{x:.1f} {unit}"
        x /= 1024
    return f"{x:.1f} TiB"

def file_hash(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def request(url: str, *, headers: dict[str,str] | None = None, method: str = "GET"):
    merged = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if headers: merged.update(headers)
    return urllib.request.Request(url, headers=merged, method=method)


def curl_executable() -> str | None:
    """System curl fallback; certificate verification remains enabled."""
    return shutil.which("curl.exe") or shutil.which("curl")


def curl_download_verified(url: str, part: Path, *, resume: bool = False) -> None:
    curl = curl_executable()
    if not curl:
        raise RuntimeError("Verified curl fallback is unavailable")
    cmd = [
        curl, "--fail", "--location", "--show-error", "--silent",
        "--retry", "3", "--retry-delay", "2",
        "--connect-timeout", "30", "--max-time", "600",
        "--user-agent", USER_AGENT,
    ]
    if resume and part.exists() and part.stat().st_size > 0:
        cmd += ["--continue-at", "-"]
    cmd += ["--output", str(part), url]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"curl TLS-verified download failed: {proc.stderr.strip()}")


def curl_json_verified(url: str) -> Any:
    curl = curl_executable()
    if not curl:
        raise RuntimeError("Verified curl fallback is unavailable")
    cmd = [
        curl, "--fail", "--location", "--show-error", "--silent",
        "--retry", "3", "--retry-delay", "2",
        "--connect-timeout", "30", "--max-time", "180",
        "--user-agent", USER_AGENT,
        "--header", "Accept: application/json", url,
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(f"curl TLS-verified JSON request failed: {proc.stderr.strip()}")
    return json.loads(proc.stdout)

def load_json(path: Path) -> dict[str,Any]:
    if not path.is_file(): return {}
    try: return json.loads(path.read_text(encoding="utf-8"))
    except Exception: return {}

def save_json_atomic(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(part, path)

def validate_resource(r: Resource, path: Path) -> None:
    if not path.is_file() or path.stat().st_size <= 0:
        raise RuntimeError(f"Missing/empty resource: {path}")
    if r.content_kind == "gzip":
        with gzip.open(path, "rb") as f:
            if not f.read(64): raise RuntimeError(f"Empty gzip payload: {path}")
    if r.expected_md5:
        got = file_hash(path, "md5").lower()
        if got != r.expected_md5.lower():
            raise RuntimeError(f"MD5 mismatch for {r.source}/{r.key}: {got}")

def download(r: Resource, output: Path, previous: dict[str,Any], force: bool) -> dict[str,Any]:
    dest = output / r.relative_path
    dest.parent.mkdir(parents=True, exist_ok=True)
    prev = previous.get("artifacts", {}).get(r.key, {}) if previous else {}
    if dest.is_file() and not force and prev:
        if prev.get("length") == dest.stat().st_size and prev.get("sha256") == file_hash(dest):
            print(f"[CURRENT] {r.source}/{r.key} | {fmt_bytes(dest.stat().st_size)}")
            return dict(prev)

    part = dest.with_name(dest.name + ".part")
    resume = part.stat().st_size if part.exists() and not force else 0
    if force and part.exists():
        part.unlink()
        resume = 0

    last = None
    hdr: dict[str,Any] = {}
    transport = "urllib"
    completed = False
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            headers = {"Range": f"bytes={resume}-"} if resume else {}
            resp = urllib.request.urlopen(request(r.url, headers=headers), timeout=DOWNLOAD_TIMEOUT)
            with resp:
                status = getattr(resp, "status", 200)
                if resume and status != 206:
                    resume = 0
                mode = "ab" if resume and status == 206 else "wb"
                if mode == "wb" and part.exists():
                    part.unlink()
                with part.open(mode) as out:
                    while True:
                        block = resp.read(CHUNK)
                        if not block:
                            break
                        out.write(block)
                hdr = dict(resp.headers.items())
            transport = "urllib"
            completed = True
            break
        except Exception as exc:
            last = exc
            if "CERTIFICATE_VERIFY_FAILED" in str(exc) and curl_executable():
                try:
                    print(
                        f"[TLS] Python certificate validation failed for {r.source}/{r.key}; "
                        "retrying with system curl and normal TLS verification.",
                        file=sys.stderr,
                    )
                    curl_download_verified(r.url, part, resume=bool(resume))
                    transport = "curl-system-trust"
                    hdr = {}
                    completed = True
                    break
                except Exception as curl_exc:
                    last = curl_exc
            print(f"[RETRY] {r.source}/{r.key} {attempt}/{MAX_ATTEMPTS}: {last}", file=sys.stderr)
            if attempt < MAX_ATTEMPTS:
                time.sleep(min(30, 2 ** attempt))

    if not completed:
        raise RuntimeError(f"Failed to download {r.source}/{r.key}: {last}")
    if not part.is_file() or part.stat().st_size <= 0:
        raise RuntimeError(f"Download completed without a usable file: {r.source}/{r.key}")

    os.replace(part, dest)
    validate_resource(r, dest)
    rec = {
        "path": r.relative_path.replace("\\","/"),
        "url": r.url,
        "length": dest.stat().st_size,
        "sha256": file_hash(dest),
        "md5": file_hash(dest,"md5") if r.expected_md5 else None,
        "etag": hdr.get("ETag"),
        "last_modified": hdr.get("Last-Modified"),
        "transport": transport,
        "tls_verification": "enabled",
        "verified_utc": utcnow(),
    }
    print(f"[SAVED]   {r.source}/{r.key} | {fmt_bytes(dest.stat().st_size)} | transport={transport}")
    return rec


def _retry_after_seconds(exc: BaseException, attempt: int) -> float:
    """Honor HTTP Retry-After; otherwise use a conservative exponential cooldown."""
    if isinstance(exc, urllib.error.HTTPError):
        raw = exc.headers.get("Retry-After") if exc.headers else None
        if raw:
            raw = str(raw).strip()
            try:
                return max(1.0, min(PANELAPP_429_MAX_WAIT_SECONDS, float(raw)))
            except ValueError:
                try:
                    from email.utils import parsedate_to_datetime
                    when = parsedate_to_datetime(raw)
                    if when.tzinfo is None:
                        when = when.replace(tzinfo=dt.timezone.utc)
                    seconds = (when - dt.datetime.now(dt.timezone.utc)).total_seconds()
                    return max(1.0, min(PANELAPP_429_MAX_WAIT_SECONDS, seconds))
                except Exception:
                    pass
    base = PANELAPP_429_DEFAULT_WAIT_SECONDS * (1.5 ** max(0, attempt - 1))
    return min(PANELAPP_429_MAX_WAIT_SECONDS, base + random.uniform(0.0, 3.0))


def _panelapp_rate_limit() -> None:
    global _panelapp_last_request_monotonic
    now = time.monotonic()
    elapsed = now - _panelapp_last_request_monotonic
    wait = PANELAPP_MIN_INTERVAL_SECONDS - elapsed
    if wait > 0:
        time.sleep(wait)
    _panelapp_last_request_monotonic = time.monotonic()


def fetch_json(url: str) -> Any:
    """PanelApp JSON fetch with proactive pacing and HTTP 429 recovery."""
    last = None
    for attempt in range(1, PANELAPP_MAX_ATTEMPTS + 1):
        _panelapp_rate_limit()
        try:
            with urllib.request.urlopen(
                request(url, headers={"Accept":"application/json"}),
                timeout=DOWNLOAD_TIMEOUT
            ) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code == 429:
                wait = _retry_after_seconds(exc, attempt)
                retry_after = exc.headers.get("Retry-After") if exc.headers else None
                suffix = f" (server Retry-After={retry_after})" if retry_after else ""
                print(
                    f"[RATE LIMIT] PanelApp HTTP 429; snapshot is preserved. "
                    f"Waiting {wait:.0f}s before retry {attempt}/{PANELAPP_MAX_ATTEMPTS}{suffix}",
                    file=sys.stderr,
                    flush=True,
                )
                time.sleep(wait)
                continue
            if attempt < PANELAPP_MAX_ATTEMPTS:
                wait = min(30.0, 2.0 ** min(attempt, 5))
                print(
                    f"[RETRY] PanelApp HTTP {exc.code}; waiting {wait:.0f}s "
                    f"before retry {attempt}/{PANELAPP_MAX_ATTEMPTS}",
                    file=sys.stderr,
                    flush=True,
                )
                time.sleep(wait)
                continue
        except Exception as exc:
            last = exc
            if "CERTIFICATE_VERIFY_FAILED" in str(exc) and curl_executable():
                try:
                    print(
                        "[TLS] Python certificate validation failed for PanelApp API; "
                        "retrying with system curl and normal TLS verification.",
                        file=sys.stderr,
                    )
                    _panelapp_rate_limit()
                    return curl_json_verified(url)
                except Exception as curl_exc:
                    last = curl_exc
            if attempt < PANELAPP_MAX_ATTEMPTS:
                wait = min(30.0, 2.0 ** min(attempt, 5))
                print(
                    f"[RETRY] PanelApp API error; waiting {wait:.0f}s "
                    f"before retry {attempt}/{PANELAPP_MAX_ATTEMPTS}: {last}",
                    file=sys.stderr,
                    flush=True,
                )
                time.sleep(wait)
                continue
        break
    raise RuntimeError(f"PanelApp API failed after {PANELAPP_MAX_ATTEMPTS} attempts: {url}: {last}")


def panelapp_snapshot(output: Path, force: bool) -> dict[str,Any]:
    """Resumable SQLite snapshot of all current human PanelApp panels."""
    pdir = output / "PanelApp"; pdir.mkdir(parents=True, exist_ok=True)
    db = pdir / "panelapp_snapshot.sqlite3"
    con = sqlite3.connect(db)
    try:
        con.execute("""CREATE TABLE IF NOT EXISTS panels(
            panel_id INTEGER PRIMARY KEY, name TEXT, version TEXT,
            fetched_utc TEXT NOT NULL, payload_json TEXT NOT NULL)""")
        con.execute("""CREATE TABLE IF NOT EXISTS metadata(
            key TEXT PRIMARY KEY, value TEXT NOT NULL)""")
        con.commit()

        url = f"{PANELAPP_BASE}/panels/?page=1"
        summaries = []
        while url:
            page = fetch_json(url)
            if isinstance(page, dict):
                results, url = page.get("results", []), page.get("next")
            elif isinstance(page, list):
                results, url = page, None
            else:
                raise RuntimeError("Unexpected PanelApp API response")
            summaries.extend(x for x in results if isinstance(x, dict))
        print(f"[PANELAPP] current panels listed: {len(summaries):,}")

        existing = {int(a): str(b or "") for a,b in con.execute("SELECT panel_id,version FROM panels")}
        if existing:
            print(f"[PANELAPP] resumable snapshot already contains {len(existing):,} panels")
        fetched = reused = 0
        current_ids = set()
        for summary in summaries:
            if summary.get("id") is None: continue
            pid = int(summary["id"]); current_ids.add(pid)
            version = str(summary.get("version") or summary.get("current_version") or "")
            if not force and existing.get(pid) == version:
                reused += 1; continue
            detail = fetch_json(f"{PANELAPP_BASE}/panels/{pid}/")
            name = str((detail.get("name") if isinstance(detail,dict) else None) or summary.get("name") or "")
            con.execute(
                "INSERT OR REPLACE INTO panels(panel_id,name,version,fetched_utc,payload_json) VALUES(?,?,?,?,?)",
                (pid, name, version, utcnow(), json.dumps(detail, ensure_ascii=False, sort_keys=True))
            )
            # Persist every successful panel so an interrupted/rate-limited run resumes
            # without losing even the last few completed requests.
            con.commit()
            fetched += 1
            if fetched % 25 == 0:
                print(f"[PANELAPP] fetched/updated {fetched:,}; reused {reused:,}", flush=True)
        con.commit()
        if current_ids:
            marks = ",".join("?" for _ in current_ids)
            con.execute(f"DELETE FROM panels WHERE panel_id NOT IN ({marks})", tuple(sorted(current_ids)))
        con.execute("INSERT OR REPLACE INTO metadata(key,value) VALUES('snapshot_utc',?)",(utcnow(),))
        con.execute("INSERT OR REPLACE INTO metadata(key,value) VALUES('api_base',?)",(PANELAPP_BASE,))
        con.commit()
        count = int(con.execute("SELECT COUNT(*) FROM panels").fetchone()[0])
    finally:
        con.close()
    if count <= 0: raise RuntimeError("PanelApp snapshot contains no panels")
    print(f"[PANELAPP] snapshot complete: {count:,} panels | fetched/updated {fetched:,} | reused {reused:,}")
    return {
        "path": str(db.relative_to(output)).replace("\\","/"),
        "length": db.stat().st_size,
        "sha256": file_hash(db),
        "panel_count": count,
        "snapshot_utc": utcnow(),
        "api_base": PANELAPP_BASE,
    }

def head_metadata(url: str) -> dict[str,Any]:
    try:
        with urllib.request.urlopen(request(url, method="HEAD"), timeout=30) as r:
            length = r.headers.get("Content-Length")
            return {
                "url": url, "final_url": r.geturl(),
                "content_length": int(length) if length and length.isdigit() else None,
                "etag": r.headers.get("ETag"),
                "last_modified": r.headers.get("Last-Modified"),
                "http_status": getattr(r,"status",200),
            }
    except Exception as exc:
        return {"url":url, "head_error":f"{type(exc).__name__}: {exc}"}

def write_attributions(output: Path, manifest: dict[str,Any]) -> None:
    lines = [
        "# Research-only knowledge attribution and restriction ledger","",
        "> **RESEARCH ONLY — NO COMMERCIAL USE.**","",
        "> PanelApp has an additional restriction: absent a separate agreement with Genomics England, "
        "it must not be used for diagnostic use, medical decision-making, or providing a healthcare service, "
        "even on a non-commercial basis.","",
        f"Snapshot installed: `{manifest['installed_utc']}`","",
    ]
    for source in sorted(manifest["sources"]):
        s = manifest["sources"][source]
        lines += [
            f"## {source}","",
            f"- License/terms: `{s['license_name']}`",
            f"- Terms URL: {s['license_url']}",
            f"- Publication scope: `{s['publication_scope']}`",
            f"- Attribution: {s['attribution']}",
        ]
        if s.get("prohibited_uses"):
            lines.append("- Prohibited without separate permission: " + ", ".join(s["prohibited_uses"]))
        if s.get("notes"): lines.append(f"- Notes: {s['notes']}")
        lines.append("")
    (output / "ATTRIBUTIONS_RESEARCH_ONLY.md").write_text("\n".join(lines), encoding="utf-8")

def install(root: Path, output: Path, force: bool) -> int:
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "research_reference_manifest.json"
    previous = load_json(manifest_path)
    if shutil.disk_usage(output).free < 1024**3:
        raise RuntimeError("At least 1 GiB free disk is required")

    print("RESEARCH-ONLY HUMAN REFERENCE KNOWLEDGE BUNDLE")
    print("  Publication:           RESEARCH ONLY")
    print("  Commercial use:        NOT ALLOWED")
    print("  Species scope:         Homo sapiens only")
    print("  PanelApp:              no diagnostic/medical-decision/healthcare use without separate GEL agreement")
    print("  CADD 79 GiB scores:    remote only; local tabix index only")
    print("  UCSC Chain Files:       NOT auto-installed; current EULA is narrower than generic non-commercial use")
    print("")

    grouped = {}
    for r in RESOURCES:
        g = grouped.setdefault(r.source, {
            "license_name": r.license_name,
            "license_class": RESEARCH_LICENSE_CLASS,
            "license_url": r.license_url,
            "publication_scope": "RESEARCH_ONLY",
            "commercial_allowed": False,
            "species_scope": "Homo sapiens",
            "attribution": r.attribution,
            "prohibited_uses": sorted(set(r.prohibited_uses)),
            "notes": r.notes,
            "artifacts": {},
        })
        prev = previous.get("sources",{}).get(r.source,{}) if previous else {}
        g["artifacts"][r.key] = download(r, output, prev, force)

    grouped["PanelApp"] = {
        "license_name": "GENOMICS-ENGLAND-PANELAPP-TERMS",
        "license_class": RESEARCH_LICENSE_CLASS,
        "license_url": PANELAPP_TERMS,
        "publication_scope": "RESEARCH_ONLY",
        "commercial_allowed": False,
        "species_scope": "Homo sapiens",
        "attribution": PANELAPP_ATTR,
        "prohibited_uses": [
            "commercial","commercial research","diagnostic use",
            "medical decision-making","providing a healthcare service"
        ],
        "notes": "Installed only for non-commercial research/hypothesis generation.",
        "artifacts": {"snapshot_sqlite": panelapp_snapshot(output, force)},
    }

    cadd_remote = output / "CADD" / "remote_score_source.json"
    save_json_atomic(cadd_remote, {
        "model":"CADD v1.7 GRCh37",
        "score_file_url":CADD_SNV_URL,
        "remote_metadata":head_metadata(CADD_SNV_URL),
        "local_tabix_md5":CADD_TBI_MD5,
        "strategy":"Future enrichment should query only required ranges/variants; do not download 79 GiB by default."
    })
    grouped["CADD"]["artifacts"]["remote_score_source"] = {
        "path":str(cadd_remote.relative_to(output)).replace("\\","/"),
        "length":cadd_remote.stat().st_size,
        "sha256":file_hash(cadd_remote),
        "created_utc":utcnow(),
    }

    manifest = {
        "schema":BUNDLE_SCHEMA, "status":"complete", "pipeline_release":PIPELINE_RELEASE,
        "installed_utc":utcnow(), "bundle_scope":"RESEARCH_ONLY",
        "commercial_allowed":False, "species_scope":"Homo sapiens",
        "report_policy":{
            "allowed_report_types":["RESEARCH"],
            "commercial_report":"BLOCK",
            "future_usage_manifest":"Dante_results/knowledge_usage_manifest.json",
            "research_banner":"RESEARCH ONLY — NOT FOR COMMERCIAL USE",
        },
        "sources":grouped,
        "optional_not_installed":{
            "OMIM":{
                "credential":"OMIM_API_KEY",
                "commercial_default":"BLOCK",
                "reason":"Separate credential and custom OMIM terms are required; OMIM is not redistributed or auto-downloaded."
            },
            "UCSC_Chain_Files":{
                "commercial_default":"BLOCK",
                "reason":(
                    "Not auto-installed. Current UCSC Chain File terms make them free only for "
                    "Independent Researchers and Nonprofit Organizations; use at the direction of "
                    "or for the benefit of a for-profit organization requires a commercial license. "
                    "GRCh37/GRCh38 rsID reconciliation will instead use the open dbSNP dual-build layer."
                ),
                "terms_url":"https://genome.ucsc.edu/license/EULA.pdf"
            }
        }
    }
    write_attributions(output, manifest)
    attr = output / "ATTRIBUTIONS_RESEARCH_ONLY.md"
    manifest["attributions"] = {
        "path":attr.name, "length":attr.stat().st_size, "sha256":file_hash(attr)
    }
    save_json_atomic(manifest_path, manifest)
    print("")
    print("[DONE] Research-only knowledge bundle installed")
    print(f"       Manifest:     {manifest_path}")
    print(f"       Restrictions: {attr}")
    print("       Sources:      " + ", ".join(sorted(grouped)))
    print("       OMIM:         optional gated connector; user API key/terms required")
    print("       UCSC chains:  intentionally excluded; dbSNP dual-build will handle target coordinate reconciliation")
    return 0

def status(output: Path) -> int:
    m = load_json(output / "research_reference_manifest.json")
    if not m: print("NOT RUN"); return 1
    if m.get("schema") != BUNDLE_SCHEMA or m.get("status") != "complete":
        print("INCOMPLETE"); return 1
    missing = []
    for source, meta in m.get("sources",{}).items():
        for key, art in meta.get("artifacts",{}).items():
            rel = art.get("path")
            if not rel: continue
            p = output / rel
            if not p.is_file() or p.stat().st_size != int(art.get("length",-1)):
                missing.append(f"{source}/{key}")
    if missing:
        print("INCOMPLETE: " + ", ".join(missing[:8])); return 1
    print("READY"); return 0

def self_test() -> int:
    assert len({(r.source,r.key) for r in RESOURCES}) == len(RESOURCES)
    assert {r.source for r in RESOURCES} == {"CADD"}
    assert next(r for r in RESOURCES if r.source=="CADD").expected_md5 == CADD_TBI_MD5
    assert "panelapp.genomicsengland.co.uk" in PANELAPP_BASE
    assert PANELAPP_MIN_INTERVAL_SECONDS >= 0.5
    assert PANELAPP_MAX_ATTEMPTS >= 8
    assert PANELAPP_429_DEFAULT_WAIT_SECONDS >= 30
    print("SELF-TEST PASS | PanelApp + CADD research bundle + rate-limit policy")
    print("TLS policy: verification always enabled; Python CA failure may fall back to system curl trust store")
    print("UCSC Chain Files: excluded from automatic bundle under current EULA")
    return 0

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("install","status","self-test"))
    ap.add_argument("--root", type=Path, default=Path.cwd())
    ap.add_argument("--output", type=Path)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    root = a.root.resolve()
    output = (a.output or root/"reference_data"/"research_reference").resolve()
    if a.mode=="self-test": return self_test()
    if a.mode=="status": return status(output)
    return install(root, output, a.force)

if __name__ == "__main__":
    try: raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted; completed files and PanelApp SQLite rows are preserved.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
