#!/usr/bin/env python3
"""Local IGV.js genome browser for the Dante WGS pipeline.

Security/privacy model:
- binds to 127.0.0.1 only;
- exposes only a fixed allowlist of project genomic files;
- never accepts arbitrary filesystem paths from HTTP requests;
- genomic data stay on localhost;
- IGV.js is downloaded once from the pinned official jsDelivr package and then served locally.

The current Stage 15 BAM is target-only. Empty regions between target loci therefore
must not be interpreted as absence of WGS coverage.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import sys
import threading
import time
import urllib.request
import webbrowser
from typing import Any

PIPELINE_RELEASE = "9.0.0"
IGV_VERSION = "3.8.5"
IGV_URL = f"https://cdn.jsdelivr.net/npm/igv@{IGV_VERSION}/dist/igv.min.js"
IGV_LICENSE = "MIT"
IGV_LICENSE_URL = "https://github.com/igvteam/igv.js/blob/master/LICENSE"
IGV_CITATION = (
    "Robinson JT, Thorvaldsdottir H, Turner D, Mesirov JP. "
    "igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer. "
    "Bioinformatics. 2023;39(1):btac830."
)
DEFAULT_PORT = 8765

TARGET_TSV_NAME = "Dante_reference_targets.tsv"
BAM_NAME = "Dante_reference_targets.bam"
CALLS_NAME = "Dante_reference_calls.vcf.gz"
CLINVAR_NAME = "clinvar.vcf.gz"
FASTA_NAME = "Homo_sapiens.GRCh37.dna.primary_assembly.fa"

INDEX_SCHEMA = "browser-rsid-index-v1"

HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Dante Local Genome Browser</title>
  <style>
    :root {
      color-scheme: light dark;
      --bg: #111318;
      --panel: #1a1e25;
      --text: #e8edf2;
      --muted: #9aa6b2;
      --border: #303843;
      --accent: #6ea8fe;
      --warn: #ffcc66;
      --ok: #7bd88f;
      --danger: #ff7b72;
    }
    @media (prefers-color-scheme: light) {
      :root {
        --bg: #f7f8fa;
        --panel: #ffffff;
        --text: #1d232a;
        --muted: #66717c;
        --border: #d8dde3;
        --accent: #0b63ce;
        --warn: #8a5a00;
        --ok: #147a2d;
        --danger: #b42318;
      }
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: system-ui, -apple-system, Segoe UI, Arial, sans-serif;
      background: var(--bg);
      color: var(--text);
    }
    header {
      padding: 12px 16px;
      border-bottom: 1px solid var(--border);
      background: var(--panel);
      position: sticky;
      top: 0;
      z-index: 20;
    }
    .row {
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
    }
    h1 {
      margin: 0 14px 0 0;
      font-size: 18px;
      font-weight: 650;
    }
    input {
      min-width: 330px;
      flex: 1;
      max-width: 720px;
      padding: 8px 10px;
      border-radius: 6px;
      border: 1px solid var(--border);
      background: var(--bg);
      color: var(--text);
      font-size: 14px;
    }
    button {
      padding: 8px 11px;
      border-radius: 6px;
      border: 1px solid var(--border);
      background: var(--panel);
      color: var(--text);
      cursor: pointer;
    }
    button:hover { border-color: var(--accent); }
    .stop { color: var(--danger); }
    .note {
      margin-top: 8px;
      font-size: 12px;
      color: var(--muted);
    }
    .note strong { color: var(--warn); }
    #status {
      padding: 7px 16px;
      font-size: 12px;
      color: var(--muted);
      border-bottom: 1px solid var(--border);
    }
    .ok { color: var(--ok); }
    .warn { color: var(--warn); }
    #search-result {
      padding: 0 16px;
      min-height: 0;
      font-size: 12px;
      color: var(--muted);
    }
    #igv-div {
      background: white;
      min-height: calc(100vh - 120px);
      padding: 4px;
    }
  </style>
</head>
<body>
<header>
  <div class="row">
    <h1>Dante Local Genome Browser</h1>
    <input id="query" placeholder="rsID or locus, e.g. rs7412 or 19:45411941-45412100">
    <button id="go">Go</button>
    <button id="stop" class="stop">Stop server</button>
  </div>
  <div class="note">
    <strong>Target-only BAM:</strong> the Stage 15 alignment shown here contains reads retained around the reference target panel.
    Empty inter-target regions do not mean that the original WGS FASTQ had no coverage there.
  </div>
</header>
<div id="status">Loading configuration…</div>
<div id="search-result"></div>
<div id="igv-div"></div>
<script src="/vendor/igv.min.js"></script>
<script>
let browser = null;

function setStatus(message, klass="") {
  const el = document.getElementById("status");
  el.textContent = message;
  el.className = klass;
}

async function init() {
  try {
    const cfg = await (await fetch("/api/config", {cache: "no-store"})).json();
    setStatus("Initializing IGV.js…");
    browser = await igv.createBrowser(document.getElementById("igv-div"), cfg.igv);
    const tracks = cfg.track_names.join(" • ");
    setStatus(`Ready • ${cfg.reference_name} • tracks: ${tracks} • rsID index: ${cfg.index_status}`, "ok");
    pollIndex();
  } catch (e) {
    console.error(e);
    setStatus("Browser initialization failed: " + e, "warn");
  }
}

async function pollIndex() {
  try {
    const s = await (await fetch("/api/status", {cache: "no-store"})).json();
    if (s.index_status === "BUILDING") {
      setStatus(`Ready • rsID index building: ${s.index_rows.toLocaleString()} rows`, "warn");
      setTimeout(pollIndex, 1500);
    } else if (s.index_status === "READY") {
      setStatus(`Ready • rsID search index: ${s.index_rows.toLocaleString()} loci • localhost only`, "ok");
    }
  } catch (_) {}
}

function normalizedCoordinateQuery(q) {
  q = q.trim();
  const m = q.match(/^chr?([0-9]+|X|Y|MT|M):([0-9,]+)(?:-([0-9,]+))?$/i);
  if (!m) return null;
  const chr = m[1].toUpperCase() === "M" ? "MT" : m[1].toUpperCase();
  const start = m[2].replaceAll(",", "");
  const end = m[3] ? m[3].replaceAll(",", "") : null;
  return end ? `${chr}:${start}-${end}` : `${chr}:${start}`;
}

async function go() {
  if (!browser) return;
  const input = document.getElementById("query");
  const q = input.value.trim();
  const result = document.getElementById("search-result");
  result.textContent = "";
  if (!q) return;

  const coord = normalizedCoordinateQuery(q);
  if (coord) {
    await browser.search(coord);
    return;
  }

  if (/^rs\d+$/i.test(q)) {
    const r = await fetch("/api/search?q=" + encodeURIComponent(q), {cache: "no-store"});
    const payload = await r.json();
    if (!r.ok) {
      result.textContent = payload.error || "rsID lookup failed";
      return;
    }
    result.textContent = `${payload.rsid} → ${payload.chromosome}:${payload.position}  REF=${payload.ref || "."}  ALT=${payload.alts || "."}`;
    await browser.search(payload.locus);
    return;
  }

  // Let IGV try any other syntax it understands.
  await browser.search(q);
}

document.getElementById("go").addEventListener("click", go);
document.getElementById("query").addEventListener("keydown", e => {
  if (e.key === "Enter") go();
});
document.getElementById("stop").addEventListener("click", async () => {
  document.getElementById("stop").disabled = true;
  setStatus("Stopping local server…");
  try { await fetch("/api/shutdown", {method: "POST"}); } catch (_) {}
  setTimeout(() => window.close(), 250);
});

init();
</script>
</body>
</html>
"""

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def ensure_igv(root: Path) -> Path:
    runtime = root / "reference_data" / "browser_runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    target = runtime / f"igv-{IGV_VERSION}.min.js"
    manifest = runtime / f"igv-{IGV_VERSION}.manifest.json"

    if target.is_file() and target.stat().st_size > 1_000_000:
        return target

    part = target.with_suffix(target.suffix + ".part")
    print(f"[BROWSER] Downloading IGV.js {IGV_VERSION} once from official jsDelivr package…")
    req = urllib.request.Request(
        IGV_URL,
        headers={"User-Agent": f"DanteWGS/{PIPELINE_RELEASE}", "Accept": "application/javascript,*/*"},
    )
    last: Exception | None = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(req, timeout=60) as response, part.open("wb") as out:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            if part.stat().st_size < 1_000_000:
                raise RuntimeError(f"Downloaded IGV.js is unexpectedly small: {part.stat().st_size} bytes")
            os.replace(part, target)
            payload = {
                "source": "igv.js",
                "version": IGV_VERSION,
                "url": IGV_URL,
                "sha256": sha256_file(target),
                "length": target.stat().st_size,
                "license": IGV_LICENSE,
                "license_url": IGV_LICENSE_URL,
                "citation": IGV_CITATION,
                "downloaded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
            manifest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            return target
        except Exception as exc:
            last = exc
            try:
                part.unlink()
            except FileNotFoundError:
                pass
            if attempt < 3:
                time.sleep(attempt * 2)
    raise RuntimeError(
        "Could not download the pinned IGV.js runtime. Internet is required once; "
        f"after that it is cached locally. Source: {IGV_URL}. Last error: {last}"
    )

def first_existing(*paths: Path) -> Path | None:
    for p in paths:
        if p.is_file() and p.stat().st_size > 0:
            return p
    return None

def resolve_artifacts(root: Path, patient_root: Path) -> dict[str, Path]:
    work = patient_root / "Dante_work"
    qc = patient_root / "Dante_qc"
    ann = root / "annotation_data"
    fasta = ann / "vep_rsid" / "fasta" / FASTA_NAME
    bam = work / BAM_NAME
    calls = work / CALLS_NAME
    clinvar = ann / "clinvar_grch37" / CLINVAR_NAME
    targets = qc / TARGET_TSV_NAME

    artifacts: dict[str, Path] = {}
    required = {
        "reference_fasta": fasta,
        "reference_fai": Path(str(fasta) + ".fai"),
        "bam": bam,
        "bai": Path(str(bam) + ".bai"),
    }
    missing = [f"{key}: {path}" for key, path in required.items() if not path.is_file()]
    if missing:
        raise RuntimeError(
            "Genome browser requires the Stage 15 target BAM and indexed GRCh37 FASTA. Missing: "
            + "; ".join(missing)
        )
    artifacts.update(required)

    if calls.is_file() and Path(str(calls) + ".tbi").is_file():
        artifacts["calls_vcf"] = calls
        artifacts["calls_tbi"] = Path(str(calls) + ".tbi")
    if clinvar.is_file() and Path(str(clinvar) + ".tbi").is_file():
        artifacts["clinvar_vcf"] = clinvar
        artifacts["clinvar_tbi"] = Path(str(clinvar) + ".tbi")
    if targets.is_file():
        artifacts["targets_tsv"] = targets
    return artifacts

def choose_initial_locus(artifacts: dict[str, Path]) -> str:
    calls = artifacts.get("calls_vcf")
    if calls:
        try:
            with gzip.open(calls, "rt", encoding="utf-8", errors="replace") as f:
                for line in f:
                    if not line.startswith("#"):
                        parts = line.rstrip("\n").split("\t")
                        chrom, pos = parts[0], int(parts[1])
                        start = max(1, pos - 80)
                        end = pos + 80
                        return f"{chrom}:{start}-{end}"
        except Exception:
            pass
    return "1:850000-870000"

class RsidIndex:
    def __init__(self, patient_root: Path, target_tsv: Path | None):
        self.root = patient_root
        self.target_tsv = target_tsv
        self.path = patient_root / "Dante_work" / "genome_browser_rsid.sqlite3"
        self.status = "UNAVAILABLE" if target_tsv is None else "PENDING"
        self.rows = 0
        self.error = ""
        self._lock = threading.Lock()

    def _source_signature(self) -> str:
        if not self.target_tsv:
            return ""
        st = self.target_tsv.stat()
        return f"{st.st_size}:{st.st_mtime_ns}"

    def is_current(self) -> bool:
        if not self.path.is_file() or not self.target_tsv:
            return False
        try:
            con = sqlite3.connect(self.path)
            meta = dict(con.execute("SELECT key,value FROM metadata"))
            count = int(con.execute("SELECT COUNT(*) FROM loci").fetchone()[0])
            con.close()
            if meta.get("schema") == INDEX_SCHEMA and meta.get("source_signature") == self._source_signature():
                self.rows = count
                return True
        except Exception:
            return False
        return False

    def start(self) -> None:
        if self.target_tsv is None:
            return
        if self.is_current():
            self.status = "READY"
            return
        self.status = "BUILDING"
        threading.Thread(target=self._build, name="rsid-index", daemon=True).start()

    def _build(self) -> None:
        assert self.target_tsv is not None
        tmp = Path(str(self.path) + ".part")
        try:
            try:
                tmp.unlink()
            except FileNotFoundError:
                pass
            con = sqlite3.connect(tmp)
            con.execute("PRAGMA journal_mode=OFF")
            con.execute("PRAGMA synchronous=OFF")
            con.execute("PRAGMA temp_store=MEMORY")
            con.execute(
                "CREATE TABLE loci(rsid TEXT PRIMARY KEY, chromosome TEXT NOT NULL, position INTEGER NOT NULL, ref TEXT, alts TEXT)"
            )
            con.execute("CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL)")
            rows = 0
            batch: list[tuple[str,str,int,str,str]] = []
            with self.target_tsv.open("rt", encoding="utf-8-sig", errors="replace", newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                fields = {str(x).lower(): x for x in (reader.fieldnames or [])}
                rs_col = fields.get("rsid")
                chr_col = fields.get("chromosome") or fields.get("chrom") or fields.get("chr")
                pos_col = fields.get("position") or fields.get("pos")
                ref_col = fields.get("ref")
                alt_col = fields.get("alts") or fields.get("alt")
                if not rs_col or not chr_col or not pos_col:
                    raise RuntimeError(f"Unexpected target panel header: {reader.fieldnames}")
                for row in reader:
                    rsid = (row.get(rs_col) or "").strip().lower()
                    if not re.fullmatch(r"rs\d+", rsid):
                        continue
                    chrom = (row.get(chr_col) or "").strip()
                    try:
                        pos = int((row.get(pos_col) or "").strip())
                    except ValueError:
                        continue
                    ref = (row.get(ref_col) or "").strip() if ref_col else ""
                    alts = (row.get(alt_col) or "").strip() if alt_col else ""
                    batch.append((rsid, chrom, pos, ref, alts))
                    if len(batch) >= 10000:
                        con.executemany(
                            "INSERT OR REPLACE INTO loci(rsid,chromosome,position,ref,alts) VALUES(?,?,?,?,?)",
                            batch,
                        )
                        rows += len(batch)
                        batch.clear()
                        with self._lock:
                            self.rows = rows
                if batch:
                    con.executemany(
                        "INSERT OR REPLACE INTO loci(rsid,chromosome,position,ref,alts) VALUES(?,?,?,?,?)",
                        batch,
                    )
                    rows += len(batch)
            con.executemany(
                "INSERT INTO metadata(key,value) VALUES(?,?)",
                [
                    ("schema", INDEX_SCHEMA),
                    ("source_signature", self._source_signature()),
                    ("source_path", str(self.target_tsv)),
                    ("rows", str(rows)),
                ],
            )
            con.commit()
            con.close()
            os.replace(tmp, self.path)
            with self._lock:
                self.rows = rows
                self.status = "READY"
        except Exception as exc:
            try:
                tmp.unlink()
            except FileNotFoundError:
                pass
            with self._lock:
                self.error = f"{type(exc).__name__}: {exc}"
                self.status = "FAILED"

    def lookup(self, rsid: str) -> dict[str, Any] | None:
        if self.status != "READY":
            return None
        con = sqlite3.connect(self.path)
        row = con.execute(
            "SELECT rsid,chromosome,position,ref,alts FROM loci WHERE rsid=?",
            (rsid.lower(),),
        ).fetchone()
        con.close()
        if not row:
            return None
        return {
            "rsid": row[0],
            "chromosome": row[1],
            "position": int(row[2]),
            "ref": row[3] or "",
            "alts": row[4] or "",
        }

def create_app(root: Path, igv_path: Path, artifacts: dict[str, Path], index: RsidIndex):
    try:
        from flask import Flask, abort, jsonify, request, send_file
    except ImportError as exc:
        raise RuntimeError(
            "Flask is not installed. Launch item 24 from start.cmd; the pipeline will install the pinned browser dependency automatically."
        ) from exc

    app = Flask(__name__)
    app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0
    allowlist = dict(artifacts)
    server_holder: dict[str, Any] = {}

    @app.after_request
    def security_headers(response):
        response.headers["Cache-Control"] = "no-store"
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["X-Frame-Options"] = "SAMEORIGIN"
        return response

    @app.get("/")
    def home():
        return HTML

    @app.get("/vendor/igv.min.js")
    def vendor_igv():
        return send_file(igv_path, mimetype="application/javascript", conditional=True, etag=True)

    @app.get("/data/<key>")
    def data(key: str):
        path = allowlist.get(key)
        if path is None:
            abort(404)
        # Flask/Werkzeug conditional file responses support byte-range requests needed by BAM/VCF/FASTA indexes.
        return send_file(path, conditional=True, etag=True, max_age=0)

    @app.get("/api/config")
    def api_config():
        tracks: list[dict[str, Any]] = [
            {
                "name": "Dante target BAM (Stage 15)",
                "type": "alignment",
                "format": "bam",
                "url": "/data/bam",
                "indexURL": "/data/bai",
                "visibilityWindow": 100000,
                "height": 420,
                "showCoverage": True,
                "displayMode": "EXPANDED",
            }
        ]
        names = ["target BAM"]
        if "calls_vcf" in artifacts:
            tracks.append(
                {
                    "name": "Dante Stage 16 genotype calls",
                    "type": "variant",
                    "format": "vcf",
                    "url": "/data/calls_vcf",
                    "indexURL": "/data/calls_tbi",
                    "visibilityWindow": 1000000,
                }
            )
            names.append("Stage 16 VCF")
        if "clinvar_vcf" in artifacts:
            tracks.append(
                {
                    "name": "ClinVar GRCh37",
                    "type": "variant",
                    "format": "vcf",
                    "url": "/data/clinvar_vcf",
                    "indexURL": "/data/clinvar_tbi",
                    "visibilityWindow": 1000000,
                }
            )
            names.append("ClinVar")
        cfg = {
            "reference": {
                "id": "GRCh37_PIPELINE",
                "name": "GRCh37 — exact pipeline FASTA",
                "fastaURL": "/data/reference_fasta",
                "indexURL": "/data/reference_fai",
            },
            "locus": choose_initial_locus(artifacts),
            "tracks": tracks,
            "showNavigation": True,
            "showRuler": True,
            "showCenterGuide": True,
            "showCursorTrackingGuide": True,
        }
        return jsonify(
            {
                "igv": cfg,
                "reference_name": "GRCh37 exact pipeline FASTA",
                "track_names": names,
                "index_status": index.status,
            }
        )

    @app.get("/api/status")
    def api_status():
        return jsonify(
            {
                "index_status": index.status,
                "index_rows": index.rows,
                "index_error": index.error,
                "localhost_only": True,
            }
        )

    @app.get("/api/search")
    def api_search():
        q = (request.args.get("q") or "").strip().lower()
        if not re.fullmatch(r"rs\d+", q):
            return jsonify({"error": "Only rsID lookup is handled by the local index. Coordinates can be entered directly in IGV."}), 400
        if index.status == "BUILDING":
            return jsonify({"error": f"rsID index is still building ({index.rows:,} rows indexed). Try again shortly."}), 503
        if index.status != "READY":
            return jsonify({"error": f"rsID index is not ready: {index.status}. {index.error}"}), 503
        hit = index.lookup(q)
        if not hit:
            return jsonify({"error": f"{q} was not found in the reference target panel."}), 404
        pos = int(hit["position"])
        hit["locus"] = f"{hit['chromosome']}:{max(1,pos-80)}-{pos+80}"
        return jsonify(hit)

    @app.post("/api/shutdown")
    def shutdown():
        server = server_holder.get("server")
        if server is not None:
            threading.Thread(target=server.shutdown, daemon=True).start()
        return jsonify({"status": "stopping"})

    app.server_holder = server_holder  # type: ignore[attr-defined]
    return app

def run_server(root: Path, patient_root: Path, port: int, open_browser: bool) -> int:
    artifacts = resolve_artifacts(root, patient_root)
    igv_path = ensure_igv(root)
    index = RsidIndex(patient_root, artifacts.get("targets_tsv"))
    index.start()

    try:
        from werkzeug.serving import make_server
    except ImportError as exc:
        raise RuntimeError("Werkzeug/Flask runtime is missing.") from exc

    app = create_app(root, igv_path, artifacts, index)

    last_error: Exception | None = None
    server = None
    selected_port = port
    for candidate in range(port, port + 20):
        try:
            server = make_server("127.0.0.1", candidate, app, threaded=True)
            selected_port = candidate
            break
        except OSError as exc:
            last_error = exc
    if server is None:
        raise RuntimeError(f"Could not bind localhost ports {port}-{port+19}: {last_error}")

    app.server_holder["server"] = server  # type: ignore[attr-defined]
    url = f"http://127.0.0.1:{selected_port}/"
    print("")
    print("DANTE LOCAL GENOME BROWSER")
    print(f"  URL:           {url}")
    print("  Bind address:  127.0.0.1 only")
    print(f"  IGV.js:        {IGV_VERSION} (cached locally)")
    print(f"  Patient root:  {patient_root}")
    print(f"  BAM:           {artifacts['bam']}")
    print("  BAM scope:     TARGET-ONLY — empty inter-target areas are not evidence of no WGS coverage")
    print("  Stop:          click 'Stop server' in the browser or press Ctrl+C")
    print("")
    if open_browser:
        threading.Timer(0.7, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    print("[BROWSER] Local server stopped.")
    return 0

def self_test() -> int:
    with tempfile_dir() as root:
        qc = root / "Dante_qc"
        work = root / "Dante_work"
        fasta_dir = root / "annotation_data" / "vep_rsid" / "fasta"
        qc.mkdir(parents=True)
        work.mkdir(parents=True)
        fasta_dir.mkdir(parents=True)
        (work / BAM_NAME).write_bytes(b"BAM")
        (work / (BAM_NAME + ".bai")).write_bytes(b"BAI")
        fasta = fasta_dir / FASTA_NAME
        fasta.write_text(">1\nACGT\n", encoding="ascii")
        Path(str(fasta) + ".fai").write_text("1\t4\t3\t4\t5\n", encoding="ascii")
        target = qc / TARGET_TSV_NAME
        target.write_text(
            "rsid\tchromosome\tposition\tref\talts\n"
            "rs1\t1\t100\tA\tG\n"
            "rs2\tX\t200\tC\tT\n",
            encoding="utf-8",
        )
        patient = root / "patients" / "Demo 5"
        # Move synthetic patient-specific directories under the patient workspace.
        patient.mkdir(parents=True, exist_ok=True)
        os.replace(work, patient / "Dante_work")
        os.replace(qc, patient / "Dante_qc")
        artifacts = resolve_artifacts(root, patient)
        assert artifacts["bam"].name == BAM_NAME
        idx = RsidIndex(patient, patient / "Dante_qc" / TARGET_TSV_NAME)
        idx._build()
        assert idx.status == "READY"
        hit = idx.lookup("rs2")
        assert hit and hit["chromosome"] == "X" and hit["position"] == 200
    print("SELF-TEST PASS | artifact allowlist + rsID SQLite index + target-only browser contract")
    return 0

class tempfile_dir:
    def __enter__(self):
        import tempfile
        self._ctx = tempfile.TemporaryDirectory()
        return Path(self._ctx.__enter__())
    def __exit__(self, *args):
        return self._ctx.__exit__(*args)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path.cwd())
    ap.add_argument("--patient-root", type=Path)
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--open-browser", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    root = args.root.resolve()
    patient_root = args.patient_root.resolve() if args.patient_root else root
    return run_server(root, patient_root, args.port, args.open_browser)

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
