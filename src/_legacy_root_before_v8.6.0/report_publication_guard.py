#!/usr/bin/env python3
"""Future PDF publication gate based on sources actually used in a report."""
from __future__ import annotations
import argparse, json, tempfile
from pathlib import Path
SCHEMA = "knowledge-usage-v1"

def check(report_type: str, usage: dict):
    if usage.get("schema") != SCHEMA:
        raise RuntimeError(f"Unsupported usage manifest schema: {usage.get('schema')}")
    blocked, notices = [], []
    for item in usage.get("sources_used", []):
        source = str(item.get("source","UNKNOWN"))
        scope = str(item.get("publication_scope","UNKNOWN")).upper()
        prohibited = [str(x).lower() for x in item.get("prohibited_uses",[])]
        if report_type == "COMMERCIAL" and scope != "COMMERCIAL_ALLOWED":
            blocked.append(f"{source}: publication_scope={scope}")
        if report_type == "RESEARCH" and scope == "UNKNOWN":
            blocked.append(f"{source}: unknown publication scope")
        if any(x in prohibited for x in (
            "diagnostic use","medical decision-making","providing a healthcare service"
        )):
            notices.append(
                f"{source}: research/hypothesis-use restriction applies; no diagnostic, "
                "medical-decision, or healthcare-service use without separate permission."
            )
    return not blocked, blocked, notices

def self_test() -> int:
    usage = {"schema":SCHEMA,"sources_used":[
        {"source":"HGNC","publication_scope":"COMMERCIAL_ALLOWED","prohibited_uses":[]},
        {"source":"PanelApp","publication_scope":"RESEARCH_ONLY",
         "prohibited_uses":["commercial","diagnostic use","medical decision-making","providing a healthcare service"]}
    ]}
    ok_c, blocked_c, _ = check("COMMERCIAL", usage)
    ok_r, blocked_r, notices_r = check("RESEARCH", usage)
    assert not ok_c and any("PanelApp" in x for x in blocked_c)
    assert ok_r and not blocked_r and notices_r
    print("SELF-TEST PASS | commercial blocks research-only; research permits with restriction notices")
    return 0

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("check","self-test"))
    ap.add_argument("--report-type", choices=("COMMERCIAL","RESEARCH"))
    ap.add_argument("--usage-manifest", type=Path)
    a = ap.parse_args()
    if a.mode=="self-test": return self_test()
    if not a.report_type or not a.usage_manifest:
        raise SystemExit("check requires --report-type and --usage-manifest")
    usage = json.loads(a.usage_manifest.read_text(encoding="utf-8"))
    ok, blocked, notices = check(a.report_type, usage)
    for n in notices: print("NOTICE:", n)
    if not ok:
        for b in blocked: print("BLOCK:", b)
        print(f"PUBLICATION BLOCKED for {a.report_type}")
        return 2
    print(f"PUBLICATION ALLOWED for {a.report_type}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
