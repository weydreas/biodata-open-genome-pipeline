"""Runtime delta for Dante WGS v9.0.0 hotfix6 / knowledge rules-v8.

This module patches the frozen v7 base module. It exists so GitHub can deploy a
small reviewed delta while preserving the already-seeded full v7 source blob.
The next materialized release may fold this delta back into the main source.
"""
from __future__ import annotations

import bisect
import gzip
import re
from collections import Counter, defaultdict
from typing import Any

RULES_SCHEMA = "dante-knowledge-rules-v8"


def install(m):
    # Save v7 implementations before patching.
    original_member_context = m.structural_member_context
    original_scan_structural = m.scan_structural
    original_apply_alphamissense = m.apply_alphamissense
    original_patient_priority = m.patient_priority
    original_write_usage_manifests = m.write_usage_manifests

    def chrom_sort_key(value: str) -> tuple[int, str]:
        chrom = m.normalize_chrom(value)
        if chrom.isdigit():
            return int(chrom), chrom
        return {"X": 23, "Y": 24, "MT": 25}.get(chrom, 99), chrom

    def parse_bnd_partner(alt: str) -> tuple[str, int] | None:
        match = re.search(r"[\[\]](?:chr)?([^:\[\]]+):(\d+)[\[\]]", alt or "", flags=re.IGNORECASE)
        if not match:
            return None
        chrom = m.normalize_chrom(match.group(1))
        if not m.canonical_chrom(chrom):
            return None
        try:
            pos = int(match.group(2))
        except ValueError:
            return None
        return chrom, pos

    def collapse_breakpoint_intervals(breakpoints, tolerance: int = 1000):
        by_chrom = defaultdict(list)
        for bp in breakpoints:
            chrom = m.normalize_chrom(str(bp.get("chrom", "")))
            if not m.canonical_chrom(chrom):
                continue
            try:
                start = int(bp.get("start", 0) or 0)
                end = int(bp.get("end", 0) or start)
            except Exception:
                continue
            if start <= 0:
                continue
            if end < start:
                start, end = end, start
            by_chrom[chrom].append((start, end))
        merged = []
        for chrom, intervals in by_chrom.items():
            current = None
            for start, end in sorted(intervals):
                if current is not None and start <= current[1] + tolerance:
                    current[1] = max(current[1], end)
                else:
                    if current is not None:
                        merged.append({"chrom": chrom, "start": current[0], "end": current[1]})
                    current = [start, end]
            if current is not None:
                merged.append({"chrom": chrom, "start": current[0], "end": current[1]})
        return sorted(merged, key=lambda bp: (chrom_sort_key(str(bp["chrom"])), int(bp["start"]), int(bp["end"])))

    def structural_member_context(row: dict[str, str]) -> dict[str, Any]:
        payload = original_member_context(row)
        chrom = m.normalize_chrom(m.first_value(row, ["SV_chrom", "chrom", "chr"], ["svchrom"]))
        try:
            start = int(float(m.first_value(row, ["SV_start", "start"], ["svstart"]) or 0))
            end = int(float(m.first_value(row, ["SV_end", "end"], ["svend"]) or start))
        except ValueError:
            start = end = 0
        alt = m.first_value(row, ["ALT", "alt"], ["alt"]) or "."
        info = m.parse_info_kv(m.first_value(row, ["INFO", "info"], ["info"]))
        partner = parse_bnd_partner(alt)
        payload.update({
            "chrom": chrom,
            "start": start,
            "end": end,
            "alt": alt,
            "mate_id": info.get("MATEID", ""),
            "partner_chrom": partner[0] if partner else ".",
            "partner_position": partner[1] if partner else 0,
        })
        return payload

    def _member_breakpoints(event):
        raw = []
        members = event.raw_context.get("structural_members", [])
        if isinstance(members, list):
            for member in members:
                if not isinstance(member, dict):
                    continue
                chrom = m.normalize_chrom(str(member.get("chrom", "")))
                try:
                    start = int(member.get("start", 0) or 0)
                    end = int(member.get("end", 0) or start)
                except Exception:
                    continue
                if m.canonical_chrom(chrom) and start > 0:
                    raw.append({"chrom": chrom, "start": start, "end": max(start, end)})
                mate_chrom = m.normalize_chrom(str(member.get("partner_chrom", "")))
                try:
                    mate_pos = int(member.get("partner_position", 0) or 0)
                except Exception:
                    mate_pos = 0
                if m.canonical_chrom(mate_chrom) and mate_pos > 0:
                    raw.append({"chrom": mate_chrom, "start": mate_pos, "end": mate_pos})
        return raw

    def scan_structural(path, kind: str, knowledge):
        result, counts = original_scan_structural(path, kind, knowledge)
        rewritten = {}
        for old_key, ev in result.items():
            sv_type = str(ev.raw_context.get("sv_type", "")).upper()
            raw = _member_breakpoints(ev)
            if sv_type in {"TRA", "BND"}:
                loci = collapse_breakpoint_intervals(raw)
                if loci:
                    first = loci[0]
                    ev.chrom = str(first["chrom"])
                    ev.position = int(first["start"])
                    ev.end = int(first["end"])
                    second = loci[1] if len(loci) > 1 else None
                    cluster_id = str(ev.raw_context.get("structural_cluster_id", old_key))
                    mate = f":{second['chrom']}:{second['start']}-{second['end']}" if second else ""
                    ev.event_key = (
                        f"GRCh37:{kind.upper()}:{ev.chrom}:{ev.position}-{ev.end}"
                        f"{mate}:{sv_type}:{m.deterministic_id('SC', cluster_id)}"
                    )
                    ev.raw_context["breakpoints"] = loci
                    if second:
                        ev.raw_context["mate_chrom"] = second["chrom"]
                        ev.raw_context["mate_start"] = second["start"]
                        ev.raw_context["mate_end"] = second["end"]
                    ev.raw_context["sv_length_abs"] = "."
            else:
                # Keep ordinary interval identity unchanged, but preserve real
                # per-caller breakpoints for GIAB QC and cross-caller merges.
                loci = collapse_breakpoint_intervals(raw, tolerance=0)
                if loci:
                    ev.raw_context["breakpoints"] = loci
            if ev.event_key in rewritten:
                # This should only happen if two independently keyed calls collapse
                # to the exact same deterministic identity, which is a data bug.
                raise RuntimeError(f"duplicate structural identity after breakpoint repair: {ev.event_key}")
            rewritten[ev.event_key] = ev
        counts["breakpoint_safe_events"] = len(rewritten)
        return rewritten, counts

    def structural_length(event) -> int:
        if str(event.raw_context.get("sv_type", "")).upper() in {"TRA", "BND"}:
            return 0
        return max(0, int(event.end) - int(event.position))

    def structural_breakpoint_positions(event) -> list[tuple[str, int]]:
        points = set()
        raw = _member_breakpoints(event)
        if not raw:
            raw = event.raw_context.get("breakpoints", []) if isinstance(event.raw_context.get("breakpoints", []), list) else []
        for bp in raw:
            chrom = m.normalize_chrom(str(bp.get("chrom", "")))
            if not m.canonical_chrom(chrom):
                continue
            for field in ("start", "end"):
                try:
                    pos = int(bp.get(field, 0) or 0)
                except Exception:
                    pos = 0
                if pos > 0:
                    points.add((chrom, pos))
        if not points:
            points.add((event.chrom, int(event.position)))
            if int(event.end) > 0:
                points.add((event.chrom, int(event.end)))
        return sorted(points, key=lambda x: (chrom_sort_key(x[0]), x[1]))

    def apply_alphamissense(events, research_root, research_manifest, image: str, work):
        counts = original_apply_alphamissense(events, research_root, research_manifest, image, work)
        promoted = 0
        for ev in events.values():
            if str(getattr(ev, "alpha_score", ".")) not in {"", "."}:
                if "K8_ALPHAMISSENSE_SUPPORT" not in ev.research_reasons:
                    promoted += 1
                ev.research_reasons.add("K8_ALPHAMISSENSE_SUPPORT")
        counts["research_overlays_materialized"] = promoted
        return counts

    def apply_giab(events, beds):
        counts = Counter()
        positions = defaultdict(list)
        pos_events = defaultdict(list)
        for ev in events.values():
            event_positions = (
                structural_breakpoint_positions(ev)
                if ev.raw_context.get("structural_members")
                else [(ev.chrom, ev.position)]
            )
            for chrom, pos in event_positions:
                positions[chrom].append(pos)
                pos_events[(chrom, pos)].append(ev)
        for chrom in positions:
            positions[chrom] = sorted(set(positions[chrom]))
        for label, path in beds:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                for raw_line in handle:
                    if raw_line.startswith("#") or not raw_line.strip():
                        continue
                    fields = raw_line.rstrip("\r\n").split("\t")
                    if len(fields) < 3:
                        continue
                    chrom = m.normalize_chrom(fields[0])
                    vals = positions.get(chrom)
                    if not vals:
                        continue
                    try:
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

    def patient_priority(reason: str, event, gene: str, profile, inheritance: str, evidence_rank: int):
        rank, relevance = original_patient_priority(reason, event, gene, profile, inheritance, evidence_rank)
        flags = [] if relevance == "DIRECT_EVIDENCE_CONTEXT" else [x for x in relevance.split(";") if x]
        if reason.startswith(("K5_", "K6_")):
            tech = m.event_technical_confidence(event)
            base_tech = tech.split("+", 1)[0]
            if base_tech != "CROSS_CALLER_CONFIRMED":
                rank = max(rank, 5)
                if "TECHNICAL_SINGLE_CALLER" not in flags:
                    flags.append("TECHNICAL_SINGLE_CALLER")
            elif "GIAB_CAUTION" in tech:
                rank = max(rank, 4)
            if base_tech.startswith("CAUTION") and "TECHNICAL_CALL_CAUTION" not in flags:
                flags.append("TECHNICAL_CALL_CAUTION")
        return rank, (";".join(flags) if flags else "DIRECT_EVIDENCE_CONTEXT")

    def write_usage_manifests(output_dir, open_rows, research_rows):
        open_path, research_path = original_write_usage_manifests(output_dir, open_rows, research_rows)
        payload = m.load_json(research_path)
        payload["note"] = (
            "Research-only publication manifest. PanelApp/AlphaMissense/CADD restrictions apply "
            "where those sources are present."
        )
        m.save_json_atomic(research_path, payload)
        return open_path, research_path

    def self_test() -> int:
        assert m.normalize_chrom("chrX") == "X"
        assert parse_bnd_partner("C]14:93713179]") == ("14", 93713179)
        assert collapse_breakpoint_intervals([
            {"chrom": "1", "start": 100, "end": 101},
            {"chrom": "14", "start": 200, "end": 201},
        ]) == [
            {"chrom": "1", "start": 100, "end": 101},
            {"chrom": "14", "start": 200, "end": 201},
        ]
        print("SELF-TEST PASS | v9.0.0 hotfix6 knowledge rules-v8 + AlphaMissense overlays + breakpoint-safe BND/TRA")
        return 0

    # Install patched globals in the v7 base module. Functions defined in the
    # base resolve these names at runtime, so the delta is used everywhere.
    m.RULES_SCHEMA = RULES_SCHEMA
    m.chrom_sort_key = chrom_sort_key
    m.parse_bnd_partner = parse_bnd_partner
    m.collapse_breakpoint_intervals = collapse_breakpoint_intervals
    m.structural_member_context = structural_member_context
    m.scan_structural = scan_structural
    m.structural_length = structural_length
    m.structural_breakpoint_positions = structural_breakpoint_positions
    m.apply_alphamissense = apply_alphamissense
    m.apply_giab = apply_giab
    m.patient_priority = patient_priority
    m.write_usage_manifests = write_usage_manifests
    m.self_test = self_test
    return m
