#!/usr/bin/env python3
"""Parse pinned bcftools norm summaries for genotyping-array QC.

The array runtime uses ``bcftools norm --check-ref x``. bcftools summary fields
have changed across releases, so this helper parses the labels and values rather
than assuming one fixed width. For bcftools 1.24 the explicit
``mismatch_removed`` counter is the authoritative REF/reference-mismatch
exclusion count. Older supported summaries expose only ``skipped``; under the
same pinned ``--check-ref x`` contract that legacy counter is retained as the
REF-mismatch exclusion count.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

SUMMARY_LINE_RE = re.compile(
    r"^\s*Lines\s+([^:\r\n]+):\s*([0-9]+(?:/[0-9]+)+)\s*$",
    re.IGNORECASE | re.MULTILINE,
)
LEGACY_LABELS = {
    ("total", "split", "realigned", "skipped"),
    ("total", "split", "joined", "realigned", "skipped"),
}


def parse_norm_summary(text: str) -> dict[str, int | str]:
    match = SUMMARY_LINE_RE.search(text)
    if not match:
        raise ValueError("bcftools norm summary line was not found")

    labels = [x.strip().lower() for x in match.group(1).split("/")]
    values = [int(x) for x in match.group(2).split("/")]
    if len(labels) != len(values):
        raise ValueError(
            f"bcftools norm summary label/value width mismatch: {len(labels)} labels vs {len(values)} values"
        )

    counters = dict(zip(labels, values))
    required = {"total", "split", "realigned"}
    missing = sorted(required - counters.keys())
    if missing:
        raise ValueError(f"bcftools norm summary missing required counters: {', '.join(missing)}")

    joined = counters.get("joined", 0)
    duplicate_removed = counters.get("dup_removed", 0)
    skipped = counters.get("skipped", 0)

    if "mismatch_removed" in counters:
        mismatch_removed = counters["mismatch_removed"]
        mismatch_source = "mismatch_removed"
        format_name = f"labeled-{len(values)}-field-mismatch_removed"
    elif tuple(labels) in LEGACY_LABELS:
        mismatch_removed = counters["skipped"]
        mismatch_source = "legacy_skipped"
        format_name = f"legacy-{len(values)}-field-skipped"
    else:
        raise ValueError(
            "bcftools norm summary has no supported REF-mismatch counter "
            f"(labels: {'/'.join(labels)})"
        )

    return {
        "schema": "dante-array-bcftools-norm-summary-v2",
        "summary_format": format_name,
        "input_records": counters["total"],
        "split_records": counters["split"],
        "joined_records": joined,
        "realigned_records": counters["realigned"],
        "reference_mismatch_records_excluded": mismatch_removed,
        "reference_mismatch_counter_source": mismatch_source,
        "duplicate_records_removed": duplicate_removed,
        "skipped_records": skipped,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)

    parse = sub.add_parser("parse-log")
    parse.add_argument("log", type=Path)
    parse.add_argument("--json", action="store_true")

    sub.add_parser("self-test")
    args = ap.parse_args()

    if args.mode == "parse-log":
        result = parse_norm_summary(args.log.read_text(encoding="utf-8", errors="replace"))
        if args.json:
            print(json.dumps(result, sort_keys=True, separators=(",", ":")))
        else:
            for key, value in result.items():
                print(f"{key}: {value}")
        return 0

    old = parse_norm_summary("Lines   total/split/realigned/skipped:\t145626/0/0/15\n")
    assert old["input_records"] == 145626
    assert old["joined_records"] == 0
    assert old["reference_mismatch_records_excluded"] == 15
    assert old["reference_mismatch_counter_source"] == "legacy_skipped"

    five = parse_norm_summary(
        "some warning\nLines   total/split/joined/realigned/skipped:\t145626/3/1/2/15\n"
    )
    assert five["input_records"] == 145626
    assert five["split_records"] == 3
    assert five["joined_records"] == 1
    assert five["realigned_records"] == 2
    assert five["reference_mismatch_records_excluded"] == 15
    assert five["reference_mismatch_counter_source"] == "legacy_skipped"

    production_124 = parse_norm_summary(
        "Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:      "
        "145626/0/0/0/15/0/0\n"
    )
    assert production_124["input_records"] == 145626
    assert production_124["reference_mismatch_records_excluded"] == 15
    assert production_124["reference_mismatch_counter_source"] == "mismatch_removed"
    assert production_124["duplicate_records_removed"] == 0
    assert production_124["skipped_records"] == 0
    assert production_124["summary_format"] == "labeled-7-field-mismatch_removed"

    try:
        parse_norm_summary("no summary here")
    except ValueError:
        pass
    else:
        raise AssertionError("missing summary must fail closed")

    try:
        parse_norm_summary("Lines total/split/realigned/dup_removed/skipped: 10/0/0/1/0\n")
    except ValueError:
        pass
    else:
        raise AssertionError("unsupported modern summary without mismatch_removed must fail closed")

    print("SELF-TEST PASS | bcftools norm labeled summaries + explicit REF-mismatch counter contract")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
