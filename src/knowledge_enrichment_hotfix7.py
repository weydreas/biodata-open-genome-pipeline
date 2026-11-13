"""Input-scope hardening over knowledge rules-v8.

For genotyping-array patients the scientific finding rules stay unchanged, but
technical confidence must not look like read-backed WGS evidence. The array
adapter has no DP/GQ/QUAL/FILTER read-level contract, so every small-variant
candidate is explicitly labelled accordingly.
"""
from __future__ import annotations

import os

import knowledge_enrichment_hotfix6 as h6

RULES_SCHEMA = h6.RULES_SCHEMA
ARRAY_KIND = "GENOTYPING_ARRAY_VCF"
ARRAY_TECHNICAL_CONFIDENCE = "ARRAY_VENDOR_CALL_NO_READ_LEVEL_QC"


def install(m):
    h6.install(m)
    original_event_technical_confidence = m.event_technical_confidence
    original_self_test = m.self_test

    def event_technical_confidence(event):
        if (
            os.environ.get("DANTE_INPUT_PROFILE_KIND", "").strip().upper() == ARRAY_KIND
            and getattr(event, "event_type", "") == "SNV_INDEL"
        ):
            suffix = "+GIAB_CAUTION" if getattr(event, "giab_flags", set()) else ""
            return ARRAY_TECHNICAL_CONFIDENCE + suffix
        return original_event_technical_confidence(event)

    def self_test() -> int:
        original = os.environ.get("DANTE_INPUT_PROFILE_KIND")
        try:
            os.environ["DANTE_INPUT_PROFILE_KIND"] = ARRAY_KIND
            event = type("E", (), {"event_type": "SNV_INDEL", "giab_flags": set(), "raw_context": {}})()
            assert event_technical_confidence(event) == ARRAY_TECHNICAL_CONFIDENCE
            event.giab_flags = {"difficult"}
            assert event_technical_confidence(event) == ARRAY_TECHNICAL_CONFIDENCE + "+GIAB_CAUTION"
        finally:
            if original is None:
                os.environ.pop("DANTE_INPUT_PROFILE_KIND", None)
            else:
                os.environ["DANTE_INPUT_PROFILE_KIND"] = original
        original_self_test()
        print("SELF-TEST PASS | input-scope technical confidence")
        return 0

    m.event_technical_confidence = event_technical_confidence
    m.self_test = self_test
    return m
