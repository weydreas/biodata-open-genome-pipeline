"""Genotype-table input-scope hardening over knowledge hotfix7.

WGS and array-VCF behavior is delegated unchanged.  For the independent
GENOTYPE_TABLE_TXT profile, small-variant findings are explicitly marked as
source genotype calls without read-level sequencing QC.
"""
from __future__ import annotations

import os

import knowledge_enrichment_hotfix7 as h7

RULES_SCHEMA = h7.RULES_SCHEMA
TXT_KIND = "GENOTYPE_TABLE_TXT"
TXT_TECHNICAL_CONFIDENCE = "GENOTYPE_TABLE_CALL_NO_READ_LEVEL_QC"


def install(m):
    h7.install(m)
    original_event_technical_confidence = m.event_technical_confidence
    original_self_test = m.self_test

    def event_technical_confidence(event):
        if (
            os.environ.get("DANTE_INPUT_PROFILE_KIND", "").strip().upper() == TXT_KIND
            and getattr(event, "event_type", "") == "SNV_INDEL"
        ):
            suffix = "+GIAB_CAUTION" if getattr(event, "giab_flags", set()) else ""
            return TXT_TECHNICAL_CONFIDENCE + suffix
        return original_event_technical_confidence(event)

    def self_test() -> int:
        original = os.environ.get("DANTE_INPUT_PROFILE_KIND")
        try:
            os.environ["DANTE_INPUT_PROFILE_KIND"] = TXT_KIND
            event = type("E", (), {"event_type": "SNV_INDEL", "giab_flags": set(), "raw_context": {}})()
            assert event_technical_confidence(event) == TXT_TECHNICAL_CONFIDENCE
            event.giab_flags = {"difficult"}
            assert event_technical_confidence(event) == TXT_TECHNICAL_CONFIDENCE + "+GIAB_CAUTION"
        finally:
            if original is None:
                os.environ.pop("DANTE_INPUT_PROFILE_KIND", None)
            else:
                os.environ["DANTE_INPUT_PROFILE_KIND"] = original
        original_self_test()
        print("SELF-TEST PASS | genotype-table TXT technical confidence + prior array/WGS scope")
        return 0

    m.event_technical_confidence = event_technical_confidence
    m.self_test = self_test
    return m
