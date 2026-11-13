#!/usr/bin/env python3
"""Versioned-output publication layer for deterministic Stage-24 genome PDFs.

Builds on hotfix9 scientific/publication semantics. The only new behavior is output
allocation: an existing report package is never overwritten. If any member of v1
exists, the next run uses _v2; if any member of _v2 exists, the next run uses _v3,
and so on. Main PDF, common-major appendix and manifest always share one version.
"""
from __future__ import annotations

import json
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from genome_report_hotfix9 import *  # noqa: F401,F403
import genome_report_hotfix9 as h9

PUBLICATION_CONTRACT = "dante-report-publication-v10"
VERSIONING_RULE = "never overwrite; allocate next free package suffix _vN (v1 is unsuffixed)"


def _version_suffix(version: int) -> str:
    if version < 1:
        raise ValueError("report version must be >= 1")
    return "" if version == 1 else f"_v{version}"


def report_package_paths(out_root: Path, base_name: str, version: int) -> tuple[Path, Path, Path]:
    suffix = _version_suffix(version)
    return (
        out_root / f"{base_name}{suffix}.pdf",
        out_root / f"{base_name}_COMMON_MAJOR_REFERENCE{suffix}.pdf",
        out_root / f"{base_name}{suffix}.manifest.json",
    )


def allocate_report_package(out_root: Path, base_name: str) -> tuple[int, Path, Path, Path]:
    """Return the first completely unused report package version.

    A partial prior run also consumes its version: if even one expected file exists,
    the next run advances to the following version. This intentionally avoids both
    overwriting and Windows file-lock failures from open PDFs.
    """
    out_root.mkdir(parents=True, exist_ok=True)
    version = 1
    while True:
        main_pdf, reference_pdf, manifest_path = report_package_paths(out_root, base_name, version)
        if not any(path.exists() for path in (main_pdf, reference_pdf, manifest_path)):
            return version, main_pdf, reference_pdf, manifest_path
        version += 1


def _write_versioned_manifest(
    *,
    output_path: Path,
    version: int,
    mode: str,
    patient: str,
    state: dict,
    root: Path,
    source_files: dict,
    research_layer_opened: bool,
    main_grouped: dict[str, list[dict]],
    common_grouped: dict[str, list[dict]],
    category_order: list[str],
    main_pdf: Path,
    reference_pdf: Path,
) -> None:
    h9._write_manifest(
        output_path=output_path,
        mode=mode,
        patient=patient,
        state=state,
        root=root,
        source_files=source_files,
        research_layer_opened=research_layer_opened,
        main_grouped=main_grouped,
        common_grouped=common_grouped,
        category_order=category_order,
        main_pdf=main_pdf,
        reference_pdf=reference_pdf,
    )
    manifest = json.loads(output_path.read_text(encoding="utf-8"))
    manifest["publication_contract"] = PUBLICATION_CONTRACT
    manifest["report_version"] = version
    manifest["version_suffix"] = _version_suffix(version)
    manifest["versioning_policy"] = {
        "rule": VERSIONING_RULE,
        "overwrite_existing_outputs": False,
        "partial_package_consumes_version": True,
    }
    output_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def self_test() -> int:
    h9.self_test()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        version, main_pdf, reference_pdf, manifest = allocate_report_package(root, "REPORT")
        assert version == 1
        assert main_pdf.name == "REPORT.pdf"
        assert reference_pdf.name == "REPORT_COMMON_MAJOR_REFERENCE.pdf"
        assert manifest.name == "REPORT.manifest.json"

        # Any existing member consumes the whole version, including a partial failed run.
        main_pdf.write_bytes(b"locked-or-existing")
        version, main_pdf, reference_pdf, manifest = allocate_report_package(root, "REPORT")
        assert version == 2
        assert main_pdf.name == "REPORT_v2.pdf"
        assert reference_pdf.name == "REPORT_COMMON_MAJOR_REFERENCE_v2.pdf"
        assert manifest.name == "REPORT_v2.manifest.json"

        manifest.write_text("{}", encoding="utf-8")
        version, main_pdf, reference_pdf, manifest = allocate_report_package(root, "REPORT")
        assert version == 3
        assert main_pdf.name == "REPORT_v3.pdf"
        assert reference_pdf.name == "REPORT_COMMON_MAJOR_REFERENCE_v3.pdf"
        assert manifest.name == "REPORT_v3.manifest.json"
    print("SELF-TEST PASS | non-overwriting report package version allocation")
    return 0


def main() -> int:
    args = h9.parse_args()
    if args.mode == "self-test":
        return self_test()
    if not args.root or not args.patient_root:
        raise RuntimeError("--root and --patient-root are required")

    root = args.root.resolve()
    patient_root = args.patient_root.resolve()
    patient = patient_root.name
    knowledge = patient_root / "Dante_results" / "knowledge_enrichment"
    state_path = knowledge / "Dante_knowledge_enrichment.state.json"
    open_path = knowledge / "Dante_knowledge_findings_open.tsv"
    research_path = knowledge / "Dante_knowledge_findings_research.tsv"
    open_usage_path = knowledge / "Dante_knowledge_usage_manifest.json"
    research_usage_path = knowledge / "Dante_knowledge_usage_research_manifest.json"

    for required in (state_path, open_path, open_usage_path):
        if not required.is_file():
            raise FileNotFoundError(required)
    state = h9.base.read_json(state_path)
    if state.get("status") != "complete" or state.get("schema") != h9.STAGE24_SCHEMA or state.get("rules_schema") != h9.RULES_SCHEMA:
        raise RuntimeError("Stage 24 is not READY under rules-v8")

    out_root = patient_root / "Dante_results" / "reports"
    out_root.mkdir(parents=True, exist_ok=True)
    stage_rows = h9.base.stage_presence(root, patient_root)
    open_rows = h9.base.read_tsv(open_path)
    open_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in open_rows:
        open_by_event[h9.base.first(row, "event_key") or h9.base.first(row, "finding_id")].append(row)

    if args.mode == "commercial":
        usage = h9.base.read_json(open_usage_path)
        h9.base.validate_usage("commercial", usage)
        all_grouped = h9.base.group_reason_events(open_rows)
        main_grouped, common_grouped = h9.split_common_major_events(all_grouped)
        version, main_pdf, reference_pdf, manifest_path = allocate_report_package(
            out_root, "Dante_genome_report_COMMERCIAL"
        )
        h9.render_report_pdf(
            main_pdf,
            title="Dante WGS - итоговый геномный отчет",
            patient=patient, mode="commercial", state=state, usage=usage, stage_rows=stage_rows,
            grouped=main_grouped, category_meta=h9.OPEN_CATEGORY_META, category_order=h9.OPEN_CATEGORY_ORDER,
            excluded_grouped=common_grouped, reference_appendix_name=reference_pdf.name,
        )
        h9.render_report_pdf(
            reference_pdf,
            title="Dante WGS - справочное приложение: common major alleles",
            patient=patient, mode="commercial", state=state, usage=usage, stage_rows=stage_rows,
            grouped=common_grouped, category_meta=h9.OPEN_CATEGORY_META, category_order=h9.OPEN_CATEGORY_ORDER,
            excluded_grouped=main_grouped, reference_appendix_name=main_pdf.name, is_reference_appendix=True,
        )
        source_files = {
            "state": {"name": state_path.name, "sha256": h9.sha256(state_path)},
            "open_findings": {"name": open_path.name, "sha256": h9.sha256(open_path)},
            "usage": {"name": open_usage_path.name, "sha256": h9.sha256(open_usage_path)},
        }
        _write_versioned_manifest(
            output_path=manifest_path, version=version, mode="COMMERCIAL", patient=patient, state=state, root=root,
            source_files=source_files, research_layer_opened=False, main_grouped=main_grouped,
            common_grouped=common_grouped, category_order=h9.OPEN_CATEGORY_ORDER,
            main_pdf=main_pdf, reference_pdf=reference_pdf,
        )
        print(f"Created COMMERCIAL report package v{version}: {main_pdf.name}")
        print(f"Created common-major reference appendix: {reference_pdf.name}")
        print(f"Created manifest: {manifest_path.name}")
        print(
            f"Publication filter {h9.COMMON_MAJOR_RULE}: main={h9._count_events(main_grouped)} events/{h9._count_rows(main_grouped)} rows; "
            f"reference={h9._count_events(common_grouped)} events/{h9._count_rows(common_grouped)} rows. "
            "Existing outputs were not overwritten; Stage 24 unchanged; OpenAI/API not used."
        )
        return 0

    if not args.accept_research_restrictions:
        if not sys.stdin.isatty():
            raise RuntimeError("Research mode requires --accept-research-restrictions")
        print("RESEARCH ONLY: non-commercial; not for diagnosis, medical decisions or healthcare-service use.")
        if input("Type RESEARCH to continue: ").strip() != "RESEARCH":
            raise RuntimeError("Research report cancelled")
    if not research_path.is_file() or not research_usage_path.is_file():
        raise FileNotFoundError("Research Stage-24 outputs are missing")

    usage = h9.base.read_json(research_usage_path)
    h9.base.validate_usage("research", usage)
    overlays = h9.base.read_tsv(research_path)
    all_grouped = h9.base.group_reason_events(overlays)
    main_grouped, common_grouped = h9.split_common_major_events(all_grouped, open_by_event)
    version, main_pdf, reference_pdf, manifest_path = allocate_report_package(
        out_root, "Dante_genome_report_RESEARCH_ONLY"
    )
    h9.render_report_pdf(
        main_pdf,
        title="Dante WGS - исследовательский геномный отчет",
        patient=patient, mode="research", state=state, usage=usage, stage_rows=stage_rows,
        grouped=main_grouped, category_meta=h9.RESEARCH_CATEGORY_META, category_order=h9.RESEARCH_CATEGORY_ORDER,
        excluded_grouped=common_grouped, reference_appendix_name=reference_pdf.name, open_by_event=open_by_event,
    )
    h9.render_report_pdf(
        reference_pdf,
        title="Dante WGS - RESEARCH ONLY: common major reference appendix",
        patient=patient, mode="research", state=state, usage=usage, stage_rows=stage_rows,
        grouped=common_grouped, category_meta=h9.RESEARCH_CATEGORY_META, category_order=h9.RESEARCH_CATEGORY_ORDER,
        excluded_grouped=main_grouped, reference_appendix_name=main_pdf.name,
        is_reference_appendix=True, open_by_event=open_by_event,
    )
    source_files = {
        "state": {"name": state_path.name, "sha256": h9.sha256(state_path)},
        "open_findings": {"name": open_path.name, "sha256": h9.sha256(open_path)},
        "research_overlays": {"name": research_path.name, "sha256": h9.sha256(research_path)},
        "usage": {"name": research_usage_path.name, "sha256": h9.sha256(research_usage_path)},
    }
    _write_versioned_manifest(
        output_path=manifest_path, version=version, mode="RESEARCH_ONLY", patient=patient, state=state, root=root,
        source_files=source_files, research_layer_opened=True, main_grouped=main_grouped,
        common_grouped=common_grouped, category_order=h9.RESEARCH_CATEGORY_ORDER,
        main_pdf=main_pdf, reference_pdf=reference_pdf,
    )
    print(f"Created RESEARCH ONLY report package v{version}: {main_pdf.name}")
    print(f"Created research common-major reference appendix: {reference_pdf.name}")
    print(f"Created manifest: {manifest_path.name}")
    print(
        f"Publication filter {h9.COMMON_MAJOR_RULE}: main={h9._count_events(main_grouped)} events/{h9._count_rows(main_grouped)} rows; "
        f"reference={h9._count_events(common_grouped)} events/{h9._count_rows(common_grouped)} rows. "
        "Existing outputs were not overwritten; Stage 24 unchanged; OpenAI/API not used."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
