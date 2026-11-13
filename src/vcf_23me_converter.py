#!/usr/bin/env python3
"""Loss-accounted single-sample VCF -> 23andMe-column-compatible SNP export.

Contract (v9.1.0):
- input: one-sample .vcf or .vcf.gz;
- output: rsid/marker, chromosome, position, genotype columns;
- only canonical A/C/G/T SNV loci are representable in this export;
- every valid CALLED canonical SNV is exported exactly at record level;
- missing rsIDs never cause silent loss: deterministic DANTE_* marker IDs are used;
- optional batched NCBI Variation Services lookup enriches missing/current rsIDs;
- source assembly coordinates are preserved. GRCh37 matches the coordinate convention
  used by 23andMe v5 raw data; GRCh38 is explicitly labelled and is never silently
  lifted to GRCh37.

The unresolved/audit TSV records every source record that cannot be represented and
every called SNP that had to use a DANTE marker. This makes loss accounting observable.
"""
from __future__ import annotations

import argparse
import gzip
import re
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

PIPELINE_RELEASE = "9.1.0"
RSID_RE = re.compile(r"(?<![A-Za-z0-9])rs\d+", re.IGNORECASE)
CANONICAL = {*(str(i) for i in range(1, 23)), "X", "Y", "M", "MT"}
BASES = {"A", "C", "G", "T"}
ASSEMBLY_ACCESSION = {
    "GRCh37": "GCF_000001405.25",
    "GRCh38": "GCF_000001405.40",
}
NCBI_SET_RSIDS = "https://api.ncbi.nlm.nih.gov/variation/v0/vcf/file/set_rsids"


@dataclass
class Record:
    chrom: str
    pos: str
    vid: str
    ref: str
    alt_field: str
    filt: str
    fmt: str
    sample: str
    raw_fields: list[str]

    @property
    def alts(self) -> list[str]:
        return self.alt_field.split(",")


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def normalize_chrom(value: str) -> str | None:
    chrom = value[3:] if value.lower().startswith("chr") else value
    if chrom not in CANONICAL:
        return None
    return "MT" if chrom == "M" else chrom


def normalize_assembly(value: str) -> str:
    text = (value or "").strip().lower()
    if any(token in text for token in ("grch37", "hg19", "b37", "gcf_000001405.25")):
        return "GRCh37"
    if any(token in text for token in ("grch38", "hg38", "gcf_000001405.40")):
        return "GRCh38"
    return "UNKNOWN"


def detect_assembly(meta_reference: str, explicit: str) -> str:
    if explicit != "auto":
        return explicit
    detected = normalize_assembly(meta_reference)
    if detected == "UNKNOWN":
        raise RuntimeError(
            "VCF assembly is unknown. Pass --assembly GRCh37 or --assembly GRCh38; "
            "the converter will not guess genomic coordinates."
        )
    return detected


def parse_gt(fmt: str, sample: str) -> str:
    keys = fmt.split(":") if fmt else []
    values = sample.split(":") if sample else []
    try:
        index = keys.index("GT")
    except ValueError:
        return "."
    return values[index] if index < len(values) else "."


def gt_indices(gt: str) -> list[int] | None:
    if gt in {"", ".", "./.", ".|."}:
        return None
    tokens = re.split(r"[/|]", gt)
    if not tokens or any(not token.isdigit() for token in tokens):
        return None
    return [int(token) for token in tokens]


def genotype_bases(ref: str, alts: list[str], indices: list[int]) -> str | None:
    # 23andMe-style raw genotype columns can represent haploid or diploid calls.
    if len(indices) not in (1, 2):
        return None
    alleles = [ref, *alts]
    if any(index < 0 or index >= len(alleles) for index in indices):
        return None
    called = [alleles[index].upper() for index in indices]
    if any(len(allele) != 1 or allele not in BASES for allele in called):
        return None
    return "".join(called)


def rsids(value: str) -> list[str]:
    found = {match.lower() for match in RSID_RE.findall(value or "")}
    return sorted(found, key=lambda item: int(item[2:]))


def fallback_marker(assembly: str, chrom: str, pos: str, ref: str, alts: list[str]) -> str:
    alleles = "_".join([ref.upper(), *[alt.upper() for alt in alts]])
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", alleles).strip("_")
    return f"DANTE_{assembly}_{chrom}_{pos}_{safe}"


def conflict_marker(base: str, ordinal: int) -> str:
    return f"{base}_CONFLICT_{ordinal}"


def ncbi_set_rsids(records: list[Record], assembly: str, timeout: int = 60) -> list[str] | None:
    """Batch-annotate VCF rows with RefSNP IDs using the public NCBI endpoint.

    NCBI's own dbSNP tutorial uses /vcf/file/set_rsids in batches. Network failure is
    non-fatal here because the lossless fallback marker contract remains available.
    """
    if not records:
        return []
    lines: list[str] = []
    for record in records:
        fields = list(record.raw_fields[:8])
        while len(fields) < 8:
            fields.append(".")
        fields[0] = record.chrom
        fields[1] = record.pos
        lines.append("\t".join(fields))
    payload = ("\n".join(lines) + "\n").encode("utf-8")
    url = NCBI_SET_RSIDS + "?assembly=" + urllib.parse.quote(ASSEMBLY_ACCESSION[assembly])
    request = urllib.request.Request(
        url,
        data=payload,
        method="POST",
        headers={
            "Content-Type": "text/plain; charset=utf-8",
            "User-Agent": f"BiodataGenomePipeline/{PIPELINE_RELEASE}",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            text = response.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError, OSError):
        return None
    rows = [line for line in text.splitlines() if line and not line.startswith("#")]
    if len(rows) != len(records):
        return None
    values: list[str] = []
    for row in rows:
        fields = row.split("\t")
        values.append(fields[2] if len(fields) >= 3 else ".")
    return values


def write_audit(
    handle: TextIO,
    record: Record,
    gt: str,
    reason: str,
    marker: str = ".",
) -> None:
    handle.write(
        "\t".join(
            [
                record.chrom,
                record.pos,
                record.ref,
                record.alt_field,
                gt or ".",
                record.vid or ".",
                marker or ".",
                reason,
            ]
        )
        + "\n"
    )


def convert(
    source: Path,
    output: Path,
    audit_path: Path,
    report: Path,
    *,
    assembly: str = "auto",
    resolve_rsids: bool = True,
    batch_size: int = 500,
    timeout: int = 60,
) -> int:
    source = source.resolve()
    output = output.resolve()
    audit_path = audit_path.resolve()
    report = report.resolve()
    for path in (output, audit_path, report):
        path.parent.mkdir(parents=True, exist_ok=True)

    output_part = output.with_name(output.name + ".part")
    audit_part = audit_path.with_name(audit_path.name + ".part")
    counts: Counter[str] = Counter()
    meta_reference = ""
    assembly_final = "UNKNOWN"
    buffered: list[tuple[Record, list[int], str, str]] = []
    seen_markers: dict[str, tuple[str, str, str]] = {}

    def emit_marker(out: TextIO, audit: TextIO, record: Record, marker: str, genotype: str) -> None:
        identity = (record.chrom, record.pos, genotype)
        existing = seen_markers.get(marker)
        if existing is None:
            seen_markers[marker] = identity
            out.write(f"{marker}\t{record.chrom}\t{record.pos}\t{genotype}\n")
            counts["output_rows"] += 1
            return
        if existing == identity:
            counts["deduplicated_output_rows"] += 1
            return
        ordinal = counts["marker_conflict_rows"] + 1
        replacement = conflict_marker(
            fallback_marker(assembly_final, record.chrom, record.pos, record.ref, record.alts),
            ordinal,
        )
        counts["marker_conflict_rows"] += 1
        seen_markers[replacement] = identity
        out.write(f"{replacement}\t{record.chrom}\t{record.pos}\t{genotype}\n")
        counts["output_rows"] += 1
        write_audit(audit, record, parse_gt(record.fmt, record.sample), "RSID_OR_MARKER_CONFLICT_EXPORTED_WITH_DANTE_MARKER", replacement)

    def flush(out: TextIO, audit: TextIO) -> None:
        nonlocal buffered
        if not buffered:
            return
        records = [item[0] for item in buffered]
        remote_ids: list[str] | None = None
        if resolve_rsids and any(not rsids(record.vid) for record in records):
            remote_ids = ncbi_set_rsids(records, assembly_final, timeout=timeout)
            if remote_ids is None:
                counts["remote_resolution_failed_batches"] += 1
            else:
                counts["remote_resolution_batches"] += 1

        for index, (record, indices, genotype, raw_gt) in enumerate(buffered):
            source_ids = rsids(record.vid)
            remote = rsids(remote_ids[index]) if remote_ids is not None else []
            markers = remote or source_ids
            if remote and not source_ids:
                counts["rsid_resolved_remote_records"] += 1
            elif source_ids:
                counts["rsid_from_source_records"] += 1

            if not markers:
                marker = fallback_marker(assembly_final, record.chrom, record.pos, record.ref, record.alts)
                markers = [marker]
                counts["fallback_marker_records"] += 1
                write_audit(audit, record, raw_gt, "NO_RSID_EXPORTED_WITH_DANTE_MARKER", marker)

            for marker in markers:
                emit_marker(out, audit, record, marker, genotype)

            counts["exported_called_snv_records"] += 1
            if "|" in raw_gt:
                counts["phased_gt_collapsed_to_unphased_columns"] += 1
            if len(indices) == 1:
                counts["haploid_called_snv_records"] += 1
            elif all(value == 0 for value in indices):
                counts["homozygous_reference_called_snv_records"] += 1
            elif len(set(indices)) == 1:
                counts["homozygous_alt_called_snv_records"] += 1
            else:
                counts["heterozygous_called_snv_records"] += 1
            if len(record.alts) > 1:
                counts["multiallelic_called_snv_records"] += 1
            if record.filt not in {"", ".", "PASS"}:
                counts["source_filter_nonpass_exported_records"] += 1
        buffered = []

    try:
        with open_text(source) as handle:
            sample_count = 0
            for raw in handle:
                if raw.startswith("##reference="):
                    meta_reference = raw.split("=", 1)[1].strip()
                    continue
                if raw.startswith("#CHROM"):
                    columns = raw.rstrip("\r\n").split("\t")
                    sample_count = max(0, len(columns) - 9)
                    if sample_count != 1:
                        raise RuntimeError(f"Expected exactly one sample column, found {sample_count}")
                    break
                if not raw.startswith("#") and raw.strip():
                    raise RuntimeError("VCF #CHROM header was not found before records")
            else:
                raise RuntimeError("VCF #CHROM header was not found")

            assembly_final = detect_assembly(meta_reference, assembly)
            with output_part.open("wt", encoding="utf-8", newline="") as out, audit_part.open(
                "wt", encoding="utf-8", newline=""
            ) as audit:
                out.write(f"# Biodata Genome Pipeline {PIPELINE_RELEASE}; VCF -> 23andMe v5-compatible columns\n")
                out.write(f"# reference_build={assembly_final}\n")
                if assembly_final == "GRCh37":
                    out.write("# coordinates=GRCh37; coordinate-compatible with 23andMe v5 raw data\n")
                else:
                    out.write("# coordinates=GRCh38; column-compatible only; no silent GRCh37 liftover was performed\n")
                out.write("# rsid\tchromosome\tposition\tgenotype\n")
                audit.write("chromosome\tposition\tref\talt\tgt\tvcf_id\texported_marker\treason\n")

                for raw in handle:
                    if raw.startswith("#") or not raw.strip():
                        continue
                    fields = raw.rstrip("\r\n").split("\t")
                    if len(fields) < 10:
                        counts["malformed_records"] += 1
                        continue
                    counts["input_records"] += 1
                    chrom_raw, pos, vid, ref, alt_field, _, filt, _, fmt, sample = fields[:10]
                    chrom = normalize_chrom(chrom_raw)
                    record = Record(
                        chrom=chrom or chrom_raw,
                        pos=pos,
                        vid=vid,
                        ref=ref.upper(),
                        alt_field=alt_field.upper(),
                        filt=filt,
                        fmt=fmt,
                        sample=sample,
                        raw_fields=fields,
                    )
                    if chrom is None:
                        counts["excluded_noncanonical_chromosome"] += 1
                        write_audit(audit, record, parse_gt(fmt, sample), "NONCANONICAL_CHROMOSOME")
                        continue
                    record.chrom = chrom
                    try:
                        if int(pos) <= 0:
                            raise ValueError
                    except ValueError:
                        counts["excluded_invalid_position"] += 1
                        write_audit(audit, record, parse_gt(fmt, sample), "INVALID_POSITION")
                        continue

                    alts = record.alts
                    if (
                        record.ref not in BASES
                        or len(record.ref) != 1
                        or not alts
                        or any(len(alt) != 1 or alt not in BASES for alt in alts)
                    ):
                        counts["excluded_non_snv_records"] += 1
                        write_audit(audit, record, parse_gt(fmt, sample), "NON_SNV")
                        continue

                    raw_gt = parse_gt(fmt, sample)
                    indices = gt_indices(raw_gt)
                    if indices is None:
                        counts["uncalled_snv_records"] += 1
                        write_audit(audit, record, raw_gt, "MISSING_OR_INVALID_GT")
                        continue
                    genotype = genotype_bases(record.ref, alts, indices)
                    if genotype is None:
                        counts["unrepresentable_called_snv_records"] += 1
                        write_audit(audit, record, raw_gt, "GT_NOT_REPRESENTABLE_IN_23ANDME_COLUMNS")
                        continue

                    counts["called_canonical_snv_records"] += 1
                    buffered.append((record, indices, genotype, raw_gt))
                    if len(buffered) >= max(1, batch_size):
                        flush(out, audit)
                flush(out, audit)

        output_part.replace(output)
        audit_part.replace(audit_path)
    except Exception:
        output_part.unlink(missing_ok=True)
        audit_part.unlink(missing_ok=True)
        raise

    called = counts["called_canonical_snv_records"]
    exported = counts["exported_called_snv_records"]
    if called != exported:
        raise RuntimeError(
            f"Loss-accounting invariant failed: called canonical SNV records={called}, exported={exported}"
        )
    if counts["output_rows"] == 0:
        output.unlink(missing_ok=True)
        raise RuntimeError("No called canonical SNP genotypes were exported")

    lines = [
        "DANTE VCF -> 23ANDME V5-COMPATIBLE SNP EXPORT",
        "=" * 78,
        f"Pipeline release: {PIPELINE_RELEASE}",
        f"Source: {source}",
        f"Assembly: {assembly_final}",
        f"Remote dbSNP rsID resolution: {'enabled' if resolve_rsids else 'disabled'}",
        "Coordinate policy: source assembly is preserved; GRCh38 is never silently lifted to GRCh37.",
        "Loss policy: every valid called canonical SNV is exported; missing rsID uses a deterministic DANTE_* marker.",
        "FILTER policy: source non-PASS calls are retained and counted; this converter does not invent a sequencing QC threshold.",
        "Phase policy: 23andMe columns do not encode phase; phased GT allele values are preserved but the separator is not.",
        "",
    ]
    for key in sorted(counts):
        lines.append(f"{key}: {counts[key]}")
    lines.extend(
        [
            f"meaningful_called_snv_loss_records: {called - exported}",
            "loss_accounting_invariant: PASS",
        ]
    )
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    convert_parser = sub.add_parser("convert")
    convert_parser.add_argument("input", type=Path)
    convert_parser.add_argument("--output", type=Path, required=True)
    convert_parser.add_argument("--audit", type=Path, required=True)
    convert_parser.add_argument("--report", type=Path, required=True)
    convert_parser.add_argument("--assembly", choices=("auto", "GRCh37", "GRCh38"), default="auto")
    convert_parser.add_argument("--no-rsid-resolution", action="store_true")
    convert_parser.add_argument("--batch-size", type=int, default=500)
    convert_parser.add_argument("--timeout", type=int, default=60)
    sub.add_parser("self-test")
    args = parser.parse_args()

    if args.mode == "convert":
        return convert(
            args.input,
            args.output,
            args.audit,
            args.report,
            assembly=args.assembly,
            resolve_rsids=not args.no_rsid_resolution,
            batch_size=max(1, args.batch_size),
            timeout=max(1, args.timeout),
        )

    import tempfile

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "sample.vcf"
        source.write_text(
            "##fileformat=VCFv4.2\n"
            "##reference=GRCh37\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"
            "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
            "1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\n"
            "X\t300\trs3\tG\tA\t.\tPASS\t.\tGT\t1\n"
            "2\t400\trs4\tA\tC,G\t.\tPASS\t.\tGT\t1/2\n",
            encoding="utf-8",
        )
        output = root / "out.txt"
        audit = root / "audit.tsv"
        report = root / "report.txt"
        convert(source, output, audit, report, resolve_rsids=False)
        text = output.read_text(encoding="utf-8")
        assert "rs1\t1\t100\tAG" in text
        assert "DANTE_GRCh37_1_200_C_T\t1\t200\tTT" in text
        assert "rs3\tX\t300\tA" in text
        assert "rs4\t2\t400\tCG" in text
        assert "meaningful_called_snv_loss_records: 0" in report.read_text(encoding="utf-8")
    print("SELF-TEST PASS | VCF -> 23andMe-compatible loss-accounted SNP export")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
