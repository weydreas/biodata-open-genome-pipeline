#!/usr/bin/env python3
"""Adapters for GRCh37 23andMe-style genotype tables.

The input contract is intentionally weaker than VCF/WGS: rows contain an ID,
coordinate and genotype but no explicit REF/ALT.  Sequence-level annotation is
therefore limited to A/C/G/T calls that can be resolved against the installed
GRCh37 FASTA.  Symbolic I/D calls are never converted into invented sequence
alleles.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

GRCH37_CONTIG_LENGTHS = {
    "1": 249250621,
    "2": 243199373,
    "3": 198022430,
    "4": 191154276,
    "5": 180915260,
    "6": 171115067,
    "7": 159138663,
    "8": 146364022,
    "9": 141213431,
    "10": 135534747,
    "11": 135006516,
    "12": 133851895,
    "13": 115169878,
    "14": 107349540,
    "15": 102531392,
    "16": 90354753,
    "17": 81195210,
    "18": 78077248,
    "19": 59128983,
    "20": 63025520,
    "21": 48129895,
    "22": 51304566,
    "X": 155270560,
    "Y": 59373566,
    "MT": 16569,
}
RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)


@dataclass(frozen=True)
class TableRow:
    marker: str
    chrom: str
    pos: int
    genotype: str
    line_no: int


@dataclass(frozen=True)
class FaiEntry:
    length: int
    offset: int
    line_bases: int
    line_width: int


def normalize_chrom(value: str) -> str:
    text = value.strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.upper() == "M":
        return "MT"
    if text.upper() in {"X", "Y", "MT"}:
        return text.upper()
    return text


def parse_rows(path: Path) -> tuple[list[str], list[TableRow], Counter[str]]:
    comments: list[str] = []
    rows: list[TableRow] = []
    counts: Counter[str] = Counter()
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if line.startswith("#"):
                comments.append(line)
                continue
            if not line.strip():
                continue
            counts["input_rows"] += 1
            fields = line.split("\t")
            if len(fields) != 4:
                counts["malformed_rows"] += 1
                continue
            marker, chrom_raw, pos_raw, genotype_raw = fields
            chrom = normalize_chrom(chrom_raw)
            try:
                pos = int(pos_raw)
            except ValueError:
                counts["invalid_position_rows"] += 1
                continue
            if chrom not in GRCH37_CONTIG_LENGTHS or pos <= 0 or pos > GRCH37_CONTIG_LENGTHS[chrom]:
                counts["invalid_coordinate_rows"] += 1
                continue
            genotype = genotype_raw.strip().upper()
            rows.append(TableRow(marker.strip() or ".", chrom, pos, genotype, line_no))
            counts["valid_coordinate_rows"] += 1
            if genotype == "--":
                counts["no_call_rows"] += 1
            elif genotype and set(genotype) <= set("ACGT") and len(genotype) in {1, 2}:
                counts["sequence_call_rows"] += 1
                counts["haploid_sequence_rows" if len(genotype) == 1 else "diploid_sequence_rows"] += 1
            elif genotype and set(genotype) <= set("ID") and len(genotype) in {1, 2}:
                counts["symbolic_indel_rows"] += 1
            else:
                counts["unsupported_genotype_rows"] += 1
    return comments, rows, counts


def atomic_gzip_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    with part.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="\n") as handle:
                handle.write(text)
    os.replace(part, path)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(text, encoding="utf-8", newline="\n")
    os.replace(part, path)


def export_23andme(source: Path, output: Path, report: Path) -> dict[str, int]:
    comments, rows, counts = parse_rows(source)
    out_lines = [
        "# Dante measured-loci genotype export from a 23andMe-style GRCh37 table",
        "# This is NOT whole-genome sequencing. Absence outside rows present in the source is UNKNOWN.",
        "# -- is preserved as a no-call. Symbolic I/D genotypes are preserved here but are not sequence-resolved for VEP.",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    for row in rows:
        out_lines.append(f"{row.marker}\t{row.chrom}\t{row.pos}\t{row.genotype}")
    atomic_gzip_text(output, "\n".join(out_lines) + "\n")
    result = dict(counts)
    result["output_rows"] = len(rows)
    result["source_comment_lines"] = len(comments)
    report_lines = [
        "DANTE GENOTYPE-TABLE DIRECT EXPORT",
        "=" * 78,
        f"Source: {source}",
        "Scope: measured rows only; no whole-genome negative evidence.",
        "No-calls and symbolic I/D genotypes are preserved exactly in this direct export.",
        "",
    ] + [f"{key}: {result[key]}" for key in sorted(result)]
    atomic_write_text(report, "\n".join(report_lines) + "\n")
    return result


def read_fai(path: Path) -> dict[str, FaiEntry]:
    entries: dict[str, FaiEntry] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        fields = raw.split("\t")
        if len(fields) < 5:
            raise ValueError(f"invalid FASTA index line: {raw}")
        chrom = normalize_chrom(fields[0])
        entries[chrom] = FaiEntry(
            length=int(fields[1]),
            offset=int(fields[2]),
            line_bases=int(fields[3]),
            line_width=int(fields[4]),
        )
    return entries


def fasta_byte_offset(entry: FaiEntry, pos: int) -> int:
    zero = pos - 1
    return entry.offset + (zero // entry.line_bases) * entry.line_width + (zero % entry.line_bases)


def fetch_reference_bases(
    fasta: Path,
    fai: Path,
    rows: Iterable[TableRow],
    block_size: int = 4 * 1024 * 1024,
) -> dict[int, str]:
    entries = read_fai(fai)
    requests: list[tuple[int, int]] = []
    for index, row in enumerate(rows):
        if not (row.genotype and set(row.genotype) <= set("ACGT") and len(row.genotype) in {1, 2}):
            continue
        entry = entries.get(row.chrom)
        if entry is None or row.pos > entry.length:
            continue
        requests.append((fasta_byte_offset(entry, row.pos), index))
    requests.sort()

    bases: dict[int, str] = {}
    with fasta.open("rb") as handle:
        cached_start = -1
        cached = b""
        for offset, index in requests:
            start = (offset // block_size) * block_size
            if not (cached_start <= offset < cached_start + len(cached)):
                handle.seek(start)
                cached = handle.read(block_size)
                cached_start = start
            relative = offset - cached_start
            if relative < 0 or relative >= len(cached):
                continue
            base = chr(cached[relative]).upper()
            bases[index] = base
    return bases


def choose_marker(ids: Iterable[str]) -> str:
    unique = sorted({x for x in ids if x and x != "."})
    rsids = [x for x in unique if RSID_RE.fullmatch(x)]
    if rsids:
        return rsids[0]
    return unique[0] if unique else "."


def vcf_header() -> list[str]:
    lines = [
        "##fileformat=VCFv4.2",
        "##source=DanteGenotypeTableAdapter",
        "##reference=GRCh37",
        '##FILTER=<ID=PASS,Description="Reference-resolved genotype-table call accepted for deterministic annotation; not read-level sequencing QC">',
        '##INFO=<ID=DANTE_GENOTYPE_TXT_SOURCE,Number=0,Type=Flag,Description="Derived from a measured 23andMe-style genotype-table row">',
        '##INFO=<ID=DANTE_SOURCE_ALIAS_COUNT,Number=1,Type=Integer,Description="Number of source marker aliases collapsed at this coordinate/genotype">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Reference-resolved genotype">',
    ]
    for chrom, length in GRCH37_CONTIG_LENGTHS.items():
        lines.append(f"##contig=<ID={chrom},length={length},assembly=GRCh37>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGENOTYPE_TABLE")
    return lines


def write_vcf_gz(path: Path, records: Iterable[str]) -> None:
    atomic_gzip_text(path, "\n".join(vcf_header() + list(records)) + "\n")


def prepare_layers(
    source: Path,
    fasta: Path,
    fai: Path,
    snp_out: Path,
    indel_out: Path,
    residual_out: Path,
    report: Path,
) -> dict[str, int]:
    if not fasta.is_file():
        raise FileNotFoundError(f"GRCh37 FASTA is missing: {fasta}")
    if not fai.is_file():
        raise FileNotFoundError(f"GRCh37 FASTA index is missing: {fai}")

    _, rows, counts = parse_rows(source)
    bases = fetch_reference_bases(fasta, fai, rows)

    # Coordinate -> event key -> source IDs. Multiple aliases with the same
    # genotype collapse; conflicting genotype assertions at one coordinate are
    # excluded rather than guessed.
    events_by_coord: dict[tuple[str, int], dict[tuple[str, tuple[str, ...], str], list[str]]] = defaultdict(lambda: defaultdict(list))
    for index, row in enumerate(rows):
        genotype = row.genotype
        if genotype == "--":
            continue
        if genotype and set(genotype) <= set("ID") and len(genotype) in {1, 2}:
            continue
        if not (genotype and set(genotype) <= set("ACGT") and len(genotype) in {1, 2}):
            continue
        ref = bases.get(index, "").upper()
        if ref not in {"A", "C", "G", "T"}:
            counts["reference_base_unavailable_rows"] += 1
            continue
        alleles = list(genotype)
        if all(allele == ref for allele in alleles):
            counts["measured_homozygous_reference_rows"] += 1
            continue
        alts = tuple(sorted({allele for allele in alleles if allele != ref}))
        if not alts:
            counts["reference_resolution_error_rows"] += 1
            continue
        indexes = {allele: i + 1 for i, allele in enumerate(alts)}
        gt_tokens = ["0" if allele == ref else str(indexes[allele]) for allele in alleles]
        gt = gt_tokens[0] if len(gt_tokens) == 1 else "/".join(gt_tokens)
        key = (ref, alts, gt)
        events_by_coord[(row.chrom, row.pos)][key].append(row.marker)
        counts["sequence_variant_source_rows"] += 1

    records: list[tuple[str, int, str]] = []
    for (chrom, pos), variants in events_by_coord.items():
        if len(variants) != 1:
            counts["conflicting_coordinates_excluded"] += 1
            counts["conflicting_source_rows_excluded"] += sum(len(ids) for ids in variants.values())
            continue
        (ref, alts, gt), ids = next(iter(variants.items()))
        counts["duplicate_alias_rows_collapsed"] += max(0, len(ids) - 1)
        if len(alts) > 1:
            counts["multiallelic_carried_sites"] += 1
        marker = choose_marker(ids)
        info = f"DANTE_GENOTYPE_TXT_SOURCE;DANTE_SOURCE_ALIAS_COUNT={len(ids)}"
        line = f"{chrom}\t{pos}\t{marker}\t{ref}\t{','.join(alts)}\t.\tPASS\t{info}\tGT\t{gt}"
        records.append((chrom, pos, line))

    chrom_order = {chrom: i for i, chrom in enumerate(GRCH37_CONTIG_LENGTHS)}
    records.sort(key=lambda item: (chrom_order[item[0]], item[1], item[2]))
    write_vcf_gz(snp_out, (line for _, _, line in records))
    write_vcf_gz(indel_out, [])
    write_vcf_gz(residual_out, [])

    counts["carried_sequence_snp_sites"] = len(records)
    counts["sequence_resolved_indel_sites"] = 0
    counts["residual_sites"] = 0
    result = dict(counts)
    report_lines = [
        "DANTE GENOTYPE-TABLE REFERENCE-RESOLVED PREPARATION",
        "=" * 78,
        f"Source: {source}",
        f"Reference: {fasta}",
        "Only measured A/C/G/T calls carrying a non-reference allele enter sequence-level annotation.",
        "Haploid A/C/G/T calls on X/Y/MT are supported as haploid VCF GT values.",
        "-- no-calls are UNKNOWN and excluded from annotation.",
        "I/D symbolic indels are preserved in direct export/QC but excluded because no inserted/deleted sequence is encoded.",
        "Duplicate aliases with an identical coordinate/genotype collapse; conflicting assertions at one coordinate are excluded.",
        "Derived FILTER=PASS means accepted source genotype after exact GRCh37 REF resolution, not read-level sequencing QC.",
        "Downstream bcftools uses --check-ref e: any adapter/reference disagreement is fatal, never strand-flipped or silently removed.",
        "",
    ] + [f"{key}: {result[key]}" for key in sorted(result)]
    atomic_write_text(report, "\n".join(report_lines) + "\n")
    return result


def _write_test_fasta(fasta: Path, fai: Path) -> None:
    sequences = {"1": "ACGTACGTACGTACGTACGTACGTACGTACGT", "X": "TGCA", "Y": "GACT", "MT": "ACGT"}
    offset = 0
    fai_lines: list[str] = []
    with fasta.open("wb") as handle:
        for name, seq in sequences.items():
            header = f">{name}\n".encode()
            handle.write(header)
            offset += len(header)
            seq_bytes = (seq + "\n").encode()
            handle.write(seq_bytes)
            fai_lines.append(f"{name}\t{len(seq)}\t{offset}\t{len(seq)}\t{len(seq)+1}")
            offset += len(seq_bytes)
    fai.write_text("\n".join(fai_lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)

    export = sub.add_parser("export-23andme")
    export.add_argument("source", type=Path)
    export.add_argument("--output", required=True, type=Path)
    export.add_argument("--report", required=True, type=Path)

    prep = sub.add_parser("prepare-layers")
    prep.add_argument("source", type=Path)
    prep.add_argument("--fasta", required=True, type=Path)
    prep.add_argument("--fai", required=True, type=Path)
    prep.add_argument("--snp-out", required=True, type=Path)
    prep.add_argument("--indel-out", required=True, type=Path)
    prep.add_argument("--residual-out", required=True, type=Path)
    prep.add_argument("--report", required=True, type=Path)

    sub.add_parser("self-test")
    args = ap.parse_args()

    if args.mode == "export-23andme":
        export_23andme(args.source, args.output, args.report)
        return 0
    if args.mode == "prepare-layers":
        prepare_layers(
            args.source,
            args.fasta,
            args.fai,
            args.snp_out,
            args.indel_out,
            args.residual_out,
            args.report,
        )
        return 0

    import tempfile

    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        source = root / "sample_normal.txt"
        source.write_text(
            "# This data file generated by 23andMe at: synthetic\n"
            "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            "rs1\t1\t1\tAA\n"       # hom-ref
            "rs2\t1\t2\tCT\n"       # ref C, carried T
            "i2\t1\t2\tCT\n"        # duplicate alias, collapses
            "rs3\t1\t3\tAA\n"       # ref G, hom-alt A
            "rs4\tX\t1\tT\n"        # haploid hom-ref
            "rs5\tX\t2\tA\n"        # haploid carried vs G
            "rs6\t1\t4\tII\n"       # unresolved symbolic indel
            "rs7\t1\t5\t--\n",      # no-call
            encoding="utf-8",
        )
        fasta = root / "ref.fa"
        fai = root / "ref.fa.fai"
        _write_test_fasta(fasta, fai)
        export_path = root / "export.txt.gz"
        export_23andme(source, export_path, root / "export_report.txt")
        with gzip.open(export_path, "rt", encoding="utf-8") as handle:
            direct = handle.read()
        assert "rs6\t1\t4\tII" in direct
        assert "rs7\t1\t5\t--" in direct

        result = prepare_layers(
            source,
            fasta,
            fai,
            root / "snp.vcf.gz",
            root / "indel.vcf.gz",
            root / "residual.vcf.gz",
            root / "prep.txt",
        )
        with gzip.open(root / "snp.vcf.gz", "rt", encoding="utf-8") as handle:
            snp_text = handle.read()
        with gzip.open(root / "indel.vcf.gz", "rt", encoding="utf-8") as handle:
            indel_text = handle.read()
        assert "rs1" not in snp_text
        assert "rs2\tC\tT" in snp_text
        assert "\tGT\t0/1" in snp_text
        assert "rs3\tG\tA" in snp_text
        assert "\tGT\t1/1" in snp_text
        assert "rs5\tG\tA" in snp_text
        assert "\tGT\t1\n" in snp_text
        assert "rs6" not in snp_text
        assert result["symbolic_indel_rows"] == 1
        assert result["duplicate_alias_rows_collapsed"] == 1
        assert result["carried_sequence_snp_sites"] == 3
        assert "\n#CHROM" in indel_text and "\n1\t" not in indel_text

    print("SELF-TEST PASS | genotype TXT passthrough + GRCh37-resolved haploid/diploid SNP adapter")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
