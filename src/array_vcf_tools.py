#!/usr/bin/env python3
"""Adapters for dense single-sample genotyping-array VCF inputs.

These adapters preserve an important boundary:
- 23andMe-style export may contain homozygous-reference calls because the array
  explicitly interrogated those loci.
- downstream variant annotation receives only carried alternate alleles. A 0/0
  array record is not a patient variant and must never enter VEP/Stage 24.
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import re
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

RSID_RE = re.compile(r"(?<![A-Za-z0-9])rs\d+", re.IGNORECASE)
CANONICAL = {*(str(i) for i in range(1, 23)), "X", "Y", "M", "MT"}
BASES = {"A", "C", "G", "T"}

# GRCh37 primary-assembly sequence lengths. Genotek array exports can omit
# ##contig declarations entirely. Once chr-prefixed record names are normalized
# to canonical GRCh37 names, htslib/bcftools requires those IDs to exist in the
# VCF header or it refuses to encode the output BCF stream.
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


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


@contextmanager
def anonymous_gzip_text(path: Path) -> Iterator[TextIO]:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                yield text


def normalize_chrom(value: str) -> str | None:
    chrom = value[3:] if value.lower().startswith("chr") else value
    if chrom not in CANONICAL:
        return None
    return "MT" if chrom == "M" else chrom


def parse_gt(fmt: str, sample: str) -> str:
    if not fmt or fmt == ".":
        return "."
    keys = fmt.split(":")
    vals = sample.split(":")
    try:
        idx = keys.index("GT")
    except ValueError:
        return "."
    return vals[idx] if idx < len(vals) else "."


def gt_indices(gt: str) -> list[int] | None:
    if gt in {"", ".", "./.", ".|."}:
        return None
    tokens = re.split(r"[/|]", gt)
    if not tokens or any(not x.isdigit() for x in tokens):
        return None
    return [int(x) for x in tokens]


def genotype_bases(ref: str, alts: list[str], indices: list[int]) -> str | None:
    alleles = [ref, *alts]
    if any(i < 0 or i >= len(alleles) for i in indices):
        return None
    called = [alleles[i].upper() for i in indices]
    if any(len(x) != 1 or x not in BASES for x in called):
        return None
    return "".join(called)


def rsids(value: str) -> list[str]:
    found = {x.lower() for x in RSID_RE.findall(value or "")}
    return sorted(found, key=lambda x: int(x[2:]))


def write_report(path: Path, title: str, counts: Counter[str], extra: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [title, "=" * 78]
    if extra:
        lines.extend(extra)
        lines.append("")
    for key in sorted(counts):
        lines.append(f"{key}: {counts[key]}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_23andme(source: Path, output: Path, report: Path) -> int:
    source = source.resolve()
    output = output.resolve()
    report = report.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    report.parent.mkdir(parents=True, exist_ok=True)
    part = output.with_name(output.name + ".part")
    counts: Counter[str] = Counter()
    sample_count = 0

    try:
        with open_text(source) as inp, anonymous_gzip_text(part) as out:
            out.write("# Dante genotyping-array export; GRCh37 forward-strand alleles\n")
            out.write("# Scope: only loci explicitly present in the source array VCF; this is not whole-genome coverage\n")
            out.write("# Source QUAL/FILTER/DP/GQ may be unavailable for vendor array calls\n")
            out.write("# rsid\tchromosome\tposition\tgenotype\n")
            for raw in inp:
                if raw.startswith("#CHROM"):
                    cols = raw.rstrip("\r\n").split("\t")
                    sample_count = max(0, len(cols) - 9)
                    if sample_count != 1:
                        raise RuntimeError(f"Expected exactly one sample column, found {sample_count}")
                    continue
                if raw.startswith("#") or not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 10:
                    counts["malformed_records"] += 1
                    continue
                counts["input_records"] += 1
                chrom_raw, pos, vid, ref, alt_field, _, _, _, fmt, sample = fields[:10]
                chrom = normalize_chrom(chrom_raw)
                if chrom is None:
                    counts["skipped_noncanonical_chromosome"] += 1
                    continue
                alts = alt_field.split(",")
                if len(ref) != 1 or any(len(a) != 1 for a in alts):
                    counts["skipped_non_snp"] += 1
                    continue
                ids = rsids(vid)
                if not ids:
                    counts["skipped_without_rsid"] += 1
                    continue
                indices = gt_indices(parse_gt(fmt, sample))
                if indices is None:
                    counts["skipped_missing_genotype"] += 1
                    continue
                genotype = genotype_bases(ref.upper(), [a.upper() for a in alts], indices)
                if genotype is None:
                    counts["skipped_invalid_genotype"] += 1
                    continue
                for rid in ids:
                    out.write(f"{rid}\t{chrom}\t{pos}\t{genotype}\n")
                    counts["output_rsid_rows"] += 1
                if all(x == 0 for x in indices):
                    counts["output_homozygous_reference_records"] += 1
                elif any(x > 0 for x in indices):
                    counts["output_alt_carrier_records"] += 1
        os.replace(part, output)
    except Exception:
        part.unlink(missing_ok=True)
        raise

    if counts["output_rsid_rows"] == 0:
        output.unlink(missing_ok=True)
        raise RuntimeError("No called rsID SNP genotypes were exported")

    write_report(
        report,
        "DANTE GENOTYPING-ARRAY 23ANDME-STYLE EXPORT",
        counts,
        [
            f"Source: {source}",
            "Interpretation boundary: homozygous-reference calls are meaningful only at explicitly interrogated array loci.",
            "Absence of a locus from this export is UNKNOWN, not homozygous reference.",
        ],
    )
    return 0


def normalized_header(lines: list[str], sample_header: str) -> list[str]:
    """Return a bcftools-safe canonical GRCh37 VCF header.

    Source array VCFs may either omit contig declarations or declare chr-prefixed
    names. Data rows are rewritten to 1..22/X/Y/MT, so the derived header must
    declare the same canonical IDs. We intentionally replace source contig lines
    with the fixed GRCh37 primary-assembly contract rather than guessing lengths.
    """
    out: list[str] = []
    inserted = False
    for line in lines:
        if line.startswith("##contig="):
            continue
        if line.startswith("#CHROM"):
            if not inserted:
                for contig, length in GRCH37_CONTIG_LENGTHS.items():
                    out.append(f"##contig=<ID={contig},length={length}>\n")
                out.append('##FILTER=<ID=PASS,Description="Accepted vendor genotype call for array adapter; source read-level FILTER unavailable">\n')
                out.append('##INFO=<ID=DANTE_ARRAY_SOURCE,Number=0,Type=Flag,Description="Record derived from genotyping-array VCF">\n')
                out.append('##INFO=<ID=DANTE_ORIGINAL_FILTER,Number=1,Type=String,Description="Original source FILTER field">\n')
                out.append('##INFO=<ID=DANTE_ARRAY_PHASE_STRIPPED,Number=0,Type=Flag,Description="Source GT used | without a documented phase-set contract; derived analytical GT is unphased">\n')
                inserted = True
            out.append(sample_header)
        else:
            out.append(line)
    return out


def _open_layer_outputs(parts: list[Path]) -> list[tuple[TextIO, gzip.GzipFile, object]]:
    handles: list[tuple[TextIO, gzip.GzipFile, object]] = []
    for part in parts:
        rawh = part.open("wb")
        gzh = gzip.GzipFile(filename="", mode="wb", fileobj=rawh, mtime=0)
        txth = io.TextIOWrapper(gzh, encoding="utf-8", newline="")
        handles.append((txth, gzh, rawh))
    return handles


def _write_layer_headers(handles: list[tuple[TextIO, gzip.GzipFile, object]], header_lines: list[str], sample_header: str) -> None:
    hdr = normalized_header(header_lines, sample_header)
    for txth, _, _ in handles:
        for line in hdr:
            txth.write(line)


def prepare_layers(source: Path, snp_out: Path, indel_out: Path, residual_out: Path, report: Path) -> int:
    source = source.resolve()
    paths = [snp_out.resolve(), indel_out.resolve(), residual_out.resolve()]
    for x in paths:
        x.parent.mkdir(parents=True, exist_ok=True)
    report = report.resolve()
    report.parent.mkdir(parents=True, exist_ok=True)

    header_lines: list[str] = []
    sample_header = ""
    counts: Counter[str] = Counter()
    parts = [x.with_name(x.name + ".part") for x in paths]
    handles: list[tuple[TextIO, gzip.GzipFile, object]] = []

    try:
        with open_text(source) as inp:
            for raw in inp:
                if raw.startswith("#"):
                    header_lines.append(raw if raw.endswith("\n") else raw + "\n")
                    if raw.startswith("#CHROM"):
                        cols = raw.rstrip("\r\n").split("\t")
                        if len(cols) != 10:
                            raise RuntimeError(f"Expected exactly one sample column, found {max(0, len(cols)-9)}")
                        sample_header = raw if raw.endswith("\n") else raw + "\n"
                    continue

                if not sample_header:
                    raise RuntimeError("VCF #CHROM header was not found before records")

                if not handles:
                    handles = _open_layer_outputs(parts)
                    _write_layer_headers(handles, header_lines, sample_header)

                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 10:
                    counts["malformed_records"] += 1
                    continue
                counts["input_records"] += 1
                chrom_raw, _, _, ref, alt_field, qual, filt, info, fmt, sample = fields
                chrom = normalize_chrom(chrom_raw)
                if chrom is None:
                    counts["skipped_noncanonical_chromosome"] += 1
                    continue
                original_gt = parse_gt(fmt, sample)
                indices = gt_indices(original_gt)
                if indices is None:
                    counts["skipped_missing_genotype"] += 1
                    continue
                if not any(x > 0 for x in indices):
                    counts["skipped_homozygous_reference"] += 1
                    continue

                alts = alt_field.split(",")
                if any(i >= len([ref, *alts]) for i in indices):
                    counts["skipped_invalid_genotype"] += 1
                    continue
                if len(ref) == 1 and alts and all(len(a) == 1 and not a.startswith("<") for a in alts):
                    target_index = 0
                    counts["carried_snp_records"] += 1
                elif alts and all(not a.startswith("<") and "[" not in a and "]" not in a for a in alts):
                    target_index = 1
                    counts["carried_indel_records"] += 1
                else:
                    target_index = 2
                    counts["carried_residual_records"] += 1

                new_info: list[str] = []
                if info and info != ".":
                    new_info.append(info)
                if "|" in original_gt:
                    new_info.append("DANTE_ARRAY_PHASE_STRIPPED")
                    fmt_keys = fmt.split(":")
                    sample_vals = sample.split(":")
                    if "GT" in fmt_keys:
                        gt_i = fmt_keys.index("GT")
                        if gt_i < len(sample_vals):
                            sample_vals[gt_i] = sample_vals[gt_i].replace("|", "/")
                            fields[9] = ":".join(sample_vals)
                            counts["phase_separators_stripped"] += 1
                new_info.extend(["DANTE_ARRAY_SOURCE", f"DANTE_ORIGINAL_FILTER={filt or '.'}"])
                fields[0] = chrom
                fields[5] = qual or "."
                fields[6] = "PASS"
                fields[7] = ";".join(new_info)
                handles[target_index][0].write("\t".join(fields) + "\n")

            if not sample_header:
                raise RuntimeError("VCF #CHROM header was not found")
            if not handles:
                handles = _open_layer_outputs(parts)
                _write_layer_headers(handles, header_lines, sample_header)

        for txth, gzh, rawh in handles:
            txth.close()
            if not gzh.closed:
                gzh.close()
            if not rawh.closed:
                rawh.close()
        handles = []
        for part, target in zip(parts, paths):
            os.replace(part, target)
    except Exception:
        for txth, gzh, rawh in handles:
            try:
                txth.close()
            except Exception:
                pass
            try:
                gzh.close()
            except Exception:
                pass
            try:
                rawh.close()
            except Exception:
                pass
        for part in parts:
            part.unlink(missing_ok=True)
        raise

    write_report(
        report,
        "DANTE GENOTYPING-ARRAY CARRIED-VARIANT PREPARATION",
        counts,
        [
            f"Source: {source}",
            "Only genotypes carrying ALT enter SNP/INDEL/residual annotation layers.",
            "0/0 loci are deliberately excluded from VEP and Stage 24; they remain available only in the direct array genotype export.",
            "Derived FILTER=PASS means accepted vendor genotype call, not read-level sequencing QC.",
            "Source | genotype separators are converted to / unless a phase-set contract exists; this adapter does not infer haplotypes.",
            "Derived VCF headers declare canonical GRCh37 primary contigs so bcftools/htslib can safely encode normalized output.",
        ],
    )
    return 0


def count_data_records(path: Path) -> int:
    count = 0
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("#") and line.strip():
                count += 1
    return count


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)

    exp = sub.add_parser("export-23andme")
    exp.add_argument("input", type=Path)
    exp.add_argument("--output", type=Path, required=True)
    exp.add_argument("--report", type=Path, required=True)

    prep = sub.add_parser("prepare-layers")
    prep.add_argument("input", type=Path)
    prep.add_argument("--snp-out", type=Path, required=True)
    prep.add_argument("--indel-out", type=Path, required=True)
    prep.add_argument("--residual-out", type=Path, required=True)
    prep.add_argument("--report", type=Path, required=True)

    sub.add_parser("self-test")
    args = ap.parse_args()

    if args.mode == "export-23andme":
        return export_23andme(args.input, args.output, args.report)
    if args.mode == "prepare-layers":
        return prepare_layers(args.input, args.snp_out, args.indel_out, args.residual_out, args.report)

    import tempfile
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        source = root / "array.vcf"
        source.write_text(
            "##fileformat=VCFv4.2\n##source=Genotek\n##reference=hg19\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n"
            "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0|0\n"
            "chr1\t200\trs2\tC\tT\t.\t.\t.\tGT\t0|1\n"
            "chr1\t300\trs3\tG\tA\t.\t.\t.\tGT\t./.\n"
            "chr2\t400\trs4\tA\tAT\t.\t.\t.\tGT\t1|1\n",
            encoding="utf-8",
        )
        out = root / "out.txt.gz"
        report = root / "export.txt"
        export_23andme(source, out, report)
        with gzip.open(out, "rt", encoding="utf-8") as handle:
            data = [x.strip() for x in handle if not x.startswith("#")]
        assert data == ["rs1\t1\t100\tAA", "rs2\t1\t200\tCT"], data

        snp = root / "snp.vcf.gz"
        indel = root / "indel.vcf.gz"
        residual = root / "residual.vcf.gz"
        prepare_layers(source, snp, indel, residual, root / "prep.txt")
        assert count_data_records(snp) == 1
        assert count_data_records(indel) == 1
        assert count_data_records(residual) == 0
        with gzip.open(snp, "rt", encoding="utf-8") as handle:
            text = handle.read()
        assert "##contig=<ID=1,length=249250621>\n" in text
        assert "##contig=<ID=2,length=243199373>\n" in text
        assert "##contig=<ID=MT,length=16569>\n" in text
        assert "##contig=<ID=chr1" not in text
        assert "\tPASS\t" in text
        assert "DANTE_ORIGINAL_FILTER=." in text
        assert "\t0|0\n" not in text
        assert "\tGT\t0/1\n" in text
    print("SELF-TEST PASS | genotyping-array export + carried-variant layers + canonical GRCh37 contig headers")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())