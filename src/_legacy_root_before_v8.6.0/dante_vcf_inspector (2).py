#!/usr/bin/env python3
"""Inspect Dante Labs VCF files without bcftools or third-party packages.

The script reads plain .vcf and compressed .vcf.gz files directly and writes
two files next to every input:

    <input_name>_header.txt
    <input_name>_stats.txt

Examples (Windows PowerShell):

    py dante_vcf_inspector.py
    py dante_vcf_inspector.py GFX0242962.raw.vcf.gz
    py dante_vcf_inspector.py --all

When the script is launched without arguments (including by double-clicking),
it shows a numbered menu of every .vcf and .vcf.gz file next to the script.
On Windows, a single digit selects the file immediately when the menu has no
more than nine entries. On an error, the console remains open until a key is
pressed.

The original VCF files are never modified.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import time
import traceback
from collections import Counter
from pathlib import Path
from typing import TextIO


RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)


def open_vcf(path: Path) -> TextIO:
    """Open a plain or gzip/bgzip-compressed VCF as text."""
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open(mode="rt", encoding="utf-8", errors="replace")


def output_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".vcf.gz"):
        return name[:-7]
    if name.lower().endswith(".vcf"):
        return name[:-4]
    return path.stem


def classify_alleles(ref: str, alt_field: str) -> set[str]:
    """Return one or more broad variant classes for a VCF record."""
    classes: set[str] = set()
    if alt_field in {"", "."}:
        return {"non-variant / ALT missing"}

    for alt in alt_field.split(","):
        if alt in {"<NON_REF>", "<*>", "*"}:
            classes.add("gVCF reference/spanning allele")
        elif alt.startswith("<") and alt.endswith(">"):
            classes.add("symbolic/SV")
        elif "[" in alt or "]" in alt:
            classes.add("breakend/SV")
        elif len(ref) == 1 and len(alt) == 1:
            classes.add("SNP")
        elif len(ref) == len(alt):
            classes.add("MNP/substitution")
        elif len(ref) < len(alt):
            classes.add("insertion")
        elif len(ref) > len(alt):
            classes.add("deletion")
        else:
            classes.add("other")
    return classes


def classify_gt(gt: str) -> str:
    """Classify a GT value while accepting phased and unphased genotypes."""
    gt = gt.strip()
    if not gt or gt in {".", "./.", ".|."}:
        return "missing"

    alleles = re.split(r"[/|]", gt)
    if not alleles or any(a == "." for a in alleles):
        return "partially missing"
    if not all(a.isdigit() for a in alleles):
        return "other/non-numeric"

    numbers = [int(a) for a in alleles]
    if all(a == 0 for a in numbers):
        return "hom_ref (0/0 or 0)"
    if len(set(numbers)) == 1 and numbers[0] > 0:
        return "hom_alt"
    if 0 in numbers and any(a > 0 for a in numbers):
        return "het ref/alt"
    return "het alt/alt or other"


def infer_build(header_lines: list[str]) -> str:
    joined = "\n".join(header_lines).lower()
    if any(token in joined for token in ("grch38", "hg38", "b38")):
        return "GRCh38/hg38 (header hint)"
    if any(token in joined for token in ("grch37", "hg19", "b37")):
        return "GRCh37/hg19 (header hint)"
    return "Not determined from header"


def inspect_vcf(path: Path, progress_every: int = 1_000_000) -> tuple[Path, Path]:
    started = time.time()
    header_lines: list[str] = []
    samples: list[str] = []
    column_header = ""

    record_count = 0
    malformed_count = 0
    known_id_records = 0
    rsid_records = 0
    biallelic_records = 0
    multiallelic_records = 0
    info_end_records = 0
    gvcf_alt_records = 0
    all_samples_hom_ref_records = 0

    chromosome_counts: Counter[str] = Counter()
    filter_counts: Counter[str] = Counter()
    variant_counts: Counter[str] = Counter()
    genotype_counts: Counter[str] = Counter()
    format_field_counts: Counter[str] = Counter()

    with open_vcf(path) as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                line = raw_line.rstrip("\r\n")
                header_lines.append(line)
                if line.startswith("#CHROM"):
                    column_header = line
                    columns = line.split("\t")
                    samples = columns[9:] if len(columns) > 9 else []
                continue

            if not raw_line.strip():
                continue

            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                malformed_count += 1
                continue

            record_count += 1
            chrom, _pos, id_field, ref, alt, _qual, filt, info = fields[:8]
            chromosome_counts[chrom] += 1

            ids = [value for value in id_field.split(";") if value and value != "."]
            if ids:
                known_id_records += 1
            if any(RSID_RE.fullmatch(value) for value in ids):
                rsid_records += 1

            alt_values = alt.split(",") if alt not in {"", "."} else []
            if len(alt_values) <= 1:
                biallelic_records += 1
            else:
                multiallelic_records += 1

            for variant_class in classify_alleles(ref, alt):
                variant_counts[variant_class] += 1

            if "<NON_REF>" in alt_values or "<*>" in alt_values:
                gvcf_alt_records += 1
            if re.search(r"(?:^|;)END=\d+(?:;|$)", info):
                info_end_records += 1

            for value in filt.split(";"):
                filter_counts[value or "."] += 1

            record_gt_classes: list[str] = []
            if len(fields) >= 9:
                format_fields = fields[8].split(":")
                for format_field in format_fields:
                    if format_field:
                        format_field_counts[format_field] += 1

                if "GT" in format_fields:
                    gt_index = format_fields.index("GT")
                    for sample_field in fields[9:]:
                        values = sample_field.split(":")
                        gt = values[gt_index] if gt_index < len(values) else "."
                        gt_class = classify_gt(gt)
                        genotype_counts[gt_class] += 1
                        record_gt_classes.append(gt_class)

            if record_gt_classes and all(
                value == "hom_ref (0/0 or 0)" for value in record_gt_classes
            ):
                all_samples_hom_ref_records += 1

            if progress_every and record_count % progress_every == 0:
                elapsed = time.time() - started
                print(
                    f"  {path.name}: read {record_count:,} records "
                    f"({elapsed:,.1f} s)",
                    flush=True,
                )

    stem = output_stem(path)
    header_path = path.with_name(f"{stem}_header.txt")
    stats_path = path.with_name(f"{stem}_stats.txt")

    header_path.write_text("\n".join(header_lines) + "\n", encoding="utf-8")

    def percent(value: int, total: int) -> str:
        return f"{(100.0 * value / total):.3f}%" if total else "n/a"

    likely_gvcf = bool(gvcf_alt_records or info_end_records)
    likely_all_sites = record_count > 0 and (
        all_samples_hom_ref_records / record_count >= 0.05
    )

    report: list[str] = [
        "DANTE VCF INSPECTION REPORT",
        "=" * 30,
        f"Input file: {path.name}",
        f"Absolute path: {path.resolve()}",
        f"File size on disk: {path.stat().st_size:,} bytes",
        f"Compression: {'gzip/bgzip' if path.name.lower().endswith('.gz') else 'none'}",
        f"Elapsed time: {time.time() - started:,.1f} seconds",
        "",
        "HEADER SUMMARY",
        "-" * 30,
        f"Header lines: {len(header_lines):,}",
        f"#CHROM line present: {'yes' if column_header else 'no'}",
        f"Samples ({len(samples)}): {', '.join(samples) if samples else '(none)'}",
        f"Genome build: {infer_build(header_lines)}",
        "",
        "RECORD SUMMARY",
        "-" * 30,
        f"Total records: {record_count:,}",
        f"Malformed/short records skipped: {malformed_count:,}",
        f"Records with any ID: {known_id_records:,} ({percent(known_id_records, record_count)})",
        f"Records with an rsID: {rsid_records:,} ({percent(rsid_records, record_count)})",
        f"Biallelic or ALT missing: {biallelic_records:,}",
        f"Multiallelic: {multiallelic_records:,}",
        f"Records with <NON_REF> or <*>: {gvcf_alt_records:,}",
        f"Records with INFO/END: {info_end_records:,}",
        f"Records where all samples are hom_ref: {all_samples_hom_ref_records:,} "
        f"({percent(all_samples_hom_ref_records, record_count)})",
        "",
        "PRELIMINARY FILE-TYPE INFERENCE",
        "-" * 30,
        f"gVCF/reference-block signs: {'present' if likely_gvcf else 'not detected'}",
        f"Substantial reference-call content: {'possible' if likely_all_sites else 'not detected'}",
        (
            "Interpretation: this may be a gVCF or contain reference blocks/calls."
            if likely_gvcf or likely_all_sites
            else "Interpretation: this looks more like a variant-only VCF than an all-sites genome."
        ),
        "This is a structural inference, not a clinical interpretation.",
        "",
        "VARIANT CLASSES",
        "-" * 30,
    ]

    report.extend(f"{key}: {value:,}" for key, value in variant_counts.most_common())
    report.extend(["", "FILTER VALUES", "-" * 30])
    report.extend(f"{key}: {value:,}" for key, value in filter_counts.most_common())
    report.extend(["", "GENOTYPE CALLS (ALL SAMPLES)", "-" * 30])
    if genotype_counts:
        report.extend(f"{key}: {value:,}" for key, value in genotype_counts.most_common())
    else:
        report.append("No GT calls found")

    report.extend(["", "FORMAT FIELDS", "-" * 30])
    if format_field_counts:
        report.extend(
            f"{key}: present in {value:,} records"
            for key, value in format_field_counts.most_common()
        )
    else:
        report.append("No FORMAT fields found")

    report.extend(["", "RECORDS BY CHROMOSOME/CONTIG", "-" * 30])
    report.extend(f"{key}: {value:,}" for key, value in chromosome_counts.items())
    report.append("")

    stats_path.write_text("\n".join(report), encoding="utf-8")
    return header_path, stats_path


def discover_vcfs(directory: Path) -> list[Path]:
    """Find supported VCF files directly beside the script."""
    paths = [
        path.resolve()
        for path in directory.iterdir()
        if path.is_file()
        and (
            path.name.lower().endswith(".vcf")
            or path.name.lower().endswith(".vcf.gz")
        )
    ]
    return sorted(set(paths), key=lambda path: path.name.lower())


def human_file_size(size: int) -> str:
    """Return a compact binary file size for the interactive menu."""
    value = float(size)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.2f} {unit}"
        value /= 1024.0
    return f"{size} B"


def vcf_role_hint(path: Path) -> str:
    """Give a filename-based hint without interpreting the file contents."""
    name = path.name.lower()
    if ".filtered.snp." in name or name.endswith(".filtered.snp.vcf"):
        return "filtered SNP"
    if ".filtered.indel." in name or name.endswith(".filtered.indel.vcf"):
        return "filtered short indels"
    if ".cnv." in name or name.endswith(".cnv.vcf"):
        return "copy-number variants (CNV)"
    if ".sv." in name or name.endswith(".sv.vcf"):
        return "structural variants (SV)"
    if ".raw." in name or name.endswith(".raw.vcf"):
        return "raw small variants"
    return "VCF"


def read_menu_choice(number_of_files: int) -> int | None:
    """Read a menu selection; return a zero-based index or None to quit."""
    if number_of_files < 1:
        return None
    if not sys.stdin or not sys.stdin.isatty():
        raise RuntimeError(
            "Interactive selection requires a console. Pass a VCF filename "
            "on the command line instead."
        )

    prompt = f"Select a file [1-{number_of_files}] or 0 to exit: "

    # With the expected five Dante files, Windows can select immediately by
    # one key. Larger menus use input() so multi-digit numbers remain possible.
    if os.name == "nt" and number_of_files <= 9:
        try:
            import msvcrt

            valid = {str(index): index - 1 for index in range(1, number_of_files + 1)}
            print(prompt, end="", flush=True)
            while True:
                key = msvcrt.getwch()
                if key in valid:
                    print(key)
                    return valid[key]
                if key in {"0", "q", "Q", "\x1b"}:
                    print(key if key != "\x1b" else "Esc")
                    return None
        except (ImportError, OSError):
            pass

    while True:
        choice = input(prompt).strip()
        if choice.lower() in {"0", "q", "quit", "exit"}:
            return None
        if choice.isdigit() and 1 <= int(choice) <= number_of_files:
            return int(choice) - 1
        print(f"Please enter a number from 1 to {number_of_files}, or 0 to exit.")


def choose_vcf_interactively(directory: Path) -> Path | None:
    """Show every discovered VCF and return the file chosen by the user."""
    paths = discover_vcfs(directory)
    if not paths:
        return None

    print("Dante VCF Inspector")
    print("VCF files found next to this script:")
    for index, path in enumerate(paths, start=1):
        print(
            f"  {index}. {path.name} "
            f"[{human_file_size(path.stat().st_size)}; {vcf_role_hint(path)}]"
        )
    print()

    selected_index = read_menu_choice(len(paths))
    return paths[selected_index] if selected_index is not None else None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create header and diagnostic statistics TXT files from .vcf or .vcf.gz "
            "without unpacking the source files."
        )
    )
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="VCF/VCF.GZ files to inspect",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Inspect every .vcf and .vcf.gz in the script directory",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Print progress every N records; use 0 to disable (default: 1000000)",
    )
    return parser.parse_args()


def pause_after_error() -> None:
    """Keep an interactive console open so the user can read the error."""
    if not sys.stdin or not sys.stdin.isatty():
        return

    print("\nPress any key to close this window...", file=sys.stderr, flush=True)
    if os.name == "nt":
        try:
            import msvcrt

            msvcrt.getwch()
            return
        except (ImportError, OSError):
            pass

    try:
        input()
    except (EOFError, KeyboardInterrupt):
        pass


def main() -> int:
    args = parse_args()
    script_directory = Path(__file__).resolve().parent

    if args.all:
        paths = discover_vcfs(script_directory)
    elif not args.files:
        discovered_paths = discover_vcfs(script_directory)
        if not discovered_paths:
            paths = []
        else:
            selected_path = choose_vcf_interactively(script_directory)
            if selected_path is None:
                print("Cancelled. No files were changed.")
                return 0
            paths = [selected_path]
    else:
        paths = [path.resolve() for path in args.files]

    if not paths:
        print("No input files were selected.")
        print("Put one or more .vcf or .vcf.gz files next to this script and launch it again.")
        print("You can also pass one file explicitly, for example:")
        print('  py dante_vcf_inspector.py "GFX0242962.raw.vcf.gz"')
        print("To inspect every VCF in the folder, run:")
        print("  py dante_vcf_inspector.py --all")
        return 2

    missing = [path for path in paths if not path.is_file()]
    if missing:
        for path in missing:
            print(f"File not found: {path}", file=sys.stderr)
        return 2

    failures = 0
    for path in paths:
        print(f"Inspecting: {path.name}", flush=True)
        try:
            header_path, stats_path = inspect_vcf(
                path, progress_every=max(0, args.progress_every)
            )
            print(f"  Created: {header_path.name}")
            print(f"  Created: {stats_path.name}")
        except Exception as exc:
            failures += 1
            print(
                f"  ERROR while reading {path.name}: "
                f"{type(exc).__name__}: {exc}",
                file=sys.stderr,
            )
            traceback.print_exc()

    if failures:
        print(f"Finished with errors in {failures} file(s).", file=sys.stderr)
        return 1

    print("Done. The original VCF files were not changed.")
    return 0


def run_cli() -> None:
    """Run the command and pause before every unsuccessful exit."""
    try:
        exit_code = main()
    except KeyboardInterrupt:
        print("\nERROR: operation cancelled by the user.", file=sys.stderr)
        exit_code = 130
    except SystemExit as exc:
        # argparse uses SystemExit for both --help (0) and syntax errors (non-zero).
        exit_code = exc.code if isinstance(exc.code, int) else 1
    except BaseException as exc:
        print(
            f"\nUNEXPECTED ERROR: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        traceback.print_exc()
        exit_code = 1

    if exit_code:
        pause_after_error()
    raise SystemExit(exit_code)


if __name__ == "__main__":
    run_cli()
