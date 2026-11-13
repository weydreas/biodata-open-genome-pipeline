from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import fastq_pairing  # noqa: E402
import reference_alignment  # noqa: E402
import reference_calling  # noqa: E402
import telomere_estimate  # noqa: E402
from reference_genotype_export import RawCall, evaluate_call  # noqa: E402


def write_fastq(path: Path, read: int, names: tuple[str, ...] = ("read1", "read2")) -> None:
    handle_context = (
        gzip.open(path, "wt", encoding="ascii", newline="\n")
        if path.name.endswith(".gz")
        else path.open("wt", encoding="ascii", newline="\n")
    )
    with handle_context as handle:
        for name in names:
            handle.write(f"@{name} {read}:N:0:1\nACGT\n+\nIIII\n")


class FastqPairingTests(unittest.TestCase):
    def test_zero_files_is_missing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(fastq_pairing.discover(Path(directory))["status"], "MISSING")

    def test_r1_without_r2_is_single_end_ready(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fastq(root / "Dante_R1_001.fastq.gz", 1)
            result = fastq_pairing.discover(root, True)
            self.assertEqual(result["status"], "SINGLE_END_READY")
            self.assertEqual(result["read_mode"], "single")
            self.assertIsNotNone(result["r1"])
            self.assertRegex(str(result["r1_content_sha256"]), r"^[0-9a-f]{64}$")
            self.assertEqual(result["r1_record_count"], 2)
            self.assertGreater(int(result["r1_decompressed_bytes"]), 0)

    def test_uncompressed_single_end_r1_full_check(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fastq(root / "Dante_R1.fastq", 1)
            result = fastq_pairing.discover(root, True)
            self.assertEqual(result["status"], "SINGLE_END_READY")
            self.assertRegex(str(result["r1_content_sha256"]), r"^[0-9a-f]{64}$")

    def test_one_valid_pair_is_ready(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fastq(root / "Dante_R1_001.fastq.gz", 1)
            write_fastq(root / "Dante_R2_001.fastq.gz", 2)
            result = fastq_pairing.discover(root, True)
            self.assertEqual(result["status"], "READY")
            self.assertEqual(result["read_mode"], "paired")
            self.assertEqual(result["r1_record_count"], 2)
            self.assertEqual(result["r2_record_count"], 2)

    def test_require_complete_accepts_single_end_r1(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fastq(root / "Dante_R1.fastq.gz", 1)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SRC / "fastq_pairing.py"),
                    "--root",
                    str(root),
                    "--require-complete",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_adding_r2_changes_mode_from_single_to_paired(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fastq(root / "Dante_R1.fastq.gz", 1)
            first = fastq_pairing.discover(root)
            self.assertEqual((first["status"], first["read_mode"]), ("SINGLE_END_READY", "single"))
            write_fastq(root / "Dante_R2.fastq.gz", 2)
            second = fastq_pairing.discover(root)
            self.assertEqual((second["status"], second["read_mode"]), ("READY", "paired"))

    def test_mate_name_mismatch_is_invalid(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fastq(root / "Dante_R1.fastq.gz", 1, ("one",))
            write_fastq(root / "Dante_R2.fastq.gz", 2, ("other",))
            result = fastq_pairing.discover(root)
            self.assertEqual(result["status"], "INVALID")
            self.assertIn("mate names differ", " ".join(result["details"]))

    def test_multiple_pairs_are_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for sample in ("A", "B"):
                write_fastq(root / f"{sample}_R1.fastq.gz", 1)
                write_fastq(root / f"{sample}_R2.fastq.gz", 2)
            self.assertEqual(fastq_pairing.discover(root)["status"], "MULTIPLE_PAIRS")

    def test_corrupt_gzip_fails_full_check(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fastq(root / "Dante_R1.fastq.gz", 1)
            write_fastq(root / "Dante_R2.fastq.gz", 2)
            r2 = root / "Dante_R2.fastq.gz"
            r2.write_bytes(r2.read_bytes()[:-5])
            self.assertEqual(fastq_pairing.discover(root, True)["status"], "INVALID")

    def test_paired_full_check_rejects_different_record_counts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_fastq(root / "Dante_R1.fastq.gz", 1, ("one", "two"))
            write_fastq(root / "Dante_R2.fastq.gz", 2, ("one",))
            result = fastq_pairing.discover(root, True)
            self.assertEqual(result["status"], "INVALID")
            self.assertIn("record counts differ", " ".join(result["details"]))


class AlignmentProgressTests(unittest.TestCase):
    def test_bwa_completed_messages_drive_exact_read_progress(self) -> None:
        tracker = reference_alignment.AlignmentProgress(400_000)
        with mock.patch.object(reference_alignment, "stage15_progress") as emit:
            tracker.handle_line(b"[M::mem_process_seqs] Processed 200000 reads in 12.0 CPU sec, 3.0 real sec\n")
            self.assertEqual(tracker.processed, 200_000)
            self.assertTrue(emit.called)
            percent = emit.call_args.args[0]
            self.assertAlmostEqual(percent, 49.0, places=3)

    def test_bwa_process_and_completed_pair_is_counted_once(self) -> None:
        tracker = reference_alignment.AlignmentProgress(400_000)
        with mock.patch.object(reference_alignment, "stage15_progress"):
            tracker.handle_line(b"[M::process] read 200000 sequences (30000000 bp)...\n")
            self.assertEqual(tracker.processed, 0)
            tracker.handle_line(b"[M::mem_process_seqs] Processed 200000 reads in 12.0 CPU sec, 3.0 real sec\n")
            self.assertEqual(tracker.processed, 200_000)

    def test_non_process_bwa_line_does_not_advance_progress(self) -> None:
        tracker = reference_alignment.AlignmentProgress(100)
        with mock.patch.object(reference_alignment, "stage15_progress") as emit:
            tracker.handle_line(b"[M::mem_pestat] mean and std.dev: (500.0, 50.0)\n")
            self.assertEqual(tracker.processed, 0)
            emit.assert_not_called()



class TargetPanelTests(unittest.TestCase):
    def test_clinvar_and_optional_panel_are_deduplicated(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            clinvar = root / "clinvar.vcf.gz"
            with gzip.open(clinvar, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                handle.write("1\t100\t.\tA\tG\t.\tPASS\tRS=1\n")
                handle.write("1\t150\t.\tAA\tA\t.\tPASS\tRS=99\n")
            panel_dir = root / "panels"
            panel_dir.mkdir()
            (panel_dir / "chip.txt").write_text(
                "# rsid chromosome position genotype\nrs1\t1\t100\tAG\nrs2\t1\t200\tCC\nrs3\t2\t300\tAA\n",
                encoding="utf-8",
            )
            (panel_dir / "conflict.txt").write_text("rs3\t2\t301\tAA\n", encoding="utf-8")
            targets = root / "targets.tsv"
            bed = root / "targets.bed"
            report = root / "report.txt"
            coordinate_conflicts = root / "coordinate_conflicts.tsv"
            genotype_conflicts = root / "genotype_conflicts.tsv"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SRC / "build_reference_targets.py"),
                    "--clinvar",
                    str(clinvar),
                    "--panel-dir",
                    str(panel_dir),
                    "--targets",
                    str(targets),
                    "--bed",
                    str(bed),
                    "--report",
                    str(report),
                    "--coordinate-conflicts",
                    str(coordinate_conflicts),
                    "--genotype-conflicts",
                    str(genotype_conflicts),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            text = targets.read_text(encoding="utf-8")
            self.assertEqual(text.count("\n"), 3)  # header + rs1 + rs2
            self.assertIn("rs1\t1\t100\tA\tG\t", text)
            self.assertIn("rs2\t1\t200\t.\t.\tchip.txt\tCC\tAVAILABLE\tchip.txt", text)
            self.assertIn("Conflicting-coordinate rsIDs excluded: 1", report.read_text(encoding="utf-8"))
            conflict_text = coordinate_conflicts.read_text(encoding="utf-8")
            self.assertIn("rs3\t2\t300\tchip.txt", conflict_text)
            self.assertIn("rs3\t2\t301\tconflict.txt", conflict_text)


    def test_conflicting_reference_source_genotypes_are_not_available_for_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            clinvar = root / "clinvar.vcf.gz"
            with gzip.open(clinvar, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                handle.write("1\t100\t.\tA\tG\t.\tPASS\tRS=10\n")
            panel_dir = root / "panels"
            panel_dir.mkdir()
            (panel_dir / "23andMe.txt").write_text("rs10\t1\t100\tAG\n", encoding="utf-8")
            (panel_dir / "Genotek.txt").write_text("rs10\t1\t100\tAA\n", encoding="utf-8")
            targets = root / "targets.tsv"
            bed = root / "targets.bed"
            report = root / "report.txt"
            coordinate_conflicts = root / "coordinate_conflicts.tsv"
            genotype_conflicts = root / "genotype_conflicts.tsv"
            completed = subprocess.run([
                sys.executable, str(SRC / "build_reference_targets.py"),
                "--clinvar", str(clinvar), "--panel-dir", str(panel_dir),
                "--targets", str(targets), "--bed", str(bed), "--report", str(report),
                "--coordinate-conflicts", str(coordinate_conflicts),
                "--genotype-conflicts", str(genotype_conflicts),
            ], capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            text = targets.read_text(encoding="utf-8")
            self.assertIn("rs10\t1\t100\tA\tG\t23andMe.txt,ClinVar,Genotek.txt\t--\tCONFLICT", text)
            conflicts = genotype_conflicts.read_text(encoding="utf-8")
            self.assertIn("rs10\t1\t100\tAA\tGenotek.txt", conflicts)
            self.assertIn("rs10\t1\t100\tAG\t23andMe.txt", conflicts)


class ReferenceExportTests(unittest.TestCase):
    def test_threshold_boundaries(self) -> None:
        high_hom = RawCall("1", 1, "A", ("G",), {"GT": "0/0", "AD": "9,1", "DP": "10", "GQ": "30"}, {})
        self.assertEqual(evaluate_call(high_hom, "A", "G").tier, "A")
        high_het = RawCall("1", 1, "A", ("G",), {"GT": "0/1", "AD": "3,7", "DP": "10", "GQ": "30"}, {})
        self.assertEqual(evaluate_call(high_het, "A", "G").tier, "A")
        supportive_hom = RawCall("1", 1, "A", ("G",), {"GT": "0/0", "AD": "4,1", "DP": "5", "GQ": "10"}, {})
        self.assertEqual(evaluate_call(supportive_hom, "A", "G").tier, "B")
        supportive_het = RawCall("1", 1, "A", ("G",), {"GT": "0/1", "AD": "3,2", "DP": "5", "GQ": "10"}, {})
        self.assertEqual(evaluate_call(supportive_het, "A", "G").tier, "B")
        hypothesis = RawCall("1", 1, "A", ("G",), {"GT": "0/0", "AD": "3,2", "DP": "5", "GQ": "5"}, {})
        self.assertEqual(evaluate_call(hypothesis, "A", "G").tier, "C")
        refq = RawCall("1", 1, "A", ("G",), {"GT": "0/0", "AD": "5,0", "DP": "5", "PL": "0,15,160"}, {})
        evaluated = evaluate_call(refq, "A", "G")
        self.assertEqual(evaluated.tier, "B")
        self.assertEqual(evaluated.quality_metric, "REFQ_PL")

    def test_strict_calls_no_calls_merge_priority_and_anonymous_gzip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            targets = root / "targets.tsv"
            targets.write_text(
                "rsid\tchromosome\tposition\tref\talts\tsources\treference_genotype\treference_genotype_status\treference_genotype_sources\n"
                "rs1\t1\t100\tA\tG\tClinVar\t--\tMISSING\t.\n"
                "rs2\t1\t200\tC\tT\tchip\tCT\tAVAILABLE\t23andMe.txt\n"
                "rs3\t1\t300\tG\tA\tchip\tGG\tAVAILABLE\tGenotek.txt\n"
                "rs4\t1\t400\tT\tC\tchip\tTT\tAVAILABLE\tAtlas.txt\n",
                encoding="utf-8",
            )
            calls = root / "calls.vcf.gz"
            with gzip.open(calls, "wt", encoding="utf-8") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDante\n")
                handle.write("1\t100\t.\tA\tG\t.\tPASS\tDP=12\tGT:AD:DP:GQ\t0/0:12,0:12:40\n")
                handle.write("1\t200\t.\tC\tT\t.\tPASS\tDP=12\tGT:AD:DP:GQ\t0/1:6,6:12:50\n")
                handle.write("1\t300\t.\tG\tA\t.\tPASS\tDP=8\tGT:AD:DP:GQ\t0/0:8,0:8:50\n")
            base = root / "base.txt.gz"
            with gzip.open(base, "wt", encoding="utf-8") as handle:
                handle.write("# rsid chromosome position genotype\n")
                handle.write("rs2\t1\t200\tTT\n")
                handle.write("rs5\t1\t500\tAG\n")
            output = root / "reference.txt.gz"
            report = root / "report.txt"
            discrepancies = root / "discrepancies.tsv"
            provenance = root / "provenance.tsv"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SRC / "reference_genotype_export.py"),
                    "--targets",
                    str(targets),
                    "--calls",
                    str(calls),
                    "--base",
                    str(base),
                    "--output",
                    str(output),
                    "--report",
                    str(report),
                    "--discrepancies",
                    str(discrepancies),
                    "--provenance",
                    str(provenance),
                    "--fastq-source",
                    "R1",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                rows = [line.strip() for line in handle if not line.startswith("#")]
            self.assertEqual(
                rows,
                ["rs1\t1\t100\tAA", "rs2\t1\t200\tTT", "rs3\t1\t300\tGG", "rs4\t1\t400\tTT", "rs5\t1\t500\tAG"],
            )
            discrepancy_text = discrepancies.read_text(encoding="utf-8")
            self.assertIn("DANTE_VCF_VS_FASTQ", discrepancy_text)
            self.assertIn("DANTE_VCF_VS_REFERENCE_SOURCE", discrepancy_text)
            provenance_text = provenance.read_text(encoding="utf-8")
            self.assertIn("rs1\t1\t100\tAA\tR1", provenance_text)
            self.assertIn("rs3\t1\t300\tGG\tR1\tGG\tB\tTIER_B_SUPPORTIVE", provenance_text)
            self.assertIn("rs4\t1\t400\tTT\tAtlas.txt", provenance_text)
            self.assertIn("Explicit no-calls (--): 0", report.read_text(encoding="utf-8"))
            raw = output.read_bytes()
            flags = raw[3]
            self.assertEqual(flags & 0x08, 0, "gzip header unexpectedly embeds a .part filename")


    def test_reference_source_only_premerge_without_fastq_calls(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            targets = root / "targets.tsv"
            targets.write_text(
                "rsid\tchromosome\tposition\tref\talts\tsources\treference_genotype\treference_genotype_status\treference_genotype_sources\n"
                "rs1\t1\t100\tA\tG\tchip\tAA\tAVAILABLE\t23andMe.txt\n"
                "rs2\t1\t200\tC\tT\tchip\tCT\tAVAILABLE\tGenotek.txt,Atlas.txt\n"
                "rs3\t1\t300\tG\tA\tchip\t--\tCONFLICT\t23andMe.txt,Genotek.txt\n",
                encoding="utf-8",
            )
            base = root / "base.txt.gz"
            with gzip.open(base, "wt", encoding="utf-8") as handle:
                handle.write("# rsid chromosome position genotype\n")
                handle.write("rs1\t1\t100\tAG\n")
                handle.write("rs4\t1\t400\tTT\n")
            output = root / "reference.txt.gz"
            report = root / "report.txt"
            discrepancies = root / "discrepancies.tsv"
            provenance = root / "provenance.tsv"
            completed = subprocess.run(
                [
                    sys.executable, str(SRC / "reference_genotype_export.py"),
                    "--targets", str(targets), "--base", str(base),
                    "--output", str(output), "--report", str(report),
                    "--discrepancies", str(discrepancies), "--provenance", str(provenance),
                    "--fastq-source", "NONE",
                ], capture_output=True, text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                rows = [line.strip() for line in handle if not line.startswith("#")]
            self.assertEqual(rows, ["rs1\t1\t100\tAG", "rs2\t1\t200\tCT", "rs3\t1\t300\t--", "rs4\t1\t400\tTT"])
            self.assertIn("DANTE_VCF_VS_REFERENCE_SOURCE", discrepancies.read_text(encoding="utf-8"))
            self.assertIn("rs2\t1\t200\tCT\tGenotek.txt,Atlas.txt", provenance.read_text(encoding="utf-8"))
            self.assertIn("FASTQ evidence source: NONE", report.read_text(encoding="utf-8"))


class TelomereEstimateTests(unittest.TestCase):
    def test_telseq_formula_contract(self) -> None:
        estimate = telomere_estimate.estimate_kb(46_000, 332_720_800)
        self.assertIsNotNone(estimate)
        self.assertAlmostEqual(float(estimate), 1.0, places=9)

    def test_telomeric_motif_count_is_non_overlapping_and_bidirectional(self) -> None:
        self.assertEqual(telomere_estimate.telomeric_motif_count(b"TTAGGG" * 8), 8)
        self.assertEqual(telomere_estimate.telomeric_motif_count(b"CCCTAA" * 9), 9)
        self.assertEqual(telomere_estimate.telomeric_motif_count(b"ACGT" * 20), 0)

    def test_single_end_fastq_scan_counts_k7_and_gc_normalizer(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastq = root / "R1.fastq.gz"
            with gzip.open(fastq, "wt", encoding="ascii", newline="\n") as handle:
                seq1 = "TTAGGG" * 10
                seq2 = "ACGT" * 15
                for i, seq in enumerate((seq1, seq2), 1):
                    handle.write(f"@r{i} 1:N:0:1\n{seq}\n+\n{'I' * len(seq)}\n")
            metrics = telomere_estimate.scan_fastq(fastq, 2, 0, 2, "R1")
            self.assertEqual(metrics["records"], 2)
            self.assertEqual(metrics["telomeric_reads_k7"], 1)
            self.assertEqual(metrics["gc_normalizer_reads"], 2)
            self.assertIsNotNone(metrics["estimate_kb_k7"])

    def test_cached_metrics_allow_report_date_change_without_rescan(self) -> None:
        cached = {
            "schema": telomere_estimate.SCHEMA,
            "read_mode": "single",
            "r1": {"content_sha256": "a" * 64},
            "r2": None,
        }
        self.assertTrue(telomere_estimate.cached_matches(cached, "single", "a" * 64, ""))
        self.assertFalse(telomere_estimate.cached_matches(cached, "paired", "a" * 64, "b" * 64))


class OrchestratorStaticTests(unittest.TestCase):
    def test_new_stage_and_version_are_present(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn("Dante WGS local pipeline v9.0.0", script)
        self.assertIn('Write-MenuItem "14" "Merge known array genotypes', script)
        self.assertIn('Write-MenuItem "18" "Run complete stage 3', script)
        self.assertIn('Write-MenuItem "19" "Estimate mean telomere length', script)
        self.assertIn('Write-MenuItem "24" "Enrich full genomic outputs', script)
        self.assertIn('Write-MenuItem "25" "Research prioritized Stage 24 findings', script)
        self.assertIn('Write-MenuItem "26" "Launch local IGV.js browser', script)
        self.assertIn('Write-MenuItem "21" "Install and run stages 1 + 2', script)
        self.assertIn('Write-MenuItem "22" "Install/update open reference knowledge bundle', script)
        self.assertIn('"open-reference" { Install-OpenReferenceKnowledge }', script)
        self.assertIn('"stage3" { Run-Stage3 }', script)
        self.assertIn('placeholder or missing; only stage 6 / item 25 is affected', script)
        self.assertIn('"SINGLE_END_READY" { return "READY: SINGLE-END R1" }', script)
        self.assertIn('schema = "reference-alignment-v3"', script)
        self.assertIn('function Resolve-DockerImage', script)
        self.assertIn('Using local approved $Label image', script)
        self.assertIn('staphb/bwa:0.7.19', script)
        self.assertIn('staphb/samtools:1.22.1', script)
        self.assertIn('staphb/bcftools:1.24', script)
        self.assertNotIn('Invoke-NativeWithRetry docker @("pull", $BcftoolsImage)', script)

    def test_alignment_r2_argument_is_optional(self) -> None:
        source = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        self.assertIn('parser.add_argument("--r2", type=Path)', source)
        self.assertNotIn('parser.add_argument("--r2", type=Path, required=True)', source)

    def test_alignment_index_uses_tool_contract_not_registry_source(self) -> None:
        source = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        self.assertIn('parser.add_argument("--bwa-contract", required=True)', source)
        self.assertIn('state.get("bwa_contract") == bwa_contract', source)
        self.assertIn('"schema": "bwa-index-v2"', source)

    def test_reference_export_streams_calls_instead_of_loading_all_into_memory(self) -> None:
        source = (SRC / "reference_genotype_export.py").read_text(encoding="utf-8")
        self.assertIn("class CallStream", source)
        self.assertNotIn("def load_calls", source)
        self.assertIn("[OVERALL] Stage {self.stage}", source)

    def test_single_end_bwa_command_contains_only_r1(self) -> None:
        command = reference_alignment.build_bwa_command(
            ROOT,
            ROOT / "index",
            "bwa:test",
            4,
            "@RG\\tID:Dante\\tSM:Dante",
            "Dante_R1.fastq",
            None,
        )
        self.assertEqual(command[-1], "/checkpoint/Dante_R1.fastq")
        self.assertNotIn("/checkpoint/Dante_R2.fastq", command)

    def test_paired_end_bwa_command_contains_both_reads(self) -> None:
        command = reference_alignment.build_bwa_command(
            ROOT,
            ROOT / "index",
            "bwa:test",
            4,
            "@RG\\tID:Dante\\tSM:Dante",
            "Dante_R1.fastq",
            "Dante_R2.fastq",
        )
        self.assertEqual(command[-2:], ["/checkpoint/Dante_R1.fastq", "/checkpoint/Dante_R2.fastq"])


    def test_runtime_tool_invocations_override_container_entrypoints(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        alignment = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        calling = (SRC / "reference_calling.py").read_text(encoding="utf-8")
        self.assertIn('"--entrypoint", "bwa", $bwaRuntimeImage', script)
        self.assertIn('"--entrypoint", "samtools", $samtoolsRuntimeImage', script)
        self.assertIn('"--entrypoint", "bcftools", $bcftoolsRuntimeImage', script)
        self.assertGreaterEqual(alignment.count('"--entrypoint"'), 5)
        self.assertIn('        \"bwa\",\n        bwa_image,', alignment)
        self.assertIn('\"samtools\",\n            samtools_image,', alignment)
        self.assertIn('\"bash\",\n        samtools_image,', alignment)
        self.assertGreaterEqual(calling.count('"--entrypoint"'), 3)
        self.assertIn('"--entrypoint", "bcftools"', calling)
        self.assertIn('"--entrypoint", "bash"', calling)
        self.assertIn('image, "index"', calling)
        self.assertIn('image, "concat", "-Oz"', calling)

    def test_powershell_tool_preflights_do_not_use_shell_command_substitution(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        alignment_block = script.split("function Assert-ReferenceAlignmentCompatibility", 1)[1].split("function Build-ReferenceTargetPanel", 1)[0]
        bcftools_block = script.split("function Assert-BcftoolsCompatibility", 1)[1].split("function Ensure-ReferenceContigTargets", 1)[0]
        self.assertNotIn('"--entrypoint", "sh"', alignment_block)
        self.assertNotIn('"--entrypoint", "sh"', bcftools_block)
        self.assertNotIn('$(bwa', alignment_block)
        self.assertNotIn('$(samtools', alignment_block)
        self.assertNotIn('$(bcftools', alignment_block + bcftools_block)
        self.assertIn('Invoke-NativeCapture "docker"', alignment_block)
        self.assertIn('Invoke-NativeCapture "docker"', bcftools_block)

    def test_version_locked_preflights_avoid_fragile_help_parsing(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn('$samVersionProbe = Invoke-NativeCapture', script)
        self.assertIn('$versionProbe = Invoke-NativeCapture', script)
        self.assertIn('"--entrypoint", "bash"', script)
        self.assertNotIn('$samHelpProbe', script)
        self.assertNotIn('$mpileupProbe', script)
        self.assertNotIn('$callProbe', script)
        self.assertNotIn('$normProbe', script)
        self.assertNotIn('bcftools norm -h 2>&1 || true', script)
        self.assertNotIn('samtools 2>&1 || true', script)

    def test_bwa_progress_parser_matches_real_mem_messages(self) -> None:
        from reference_alignment import AlignmentProgress
        match = AlignmentProgress.PROCESS_RE.search("[M::mem_process_seqs] Processed 196608 reads in 12.3 CPU sec, 3.4 real sec")
        self.assertIsNotNone(match)
        self.assertEqual(match.group(1), "196608")


    def test_max_performance_profile_replaces_old_four_thread_cap(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn('[string]$PerformanceProfile = "auto"', script)
        self.assertIn('Set-DockerResourceProfile', script)
        self.assertIn('[Math]::Floor($pipelineCpu * $ratio)', script)
        self.assertNotIn('[Math]::Min(4, [Environment]::ProcessorCount - 1)', script)
        launcher = (ROOT / "start.cmd").read_text(encoding="utf-8")
        self.assertIn('-PerformanceProfile max', launcher)

    def test_alignment_has_separate_bwa_sort_and_index_resources(self) -> None:
        source = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        self.assertIn('parser.add_argument("--bwa-threads", type=int)', source)
        self.assertIn('parser.add_argument("--sort-threads", type=int)', source)
        self.assertIn('parser.add_argument("--sort-memory", default="512M")', source)
        self.assertIn('parser.add_argument("--index-threads", type=int)', source)
        self.assertIn('-@ {sort_threads} -m {sort_memory}', source)
        self.assertIn('str(index_threads)', source)

    def test_resource_settings_are_not_added_to_alignment_analytical_identity(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        state_block = script.split('schema = "reference-alignment-v3"', 1)[1].split('Set-AtomicText "$bam.state.json"', 1)[0]
        self.assertNotIn('BwaThreads', state_block)
        self.assertNotIn('SamtoolsSortThreads', state_block)
        self.assertNotIn('SamtoolsSortMemory', state_block)

    def test_stage4_telomere_is_local_and_r1_capable(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        tel = (SRC / "telomere_estimate.py").read_text(encoding="utf-8")
        self.assertIn('"telomere" { Run-TelomereEstimate }', script)
        self.assertIn('Write-Host "STAGE 4 - TELOMERE CALCULATION FROM WGS FASTQ"', script)
        self.assertIn('parser.add_argument("--r2", type=Path)', tel)
        self.assertIn('DEFAULT_K = 7', tel)
        self.assertIn('GENOME_LENGTH_AT_TEL_GC = 332_720_800', tel)
        self.assertIn('TELOMERE_ENDS = 46', tel)

    def test_telomere_report_is_txt_and_metrics_are_cached_separately(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn('$TelomereReportName = "Dante_telomere_report.txt"', script)
        self.assertIn('$TelomereMetricsName = "Dante_telomere_metrics.json"', script)
        self.assertIn('reuse the telomere metrics without rereading FASTQ', script)

    def test_no_backslash_line_continuations(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertFalse(any(line.rstrip().endswith("\\") for line in script.splitlines()))




class ResumableAlignmentV8Tests(unittest.TestCase):
    def test_alignment_uses_resumable_chunk_contract_and_no_full_name_sort(self) -> None:
        script = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        self.assertIn("reference-alignment-checkpoints-v1", script)
        self.assertIn("samtools collate -O -u", script)
        self.assertNotIn("samtools sort -T /work/Dante_work/reference_sort_tmp/name -n", script)
        self.assertIn("samtools merge -u -c -p", script)
        self.assertIn("samtools markdup", script)
        self.assertIn("CHECKPOINT SAVED", script)
        self.assertIn("SAFE STOP", script)

    def test_single_end_path_skips_collate_and_paired_path_uses_collate(self) -> None:
        script = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        self.assertIn("if r2_chunk is None:", script)
        self.assertIn("Single-end duplicate matching only needs coordinate", script)
        self.assertIn("Paired-end fixmate requires mates to be grouped", script)

    def test_max_profile_targets_full_bwa_threads_on_twelve_cpu_host(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn('$PerformanceProfile -eq "max"', script)
        self.assertIn('[Math]::Min(12, $usableCpu)', script)
        self.assertIn('DANTE_REFERENCE_CHECKPOINT_DIR', script)
        self.assertIn('"--chunk-records", "10000000"', script)

    def test_stage15_measured_home_computer_baseline_is_locked(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        alignment = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        self.assertIn('$sort = 4', script)
        self.assertIn('$sortMemory = "1536M"', script)
        self.assertIn('"--chunk-records", "10000000"', script)
        self.assertIn('DOUBLE BUFFER | stage first', alignment)

    def test_checkpoint_manifest_identity_is_strict(self) -> None:
        identity = {
            "schema": reference_alignment.CHECKPOINT_SCHEMA,
            "read_mode": "single",
            "r1_sha256": "a" * 64,
            "r2_sha256": "",
            "r1_records": 10,
            "r2_records": 0,
            "fasta_sha256": "b" * 64,
            "bwa_contract": "bwa",
            "samtools_contract": "samtools",
            "chunk_records": 5,
            "checkpoint_format": "coordinate-bam-v1",
        }
        manifest = {"schema": reference_alignment.CHECKPOINT_SCHEMA, "identity": dict(identity)}
        self.assertTrue(reference_alignment.manifest_matches(manifest, identity))
        changed = dict(identity)
        changed["r1_sha256"] = "c" * 64
        self.assertFalse(reference_alignment.manifest_matches(manifest, changed))

    def test_valid_completed_chunks_requires_contiguous_atomic_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            a = root / "chunk_00000.bam"
            b = root / "chunk_00001.bam"
            a.write_bytes(b"abc")
            b.write_bytes(b"defg")
            manifest = {
                "completed_chunks": [
                    {"index": 0, "bam_name": a.name, "bam_length": 3, "records_per_mate": 5},
                    {"index": 1, "bam_name": b.name, "bam_length": 4, "records_per_mate": 5},
                ]
            }
            rows = reference_alignment.valid_completed_chunks(manifest, root)
            self.assertEqual(len(rows), 2)
            b.write_bytes(b"broken")
            rows = reference_alignment.valid_completed_chunks(manifest, root)
            self.assertEqual(len(rows), 1)



class ResumableCallingV82Tests(unittest.TestCase):
    def test_target_bed_chunking_is_contiguous_and_resumable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / "targets.bed"
            bed.write_text("".join(f"1\t{i}\t{i+1}\n" for i in range(7)), encoding="utf-8")
            self.assertEqual(reference_calling.count_target_sites(bed), 7)
            with bed.open("rt", encoding="utf-8") as handle:
                first = root / "first.bed"
                second = root / "second.bed"
                self.assertEqual(reference_calling.write_target_chunk(handle, first, 3), 3)
                self.assertEqual(reference_calling.write_target_chunk(handle, second, 2), 2)
                self.assertEqual(reference_calling.skip_target_sites(handle, 2), 2)
            self.assertEqual(len(first.read_text().splitlines()), 3)
            self.assertEqual(len(second.read_text().splitlines()), 2)

    def test_valid_call_checkpoints_require_contiguous_target_boundaries(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = []
            start = 0
            for index, count in enumerate((3, 2)):
                vcf = root / f"call_chunk_{index:05d}.vcf.gz"
                tbi = Path(str(vcf) + ".tbi")
                vcf.write_bytes(b"vcf" + bytes([index]))
                tbi.write_bytes(b"tbi" + bytes([index]))
                rows.append({
                    "index": index,
                    "target_start": start,
                    "target_count": count,
                    "vcf_name": vcf.name,
                    "vcf_length": vcf.stat().st_size,
                    "tbi_length": tbi.stat().st_size,
                })
                start += count
            manifest = {"completed_chunks": rows}
            self.assertEqual(len(reference_calling.valid_completed_chunks(manifest, root)), 2)
            rows[1]["target_start"] = 99
            self.assertEqual(len(reference_calling.valid_completed_chunks(manifest, root)), 1)

    def test_call_checkpoint_identity_excludes_chunk_size(self) -> None:
        source = (SRC / "reference_calling.py").read_text(encoding="utf-8")
        identity_block = source.split("def expected_identity", 1)[1].split("def manifest_matches", 1)[0]
        self.assertNotIn("chunk_sites", identity_block)
        self.assertIn('"checkpoint_format": "bgzip-vcf-tbi-known-alleles-coordinate-safe-v4"', identity_block)

    def test_stage16_contract_has_real_progress_resume_disk_guard_and_concat(self) -> None:
        source = (SRC / "reference_calling.py").read_text(encoding="utf-8")
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn("reference-call-checkpoints-v4-known-alleles-ploidy-coordinate-safe", source)
        self.assertIn("CHECKPOINT SAVED", source)
        self.assertIn("target loci preserved", source)
        self.assertIn("SAFE STOP", source)
        self.assertIn('"concat", "-Oz"', source)
        self.assertIn("TransientDockerError", source)
        self.assertIn("wait_for_docker_engine(120)", source)
        self.assertIn("DANTE_REFERENCE_CALL_CHECKPOINT_DIR", script)
        self.assertIn("DANTE_REFERENCE_CALL_CHUNK_SITES", script)
        self.assertIn('"--chunk-sites", "$callChunkSites"', script)
        self.assertIn('-WatchDirectory $callCheckpointDirectory', script)
        self.assertIn('calling_strategy = "resumable-v4-checkpoints-known-alleles-ploidy-plus-PL-reference-quality-v5"', script)

    def test_legacy_reference_calls_v1_is_intentionally_stale(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertNotIn('$state.schema -eq "reference-calls-v1"', script)
        self.assertIn('$state.schema -eq "reference-calls-v5-known-alleles-ploidy-refq"', script)
        self.assertIn('schema = "reference-calls-v5-known-alleles-ploidy-refq"', script)



    def test_stage16_uses_known_alleles_gq_and_ploidy_contract(self) -> None:
        source = (SRC / "reference_calling.py").read_text(encoding="utf-8")
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn("-C alleles -T", source)
        self.assertIn("-f GQ", source)
        self.assertIn("--ploidy-file", source)
        self.assertIn("sample_sex.tsv", source)
        self.assertIn("hom_ref_reference_quality_fraction", source)
        self.assertIn("DANTE_SAMPLE_SEX", script)
        self.assertIn("STALE: PRE-REFQ CALL MODEL", script)
        self.assertIn("STALE: PRE-COORDINATE-SAFE CHECKPOINTS", script)
        self.assertTrue((ROOT / "reference_data" / "GRCh37_bcftools_ploidy.tsv").is_file())

    def test_target_locus_stream_unions_known_alts_and_falls_back_to_all_snps(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "ref.fa"
            fasta.write_text(">1\nACGTACGTACGT\n", encoding="ascii")
            (root / "ref.fa.fai").write_text("1\t12\t3\t12\t13\n", encoding="ascii")
            targets = root / "targets.tsv"
            targets.write_text(
                "rsid\tchromosome\tposition\tref\talts\tsources\treference_genotype\treference_genotype_status\treference_genotype_sources\n"
                "rs1\t1\t2\tC\tT\tx\t--\tMISSING\t.\n"
                "rs2\t1\t2\tC\tA\tx\t--\tMISSING\t.\n"
                "rs3\t1\t3\t.\t.\tx\t--\tMISSING\t.\n",
                encoding="utf-8",
            )
            with reference_calling.FastaIndex(fasta) as ref:
                loci = list(reference_calling.iter_target_loci(targets, ref))
            self.assertEqual(len(loci), 2)
            self.assertEqual(loci[0].ref, "C")
            self.assertEqual(set(loci[0].alts), {"A", "T"})
            self.assertFalse(loci[0].fallback_all_snps)
            self.assertEqual(loci[1].ref, "G")
            self.assertEqual(set(loci[1].alts), {"A", "C", "T"})
            self.assertTrue(loci[1].fallback_all_snps)

    def test_reference_calling_print_contract_identifies_current_release(self) -> None:
        self.assertEqual(reference_calling.PIPELINE_RELEASE, "9.0.0")
        self.assertEqual(reference_calling.CALLING_SCHEMA, "reference-calls-v5-known-alleles-ploidy-refq")
        self.assertEqual(reference_calling.CHECKPOINT_SCHEMA, "reference-call-checkpoints-v4-known-alleles-ploidy-coordinate-safe")

    def test_non_acgt_grch37_targets_are_skipped_and_reported_not_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "ref.fa"
            fasta.write_text(">1\nACNTACGT\n", encoding="ascii")
            (root / "ref.fa.fai").write_text("1\t8\t3\t8\t9\n", encoding="ascii")
            targets = root / "targets.tsv"
            targets.write_text(
                "rsid\tchromosome\tposition\tref\talts\tsources\treference_genotype\treference_genotype_status\treference_genotype_sources\n"
                "rs1\t1\t2\tC\tT\tx\t--\tMISSING\t.\n"
                "rsN\t1\t3\t.\t.\tx\t--\tMISSING\t.\n"
                "rs2\t1\t4\tT\tC\tx\t--\tMISSING\t.\n",
                encoding="utf-8",
            )
            report = root / "uncallable.tsv"
            with reference_calling.FastaIndex(fasta) as ref:
                loci = list(reference_calling.iter_target_loci(targets, ref))
                self.assertEqual(len(loci), 3)
                self.assertTrue(loci[0].callable)
                self.assertFalse(loci[1].callable)
                self.assertEqual(loci[1].ref, "N")
                self.assertEqual(loci[1].uncallable_reason, "NON_ACGT_REFERENCE:N")
                self.assertTrue(loci[2].callable)
                stats = reference_calling.scan_target_reference_contract(targets, ref, report)

            self.assertEqual(stats["panel_loci"], 3)
            self.assertEqual(stats["callable_loci"], 2)
            self.assertEqual(stats["skipped_non_acgt_reference_loci"], 1)
            self.assertEqual(stats["uncallable_reference_symbols"], {"N": 1})
            report_text = report.read_text(encoding="utf-8")
            self.assertIn("1\t3\tN\tNON_ACGT_REFERENCE:N", report_text)

    def test_sex_auto_inference_has_ambiguous_guard(self) -> None:
        targets = {"1": 1000, "X": 100, "Y": 50}
        male = reference_calling.infer_sex_from_normalized_depth(
            {"1": 10000, "X": 500, "Y": 250}, targets
        )
        female = reference_calling.infer_sex_from_normalized_depth(
            {"1": 10000, "X": 1000, "Y": 5}, targets
        )
        ambiguous = reference_calling.infer_sex_from_normalized_depth(
            {"1": 10000, "X": 800, "Y": 250 // 4}, targets
        )
        self.assertEqual(male["inferred"], "M")
        self.assertEqual(female["inferred"], "F")
        self.assertEqual(ambiguous["inferred"], "AMBIGUOUS")

    def test_call_contract_requires_hom_ref_gq_and_enforces_ploidy(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            good = root / "good.vcf.gz"
            with gzip.open(good, "wt", encoding="utf-8") as out:
                out.write("##fileformat=VCFv4.2\n")
                out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
                for pos in range(1, 101):
                    out.write(f"1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT:GQ:DP:AD\t0/0:60:20:20,0\n")
                out.write("X\t3000000\t.\tA\tC\t.\tPASS\t.\tGT:GQ:DP:AD\t0:60:20:20,0\n")
                out.write("Y\t1000000\t.\tA\tC\t.\tPASS\t.\tGT:GQ:DP:AD\t0:60:20:20,0\n")
                out.write("MT\t100\t.\tA\tC\t.\tPASS\t.\tGT:GQ:DP:AD\t0:60:20:20,0\n")
            qc = reference_calling.validate_call_contract(good, "M")
            self.assertEqual(qc["hom_ref_with_native_gq"], qc["hom_ref_records"])
            self.assertEqual(qc["ploidy_mismatch_records"], 0)

            bad = root / "bad.vcf.gz"
            with gzip.open(bad, "wt", encoding="utf-8") as out:
                out.write("##fileformat=VCFv4.2\n")
                out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
                for pos in range(1, 101):
                    out.write(f"1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT:DP:AD\t0/0:20:20,0\n")
            with self.assertRaisesRegex(RuntimeError, "native GQ or a PL-derived reference quality"):
                reference_calling.validate_call_contract(bad, "M")

    def test_reference_export_rejects_alt_not_known_for_rsid(self) -> None:
        from reference_genotype_export import evaluate_call as evaluate
        call = RawCall(
            "1", 100, "A", ("C", "G"),
            {"GT": "0/2", "DP": "20", "GQ": "60", "AD": "10,0,10"},
            {},
        )
        result = evaluate(call, "A", "C")
        self.assertEqual(result.tier, "NONE")
        self.assertEqual(result.reason, "ALT_NOT_IN_TARGET")




class OpenReferenceBundleTests(unittest.TestCase):
    def test_open_reference_catalog_is_license_allowlisted_and_complete(self) -> None:
        import open_reference_manager as orm
        self.assertEqual(orm.self_test(), 0)
        sources = {resource.source for resource in orm.RESOURCES}
        self.assertEqual(
            {"HGNC", "dbSNP", "GIAB", "GenCC", "ClinGen", "Mondo", "HPO"},
            sources,
        )
        self.assertTrue(all(r.license_class in orm.LICENSE_CLASSES for r in orm.RESOURCES))

    def test_hpo_is_ontology_only_and_full_dbsnp_vcfs_are_not_bundled(self) -> None:
        import open_reference_manager as orm
        hpo = [r for r in orm.RESOURCES if r.source == "HPO"]
        self.assertEqual([r.key for r in hpo], ["ontology_obo"])
        full_vcf_names = [
            r.relative_path for r in orm.RESOURCES
            if r.source == "dbSNP" and r.relative_path.endswith(".gz")
        ]
        self.assertEqual(full_vcf_names, [])
        self.assertIn("GRCh37.p13", orm.DBSNP_REMOTE_VCFS)
        self.assertIn("GRCh38.p14", orm.DBSNP_REMOTE_VCFS)

    def test_alphamissense_is_research_only_not_open(self) -> None:
        import open_reference_manager as orm
        import research_reference_manager as rrm
        self.assertFalse(any(r.source == "AlphaMissense" for r in orm.RESOURCES))
        alpha = next(r for r in rrm.RESOURCES if r.source == "AlphaMissense")
        self.assertEqual(alpha.expected_md5, rrm.ALPHAMISSENSE_MD5)
        self.assertIn("NC-SA", alpha.license_name)
        self.assertIn("commercial", alpha.prohibited_uses)

    def test_research_bundle_reuses_verified_migrated_artifact_without_network(self) -> None:
        import hashlib
        import research_reference_manager as rrm
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rel = "AlphaMissense/tiny.tsv.gz"
            dest = root / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(dest, "wb") as out:
                out.write(b"# licence\nchr1\t1\tA\tC\n")
            md5 = hashlib.md5(dest.read_bytes()).hexdigest()
            resource = rrm.Resource(
                "AlphaMissense", "tiny", "https://example.invalid/never-fetch", rel,
                "CC-BY-NC-SA-4.0", rrm.ALPHAMISSENSE_LICENSE_URL, "test",
                content_kind="gzip", expected_md5=md5, prohibited_uses=("commercial",),
            )
            with mock.patch.object(rrm.urllib.request, "urlopen", side_effect=AssertionError("network must not be used")):
                rec = rrm.download(resource, root, previous={}, force=False)
            self.assertEqual(rec["transport"], "reused-verified-local")
            self.assertEqual(rec["md5"], md5)
            self.assertEqual(rec["length"], dest.stat().st_size)


class GenomewideVUSV8Tests(unittest.TestCase):
    def test_chess_reference_is_embedded_and_whole_genome_rule_is_present(self) -> None:
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        shortlist = (SRC / "build_shortlist.py").read_text(encoding="utf-8")
        chess = ROOT / "reference_data" / "CHESS_2.2_genes.tsv"
        self.assertTrue(chess.is_file())
        self.assertIn("deterministic-shortlist-v3-chess-genomewide", shortlist)
        self.assertIn("classify_genomewide_uncertain", shortlist)
        self.assertIn("Dante_genomewide_uncertain_coding.tsv", script)
        self.assertIn("genome-wide VUS/rare coding sweep", script)
        self.assertNotIn("priority_genes.txt", script)

    def test_genomewide_rule_does_not_require_chess_membership_or_damaging_predictor(self) -> None:
        import build_shortlist
        row = {
            "alt_carried": "yes", "vcf_filter": "PASS", "symbol": "NOT_IN_CHESS",
            "all_symbols": "NOT_IN_CHESS", "clinvar_current_significance": "Uncertain significance",
            "clin_sig": ".", "clinvar_current_review_status": ".", "max_af": "0.2",
            "impact": "MODERATE", "all_consequences": "missense_variant", "sift": "tolerated", "polyphen": "benign",
        }
        result = build_shortlist.classify_genomewide_uncertain(row, 0.01)
        self.assertIsNotNone(result)
        self.assertEqual(result[1], "GENOME_VUS_MISSENSE")

    def test_unclassified_common_missense_is_not_promoted_but_rare_is(self) -> None:
        import build_shortlist
        row = {
            "alt_carried": "yes", "vcf_filter": "PASS", "symbol": "GENE1", "all_symbols": "GENE1",
            "clinvar_current_significance": ".", "clin_sig": ".", "clinvar_current_review_status": ".",
            "impact": "MODERATE", "all_consequences": "missense_variant", "sift": "tolerated", "polyphen": "benign",
        }
        common = dict(row, max_af="0.05")
        rare = dict(row, max_af="0.005")
        self.assertIsNone(build_shortlist.classify_genomewide_uncertain(common, 0.01))
        result = build_shortlist.classify_genomewide_uncertain(rare, 0.01)
        self.assertIsNotNone(result)
        self.assertEqual(result[1], "GENOME_UNCLASSIFIED_MISSENSE")

    def test_stage15_v8_has_orphan_container_labels_and_cleanup(self) -> None:
        alignment = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        script = (SRC / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn('biodata.dante.stage=15', alignment)
        self.assertIn('cleanup_stale_stage15_containers()', alignment)
        self.assertIn('function Remove-Stage15OrphanContainers', script)
        self.assertIn('Remove-Stage15OrphanContainers', script)

    def test_activity_watch_status_is_not_repeated_unchanged(self) -> None:
        activity = (SRC / "dante_activity.py").read_text(encoding="utf-8")
        self.assertIn("last_watch_status", activity)
        self.assertIn("if watched == last_watch_status", activity)


    def test_v801_transient_docker_recovery_contract(self):
        script = (SRC / "reference_alignment.py").read_text(encoding="utf-8")
        self.assertIn("class TransientDockerError", script)
        self.assertIn("wait_for_docker_engine(120)", script)
        self.assertIn("max_chunk_attempts = 3", script)
        self.assertIn("DOUBLE BUFFER | stage first", script)
        self.assertIn("staged FASTQ + completed checkpoints preserved", script)
        self.assertIn("bwa_code == 125", script)
        self.assertIn("completed checkpoints are preserved", script)

if __name__ == "__main__":
    unittest.main()
