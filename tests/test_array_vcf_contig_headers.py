from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import array_vcf_tools as avt  # noqa: E402


class ArrayVcfContigHeaderRegressionTests(unittest.TestCase):
    def test_genotek_vcf_without_contigs_gets_canonical_grch37_header(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "genotek.vcf"
            source.write_text(
                "##fileformat=VCFv4.2\n"
                "##source=Genotek\n"
                "##reference=hg19\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                "chr1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0|1\n"
                "chr2\t200\trs2\tC\tT\t.\t.\t.\tGT\t1|1\n"
                "chrM\t300\trs3\tA\tG\t.\t.\t.\tGT\t0|1\n",
                encoding="utf-8",
            )
            snp = root / "snp.vcf.gz"
            indel = root / "indel.vcf.gz"
            residual = root / "residual.vcf.gz"
            avt.prepare_layers(source, snp, indel, residual, root / "report.txt")

            with gzip.open(snp, "rt", encoding="utf-8") as handle:
                lines = handle.readlines()

        contigs = {
            line.split("ID=", 1)[1].split(",", 1)[0]
            for line in lines
            if line.startswith("##contig=<ID=")
        }
        data_contigs = {
            line.split("\t", 1)[0]
            for line in lines
            if line and not line.startswith("#")
        }
        self.assertTrue({"1", "2", "MT"}.issubset(contigs))
        self.assertTrue(data_contigs.issubset(contigs))
        self.assertFalse(any(x.startswith("chr") for x in contigs))
        self.assertIn("##contig=<ID=1,length=249250621>\n", lines)
        self.assertIn("##contig=<ID=MT,length=16569>\n", lines)

    def test_source_chr_contig_declarations_are_replaced_not_duplicated(self) -> None:
        source_header = [
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=249250621>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n",
        ]
        header = avt.normalized_header(source_header, source_header[-1])
        text = "".join(header)
        self.assertNotIn("ID=chr1", text)
        self.assertEqual(text.count("##contig=<ID=1,length=249250621>"), 1)


if __name__ == "__main__":
    unittest.main()
