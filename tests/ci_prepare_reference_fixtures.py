from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    refs = ROOT / "reference_data"
    refs.mkdir(exist_ok=True)
    (refs / "CHESS_2.2_genes.tsv").write_text(
        "Gene_Name\tDescription\tGene_Type\tGFF_ID\tLocation\tDatabase\tRefSeq_GeneID\n"
        "CI_TEST\tCI-only synthetic fixture\tprotein_coding\tCI_TEST\t1:1-10\tCI\t-\n",
        encoding="utf-8",
        newline="\n",
    )
    (refs / "GRCh37_bcftools_ploidy.tsv").write_text(
        "# CI-only presence fixture; production file is local reference data.\n",
        encoding="utf-8",
        newline="\n",
    )
    print("CI reference presence fixtures prepared without modifying tracked source files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
