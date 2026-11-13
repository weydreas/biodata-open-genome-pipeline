from __future__ import annotations

import csv
import collections
import gzip
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

import knowledge_enrichment as ke  # noqa: E402
import variant_research_agent as vra  # noqa: E402


def dump_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def artifact(path: str) -> dict:
    return {"path": path}


def make_bundles(root: Path) -> tuple[Path, Path]:
    open_root = root / "open"
    research_root = root / "research"
    open_root.mkdir()
    research_root.mkdir()

    hgnc = open_root / "hgnc.tsv"
    hgnc.write_text(
        "hgnc_id\tsymbol\tname\talias_symbol\tprev_symbol\n"
        "HGNC:1\tTESTAR\tAR test gene\tOLDAR\t.\n"
        "HGNC:2\tTESTAD\tAD test gene\t.\t.\n",
        encoding="utf-8",
    )
    gencc = open_root / "gencc.tsv"
    gencc.write_text(
        "gene_symbol\tdisease_curie\tdisease_title\tclassification_title\tmoi_title\n"
        "TESTAR\tMONDO:0000001\tAR disorder\tDefinitive\tAutosomal recessive\n"
        "TESTAD\tMONDO:0000002\tAD disorder\tStrong\tAutosomal dominant\n",
        encoding="utf-8",
    )
    clingen = open_root / "clingen.csv"
    clingen.write_text(
        "Gene Symbol,Disease ID,Disease Label,Classification,Mode of Inheritance\n"
        "TESTAD,MONDO:0000002,AD disorder,Definitive,Autosomal dominant\n",
        encoding="utf-8",
    )
    dosage = open_root / "dosage.csv"
    dosage.write_text(
        "Gene Symbol,Haploinsufficiency Score,Triplosensitivity Score\n"
        "TESTAD,3,0\n",
        encoding="utf-8",
    )
    actionability = open_root / "actionability.json"
    dump_json(actionability, {"records": [{"condition": "AD disorder", "gene": {"symbol": "TESTAD"}}]})
    mondo = open_root / "mondo.obo"
    mondo.write_text(
        "[Term]\nid: MONDO:0000001\nname: AR disorder\n\n"
        "[Term]\nid: MONDO:0000002\nname: AD disorder\n",
        encoding="utf-8",
    )
    hpo = open_root / "hp.obo"
    hpo.write_text("[Term]\nid: HP:0000001\nname: All\n", encoding="utf-8")

    panel_db = research_root / "panelapp.sqlite3"
    con = sqlite3.connect(panel_db)
    con.execute("CREATE TABLE panels(panel_id INTEGER PRIMARY KEY, name TEXT, version TEXT, fetched_utc TEXT, payload_json TEXT)")
    con.execute(
        "INSERT INTO panels VALUES(?,?,?,?,?)",
        (
            1,
            "Synthetic panel",
            "1.0",
            "2026-01-01T00:00:00Z",
            json.dumps({
                "genes": [
                    {
                        "entity_name": "TESTAD",
                        "confidence_level": "3",
                        "mode_of_inheritance": "MONOALLELIC",
                        "phenotypes": ["AD disorder"],
                    }
                ]
            }),
        ),
    )
    con.commit()
    con.close()

    open_manifest = {
        "schema": "open-reference-bundle-v2",
        "status": "complete",
        "installed_utc": "2026-01-01T00:00:00Z",
        "sources": {
            "HGNC": {"artifacts": {"complete": artifact("hgnc.tsv")}},
            "GenCC": {"artifacts": {"submissions_new_tsv": artifact("gencc.tsv")}},
            "ClinGen": {"artifacts": {
                "gene_disease_validity": artifact("clingen.csv"),
                "dosage_genes_regions": artifact("dosage.csv"),
                "actionability_adult": artifact("actionability.json"),
            }},
            "Mondo": {"artifacts": {"ontology_obo": artifact("mondo.obo")}},
            "HPO": {"artifacts": {"ontology_obo": artifact("hp.obo")}},
        },
    }
    research_manifest = {
        "schema": "research-reference-bundle-v2",
        "status": "complete",
        "installed_utc": "2026-01-01T00:00:00Z",
        "sources": {
            "PanelApp": {"artifacts": {"snapshot_sqlite": artifact("panelapp.sqlite3")}},
            "AlphaMissense": {"artifacts": {}},
        },
    }
    dump_json(open_root / "open_reference_manifest.json", open_manifest)
    dump_json(research_root / "research_reference_manifest.json", research_manifest)
    return open_root, research_root


def write_catalog(path: Path) -> None:
    fields = [
        "variant_key", "chrom", "position", "ref", "alt", "alt_carried", "vcf_filter",
        "source_layer", "GT", "symbol", "all_symbols", "all_consequences", "impact",
        "max_af", "rsid", "clinvar_exact_match", "clinvar_current_rsid",
        "clinvar_current_significance", "clinvar_current_review_status", "clinvar_review_stars",
        "clinvar_current_conditions", "clinvar_current_geneinfo", "hgvsc", "hgvsp", "sift", "polyphen", "DP", "GQ", "AD",
    ]
    rows = [
        # Two rare heterozygous LoF variants in a definitive AR gene.
        ["GRCh37:1:100:A:T", "1", "100", "A", "T", "yes", "PASS", "DANTE", "0/1", "OLDAR", "OLDAR", "stop_gained", "HIGH", "0.0001", "rs100", "no", ".", ".", ".", "0", ".", ".", ".", ".", ".", ".", "30", "99", "15,15"],
        ["GRCh37:1:200:G:C", "1", "200", "G", "C", "yes", "PASS", "DANTE", "0/1", "TESTAR", "TESTAR", "frameshift_variant", "HIGH", "0.0002", "rs200", "no", ".", ".", ".", "0", ".", ".", ".", ".", ".", ".", "30", "99", "16,14"],
        # Rare missense in actionable / PanelApp AD gene.
        ["GRCh37:2:300:C:T", "2", "300", "C", "T", "yes", "PASS", "DANTE", "0/1", "TESTAD", "TESTAD", "missense_variant", "MODERATE", "0.0001", "rs300", "no", ".", ".", ".", "0", ".", ".", ".", ".", ".", ".", "30", "99", "15,15"],
        # Common missense should not be materialized merely because the gene is interesting.
        ["GRCh37:2:400:A:G", "2", "400", "A", "G", "yes", "PASS", "DANTE", "0/1", "TESTAD", "TESTAD", "missense_variant", "MODERATE", "0.25", "rs400", "no", ".", ".", ".", "0", ".", ".", ".", ".", ".", ".", "30", "99", "15,15"],
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


class KnowledgeEnrichmentTests(unittest.TestCase):
    def test_reference_index_canonicalizes_and_compiles_open_and_research_layers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            knowledge = ke.GeneKnowledge(db)
            self.assertEqual(knowledge.canonicalize("OLDAR")[0], "TESTAR")
            self.assertTrue(knowledge.profile("TESTAR").strong_relations)
            profile = knowledge.profile("TESTAD")
            self.assertTrue(profile.actionable)
            self.assertTrue(profile.dosage)
            self.assertEqual(profile.panelapp[0]["panel_name"], "Synthetic panel")

    def test_full_catalog_scan_filters_and_marks_unphased_ar_pair(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            knowledge = ke.GeneKnowledge(db)
            catalog = root / "catalog.tsv.gz"
            write_catalog(catalog)
            events, counts = ke.scan_small_catalog(catalog, knowledge, 0.01, progress_every=0)
            self.assertEqual(counts["catalog_rows"], 4)
            self.assertEqual(len(events), 3)
            self.assertNotIn("GRCh37:2:400:A:G", events)
            self.assertEqual(ke.mark_possible_biallelic(events, knowledge, 0.01), 2)
            self.assertTrue(events["GRCh37:1:100:A:T"].possible_biallelic_unphased)
            missense = events["GRCh37:2:300:C:T"]
            self.assertIn("K3_RARE_MISSENSE_VALIDATED_GENE", missense.open_reasons)
            self.assertIn("K4_CLINGEN_ACTIONABLE", missense.open_reasons)
            self.assertIn("K7_PANELAPP_RESEARCH", missense.research_reasons)
            self.assertIn("K8_CADD_SUPPORT_TARGET", missense.research_reasons)

    def test_open_and_research_rows_are_physically_separated(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            knowledge = ke.GeneKnowledge(db)
            catalog = root / "catalog.tsv.gz"
            write_catalog(catalog)
            events, _ = ke.scan_small_catalog(catalog, knowledge, 0.01, progress_every=0)
            missense = events["GRCh37:2:300:C:T"]
            missense.alpha_score = "0.98"
            missense.alpha_class = "likely_pathogenic"
            missense.research_reasons.add("K8_ALPHAMISSENSE_SUPPORT")
            missense.cadd_raw = "4.2"
            missense.cadd_phred = "27.1"
            open_rows = ke.open_findings(events, knowledge)
            research_rows = ke.research_overlays(events, open_rows, knowledge)
            self.assertFalse(any("AlphaMissense" in row["evidence_sources"] for row in open_rows if row["event_key"] == missense.event_key))
            self.assertNotIn("alpha_missense_pathogenicity", open_rows[0])
            self.assertNotIn("cadd_phred", open_rows[0])
            self.assertNotIn("panelapp_panels", open_rows[0])
            k7 = next(row for row in research_rows if row["reason_code"] == "K7_PANELAPP_RESEARCH")
            ka = next(row for row in research_rows if row["reason_code"] == "K8_ALPHAMISSENSE_SUPPORT")
            kc = next(row for row in research_rows if row["reason_code"] == "K8_CADD_SUPPORT_TARGET")
            self.assertEqual(k7["cadd_phred"], ".")
            self.assertNotEqual(k7["panelapp_panels"], ".")
            self.assertEqual(ka["alpha_missense_pathogenicity"], "0.98")
            self.assertEqual(ka["cadd_phred"], ".")
            self.assertEqual(kc["panelapp_panels"], ".")
            self.assertEqual(kc["cadd_phred"], "27.1")

    def test_duplicate_gencc_clingen_assertions_collapse_to_one_semantic_finding(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            con = sqlite3.connect(db)
            try:
                con.execute(
                    "INSERT OR IGNORE INTO gene_disease(gene,hgnc_id,disease_id,disease_name,classification,inheritance,source) VALUES(?,?,?,?,?,?,?)",
                    ("TESTAD", "HGNC:2", "MONDO:0000002", "AD disorder", "Strong", "Autosomal dominant", "ClinGen"),
                )
                con.commit()
            finally:
                con.close()
            knowledge = ke.GeneKnowledge(db)
            catalog = root / "catalog.tsv.gz"
            write_catalog(catalog)
            events, _ = ke.scan_small_catalog(catalog, knowledge, 0.01, progress_every=0)
            open_rows = ke.open_findings(events, knowledge)
            k3 = [
                row for row in open_rows
                if row["event_key"] == "GRCh37:2:300:C:T" and row["reason_code"] == "K3_RARE_MISSENSE_VALIDATED_GENE"
            ]
            self.assertEqual(len(k3), 1)
            self.assertEqual(k3[0]["gene_disease_classification"], "DEFINITIVE")
            self.assertIn("GenCC", k3[0]["evidence_sources"])
            self.assertIn("ClinGen", k3[0]["evidence_sources"])

    def test_diseases_without_ids_do_not_collide_in_finding_identity(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            con = sqlite3.connect(db)
            try:
                con.execute("DELETE FROM gene_disease WHERE gene='TESTAD'")
                con.executemany(
                    "INSERT INTO gene_disease(gene,hgnc_id,disease_id,disease_name,classification,inheritance,source) VALUES(?,?,?,?,?,?,?)",
                    [
                        ("TESTAD", "HGNC:2", "", "Disease Alpha", "Strong", "Autosomal dominant", "GenCC"),
                        ("TESTAD", "HGNC:2", "", "Disease Beta", "Strong", "Autosomal dominant", "ClinGen"),
                    ],
                )
                con.commit()
            finally:
                con.close()
            knowledge = ke.GeneKnowledge(db)
            catalog = root / "catalog.tsv.gz"
            write_catalog(catalog)
            events, _ = ke.scan_small_catalog(catalog, knowledge, 0.01, progress_every=0)
            rows = [
                row for row in ke.open_findings(events, knowledge)
                if row["event_key"] == "GRCh37:2:300:C:T" and row["reason_code"] == "K3_RARE_MISSENSE_VALIDATED_GENE"
            ]
            self.assertEqual({row["disease_name"] for row in rows}, {"Disease Alpha", "Disease Beta"})
            self.assertEqual(len({row["finding_id"] for row in rows}), 2)


    def test_inheritance_bucket_is_idempotent_for_canonical_tokens(self) -> None:
        self.assertEqual(ke.inheritance_bucket("AR"), "AR")
        self.assertEqual(ke.inheritance_bucket("AD"), "AD")
        self.assertEqual(ke.inheritance_bucket("X_LINKED"), "X_LINKED")
        self.assertEqual(ke.inheritance_bucket("MITOCHONDRIAL"), "MITOCHONDRIAL")
        self.assertEqual(ke.inheritance_bucket("UNKNOWN"), "UNKNOWN")
        self.assertEqual(ke.inheritance_bucket("X-linked recessive"), "X_LINKED")

    def test_unknown_moi_is_absorbed_when_same_disease_has_one_specific_moi(self) -> None:
        profile = ke.GeneProfile("XGENE", "HGNC:999")
        profile.relations = [
            ke.GeneDisease("XGENE", "HGNC:999", "MONDO:9", "X disorder", "Definitive", "X-linked", "GenCC"),
            ke.GeneDisease("XGENE", "HGNC:999", "MONDO:9", "X disorder", "Strong", "", "ClinGen"),
        ]
        event = ke.PatientEvent("SNV_INDEL", "GRCh37:X:1:A:T", {"DANTE"}, "X", 1, 1)
        rows = ke.relation_rows(event, "XGENE", profile)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].inheritance, "X_LINKED")
        self.assertEqual(set(rows[0].source.split(";")), {"ClinGen", "GenCC"})

    def test_open_findings_do_not_duplicate_xlinked_plus_unknown_assertions(self) -> None:
        class MiniKnowledge:
            def __init__(self) -> None:
                self.p = ke.GeneProfile("XGENE", "HGNC:999")
                self.p.relations = [
                    ke.GeneDisease("XGENE", "HGNC:999", "MONDO:9", "X disorder", "Definitive", "X-linked", "GenCC"),
                    ke.GeneDisease("XGENE", "HGNC:999", "MONDO:9", "X disorder", "Definitive", "", "ClinGen"),
                ]
            def profile(self, symbol: str):
                return self.p if symbol == "XGENE" else None

        event = ke.PatientEvent(
            "SNV_INDEL", "GRCh37:X:123:C:G", {"DANTE"}, "X", 123, 123,
            ref="C", alt="G", genotype="1", zygosity="HEMIZYGOUS_ALT", genes={"XGENE"},
            primary_gene="XGENE", consequence="missense_variant", impact="MODERATE", max_af="0.0001",
        )
        event.open_reasons.add("K3_RARE_MISSENSE_VALIDATED_GENE")
        rows = ke.open_findings({event.event_key: event}, MiniKnowledge())
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["inheritance"], "X_LINKED")
        self.assertIn("GenCC", rows[0]["evidence_sources"])
        self.assertIn("ClinGen", rows[0]["evidence_sources"])



    def test_clinvar_bucket_does_not_promote_conflicts_or_association_not_found(self) -> None:
        self.assertEqual(ke.clinvar_bucket("Conflicting_classifications_of_pathogenicity"), "CONFLICT")
        self.assertEqual(ke.clinvar_bucket("Pathogenic/Likely_pathogenic"), "PATHOGENIC")
        self.assertEqual(ke.clinvar_bucket("association_not_found"), "OTHER")
        self.assertEqual(ke.clinvar_bucket("risk_factor"), "ASSOCIATION")

    def test_event_level_gene_evidence_does_not_leak_to_secondary_genes(self) -> None:
        class MiniKnowledge:
            def __init__(self) -> None:
                self.main = ke.GeneProfile("TESTAD", "HGNC:2")
                self.main.relations = [
                    ke.GeneDisease("TESTAD", "HGNC:2", "MONDO:2", "AD disorder", "Definitive", "Autosomal dominant", "GenCC")
                ]
                self.main.actionable = True
                self.other = ke.GeneProfile("OTHER", "HGNC:3")
            def profile(self, symbol: str):
                return {"TESTAD": self.main, "OTHER": self.other}.get(symbol)

        event = ke.PatientEvent(
            "SNV_INDEL", "GRCh37:2:300:C:T", {"DANTE"}, "2", 300, 300,
            ref="C", alt="T", genotype="0/1", zygosity="HET_ALT",
            genes={"TESTAD", "OTHER"}, primary_gene="TESTAD",
            consequence="missense_variant", impact="MODERATE", max_af="0.0001",
        )
        event.open_reasons.update({"K3_RARE_MISSENSE_VALIDATED_GENE", "K4_CLINGEN_ACTIONABLE"})
        rows = ke.open_findings({event.event_key: event}, MiniKnowledge())
        self.assertEqual({row["gene"] for row in rows}, {"TESTAD"})
        self.assertTrue(all(row["gene_disease_classification"] != "." for row in rows if row["reason_code"] == "K3_RARE_MISSENSE_VALIDATED_GENE"))
        self.assertTrue(all(row["actionability"] == "yes" for row in rows if row["reason_code"] == "K4_CLINGEN_ACTIONABLE"))

    def test_exact_clinvar_findings_are_scoped_to_clinvar_geneinfo(self) -> None:
        class MiniKnowledge:
            def __init__(self) -> None:
                self.a = ke.GeneProfile("CLINGENE", "HGNC:10")
                self.b = ke.GeneProfile("VEPNEIGHBOR", "HGNC:11")
            def profile(self, symbol: str):
                return {"CLINGENE": self.a, "VEPNEIGHBOR": self.b}.get(symbol)

        event = ke.PatientEvent(
            "SNV_INDEL", "GRCh37:1:10:A:T", {"DANTE"}, "1", 10, 10,
            ref="A", alt="T", genes={"CLINGENE", "VEPNEIGHBOR"}, primary_gene="VEPNEIGHBOR",
            clinvar_genes={"CLINGENE"}, clinvar_significance="Pathogenic", clinvar_conditions="Synthetic disorder",
        )
        event.open_reasons.add("K1_EXACT_CLINVAR_PLP")
        rows = ke.open_findings({event.event_key: event}, MiniKnowledge())
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["gene"], "CLINGENE")

    def test_dosage_region_requires_direction_and_substantial_coverage(self) -> None:
        region = ke.DosageRegion("ISCA:1", "Synthetic region", "1", 1000, 1999, "3", "0")
        self.assertTrue(ke.dosage_region_supports(region, "DEL"))
        self.assertFalse(ke.dosage_region_supports(region, "DUP"))
        self.assertLess(ke.region_coverage_fraction(1050, 1100, region), 0.10)
        self.assertGreaterEqual(ke.region_coverage_fraction(1000, 1599, region), 0.50)

    def test_structural_acmg_class_controls_k6_promotion(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            knowledge = ke.GeneKnowledge(db)
            structural = root / "sv.tsv"
            structural.write_text(
                "SV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ID\tGene_name\tACMG_class\n"
                "2\t1000\t1100\tDEL\tbenign\tTESTAD\t1\n"
                "2\t2000\t2100\tDEL\tvus\tTESTAD\t3\n"
                "2\t3000\t3100\tDEL\tpath\tTESTAD\t5\n",
                encoding="utf-8",
            )
            events, _ = ke.scan_structural(structural, "sv", knowledge)
            benign = next(ev for ev in events.values() if "benign" in ev.raw_context.get("structural_cluster_id", ""))
            self.assertFalse(benign.open_reasons)
            vus = next(ev for ev in events.values() if "vus" in ev.raw_context.get("structural_cluster_id", ""))
            path = next(ev for ev in events.values() if "path" in ev.raw_context.get("structural_cluster_id", ""))
            self.assertIn("K6_STRUCTURAL_GENE_CANDIDATE", vus.open_reasons)
            self.assertNotIn("K6_HIGH_EVIDENCE_SV_GENE", vus.open_reasons)
            self.assertIn("K6_HIGH_EVIDENCE_SV_GENE", path.open_reasons)

    def test_structural_split_rows_and_reciprocal_manta_mates_collapse(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            knowledge = ke.GeneKnowledge(db)
            structural = root / "sv.tsv"
            structural.write_text(
                "SV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ID\tID\tINFO\tFORMAT\tGFX0242962\tGene_name\tACMG_class\tAnnotation_mode\tFILTER\tQUAL\n"
                "2\t1000\t5000000\tDEL\ta\tMantaBND:1:0\tSVTYPE=BND;MATEID=MantaBND:1:1\tGT:FT:GQ:PR:SR\t0/1:PASS:99:10,8:10,8\tTESTAD\t5\tfull\tPASS\t200\n"
                "2\t1000\t5000000\tDEL\ta\tMantaBND:1:0\tSVTYPE=BND;MATEID=MantaBND:1:1\tGT:FT:GQ:PR:SR\t0/1:PASS:99:10,8:10,8\tTESTAD\tfull=5\tsplit\tPASS\t200\n"
                "2\t1000\t5000012\tDEL\tb\tMantaBND:1:1\tSVTYPE=BND;MATEID=MantaBND:1:0\tGT:FT:GQ:PR:SR\t0/1:PASS:99:10,8:10,8\tTESTAD\t5\tfull\tPASS\t200\n",
                encoding="utf-8",
            )
            events, counts = ke.scan_structural(structural, "sv", knowledge)
            self.assertEqual(counts["split_rows_skipped"], 1)
            self.assertEqual(counts["reciprocal_or_duplicate_full_rows_collapsed"], 1)
            self.assertEqual(len(events), 1)

    def test_large_unconfirmed_manta_del_is_demoted_but_cross_caller_is_confirmed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            knowledge = ke.GeneKnowledge(db)
            sv_path = root / "sv.tsv"
            cnv_path = root / "cnv.tsv"
            header = "SV_chrom\tSV_start\tSV_end\tSV_type\tAnnotSV_ID\tID\tINFO\tFORMAT\tGFX0242962\tGene_name\tACMG_class\tAnnotation_mode\tFILTER\tQUAL\n"
            sv_path.write_text(
                header +
                "2\t1000\t5001000\tDEL\tmanta\tMantaDEL:1\tEND=5001000;SVTYPE=DEL\tGT:FT:GQ:PR:SR\t0/1:PASS:99:10,8:10,8\tTESTAD\t5\tfull\tPASS\t200\n",
                encoding="utf-8",
            )
            cnv_path.write_text(
                header +
                "2\t1200\t5000800\tDEL\tdragen\tDRAGEN:LOSS:2\tSVTYPE=CNV;END=5000800\tGT:SM:CN:BC:PE\t0/1:0.5:1:100:10,10\tTESTAD\t5\tfull\tPASS\t80\n",
                encoding="utf-8",
            )
            sv_events, _ = ke.scan_structural(sv_path, "SV", knowledge)
            cnv_events, _ = ke.scan_structural(cnv_path, "CNV", knowledge)
            # Without independent support, the multi-megabase Manta call is not high-confidence.
            lone, counts = ke.merge_cross_caller_structural({}, dict(sv_events))
            only = next(iter(lone.values()))
            self.assertTrue(ke.event_technical_confidence(only).startswith("CAUTION_LARGE_UNCONFIRMED"))
            self.assertNotIn("K6_HIGH_EVIDENCE_SV_GENE", only.open_reasons)

            sv_events, _ = ke.scan_structural(sv_path, "SV", knowledge)
            merged, counts = ke.merge_cross_caller_structural(cnv_events, sv_events)
            only = next(iter(merged.values()))
            self.assertEqual(ke.event_technical_confidence(only), "CROSS_CALLER_CONFIRMED")
            self.assertEqual(counts["cross_caller_merged"], 1)

    def test_biallelic_flag_is_gene_scoped_and_requires_rare_variants(self) -> None:
        class MiniKnowledge:
            def __init__(self) -> None:
                self.ar = ke.GeneProfile("ARGENE", "HGNC:20")
                self.ar.relations = [ke.GeneDisease("ARGENE", "HGNC:20", "MONDO:20", "AR disease", "Definitive", "Autosomal recessive", "GenCC")]
            def profile(self, symbol: str):
                return self.ar if symbol == "ARGENE" else None

        e1 = ke.PatientEvent("SNV_INDEL", "v1", {"DANTE"}, "1", 1, 1, zygosity="HET_ALT", genes={"ARGENE", "OTHER"}, primary_gene="ARGENE", consequence="missense_variant", max_af="0.0001")
        e2 = ke.PatientEvent("SNV_INDEL", "v2", {"DANTE"}, "1", 2, 2, zygosity="HET_ALT", genes={"ARGENE"}, primary_gene="ARGENE", consequence="missense_variant", max_af="0.0002")
        common = ke.PatientEvent("SNV_INDEL", "v3", {"DANTE"}, "1", 3, 3, zygosity="HET_ALT", genes={"ARGENE"}, primary_gene="ARGENE", consequence="missense_variant", max_af="0.25")
        events = {e.event_key: e for e in (e1, e2, common)}
        self.assertEqual(ke.mark_possible_biallelic(events, MiniKnowledge(), 0.01), 2)
        self.assertEqual(e1.biallelic_genes, {"ARGENE"})
        self.assertEqual(e2.biallelic_genes, {"ARGENE"})
        self.assertFalse(common.possible_biallelic_unphased)

    def test_region_writer_supports_chr_prefixed_tabix_resources(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "regions.tsv"
            event = ke.PatientEvent("SNV_INDEL", "v", {"DANTE"}, "X", 123, 123, ref="A", alt="T")
            self.assertEqual(ke.write_regions([event], path, only_snv=True, chrom_prefix="chr"), 1)
            self.assertEqual(path.read_text(encoding="ascii"), "chrX\t123\t123\n")

    def test_research_usage_manifest_does_not_relabel_open_source_as_research_only(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_rows = [{"evidence_sources": "ClinVar;GIAB", "disease_id": "."}]
            research_rows = [{"evidence_sources": "GIAB;PanelApp"}]
            _, research_path = ke.write_usage_manifests(root, open_rows, research_rows)
            data = json.loads(research_path.read_text(encoding="utf-8"))
            records = data["sources_used"]
            giab = [r for r in records if r["source"] == "GIAB"]
            panel = [r for r in records if r["source"] == "PanelApp"]
            self.assertEqual(len(giab), 1)
            self.assertEqual(giab[0]["publication_scope"], "COMMERCIAL_ALLOWED")
            self.assertEqual(len(panel), 1)
            self.assertEqual(panel[0]["publication_scope"], "RESEARCH_ONLY")


    def test_missing_af_is_unknown_not_rare(self) -> None:
        self.assertEqual(ke.rarity_status(".", 0.01), "UNKNOWN")
        self.assertEqual(ke.rarity_status("", 0.01), "UNKNOWN")
        self.assertFalse(ke.is_rare(".", 0.01))
        self.assertTrue(ke.is_rare("0.001", 0.01))
        self.assertFalse(ke.is_rare("0.1", 0.01))

    def test_clinvar_mixed_pathogenic_protective_is_not_clean_plp(self) -> None:
        self.assertEqual(ke.clinvar_bucket("Likely_pathogenic|protective"), "MIXED")
        self.assertEqual(ke.clinvar_bucket("Pathogenic|risk_factor"), "MIXED")
        self.assertEqual(ke.clinvar_bucket("Pathogenic/Likely_pathogenic"), "PATHOGENIC")

    def test_primary_assembly_contig_wins_over_alt_scaffolds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            report = Path(directory) / "assembly_report.txt"
            report.write_text(
                "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
                "1\tassembled-molecule\t1\tChromosome\tCM000663.1\t=\tNC_000001.10\tPrimary Assembly\t249250621\tchr1\n"
                "HSCHR1_1_CTG31\talt-scaffold\t1\tChromosome\tGL000191.1\t<>\tNT_167207.1\tALT_REF_LOCI_1\t106433\tna\n"
                "X\tassembled-molecule\tX\tChromosome\tCM000685.1\t=\tNC_000023.10\tPrimary Assembly\t155270560\tchrX\n",
                encoding="utf-8",
            )
            mapping = ke.parse_assembly_report(report)
            self.assertEqual(mapping["1"], "NC_000001.10")
            self.assertEqual(mapping["X"], "NC_000023.10")

    def test_exact_clinvar_benign_blocks_gene_level_promotion(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            knowledge = ke.GeneKnowledge(db)
            catalog = root / "catalog.tsv.gz"
            fields = [
                "variant_key", "chrom", "position", "ref", "alt", "alt_carried", "vcf_filter",
                "source_layer", "GT", "symbol", "all_symbols", "consequence", "all_consequences", "impact",
                "max_af", "rsid", "clinvar_exact_match", "clinvar_current_rsid",
                "clinvar_current_significance", "clinvar_current_review_status", "clinvar_review_stars",
                "clinvar_current_conditions", "clinvar_current_geneinfo", "hgvsc", "hgvsp", "sift", "polyphen", "DP", "GQ", "AD",
            ]
            row = [
                "GRCh37:2:500:C:T", "2", "500", "C", "T", "yes", "PASS", "DANTE", "0/1",
                "TESTAD", "TESTAD", "missense_variant", "missense_variant", "MODERATE", "0.0001", "rs500",
                "yes", "rs500", "Benign", "criteria_provided,_single_submitter", "1", "Synthetic", "TESTAD:2",
                ".", ".", ".", ".", "30", "99", "15,15",
            ]
            with gzip.open(catalog, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(fields)
                writer.writerow(row)
            events, counts = ke.scan_small_catalog(catalog, knowledge, 0.01, progress_every=0)
            self.assertEqual(counts["catalog_rows"], 1)
            self.assertEqual(events, {})

    def test_selected_gene_consequence_does_not_inherit_other_transcript_lof(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            open_root, research_root = make_bundles(root)
            db = ke.build_reference_index(open_root, research_root, root / "index")
            knowledge = ke.GeneKnowledge(db)
            catalog = root / "catalog.tsv.gz"
            fields = [
                "variant_key", "chrom", "position", "ref", "alt", "alt_carried", "vcf_filter",
                "source_layer", "GT", "symbol", "all_symbols", "consequence", "all_consequences", "impact",
                "max_af", "rsid", "clinvar_exact_match", "clinvar_current_rsid",
                "clinvar_current_significance", "clinvar_current_review_status", "clinvar_review_stars",
                "clinvar_current_conditions", "clinvar_current_geneinfo", "hgvsc", "hgvsp", "sift", "polyphen", "DP", "GQ", "AD",
            ]
            row = [
                "GRCh37:2:600:C:T", "2", "600", "C", "T", "yes", "PASS", "DANTE", "0/1",
                "TESTAD", "TESTAD;OTHER", "synonymous_variant", "synonymous_variant;stop_gained", "LOW", "0.0001", "rs600",
                "no", ".", ".", ".", "0", ".", ".", ".", ".", ".", ".", "30", "99", "15,15",
            ]
            with gzip.open(catalog, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(fields)
                writer.writerow(row)
            events, _ = ke.scan_small_catalog(catalog, knowledge, 0.01, progress_every=0)
            self.assertEqual(events, {})

    def test_cadd_tabix_failure_uses_atomic_api_fallback(self) -> None:
        event = ke.PatientEvent(
            "SNV_INDEL", "GRCh37:2:300:C:T", {"DANTE"}, "2", 300, 300,
            ref="C", alt="T", consequence="missense_variant", max_af="0.001",
        )
        event.research_reasons.add("K8_CADD_SUPPORT_TARGET")
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            with mock.patch.object(ke, "docker_tabix", side_effect=RuntimeError("http2 failure")), \
                 mock.patch.object(ke, "_cadd_api_lookup", return_value=({event.event_key: ("1.5", "20.0")}, collections.Counter({"api_requests_ok": 1, "matched": 1}))):
                counts = ke.apply_cadd({event.event_key: event}, work, {}, "image", work)
            self.assertEqual(event.cadd_raw, "1.5")
            self.assertEqual(event.cadd_phred, "20.0")
            self.assertEqual(counts["transport_api_fallback"], 1)

    def test_powershell_ready_contract_matches_python_rules_schema(self) -> None:
        ps1 = (ROOT / "src" / "dante_wgs_pipeline.ps1").read_text(encoding="utf-8")
        self.assertIn(f'$state.rules_schema -ne "{ke.RULES_SCHEMA}"', ps1)

    def test_stage25_deduplicates_disease_rows_by_event(self) -> None:
        rows = [
            {"event_key": "v1", "priority_rank": "3", "gene": "A", "position": "10", "finding_id": "f1", "disease_id": "MONDO:1"},
            {"event_key": "v1", "priority_rank": "3", "gene": "A", "position": "10", "finding_id": "f2", "disease_id": "MONDO:2"},
            {"event_key": "v2", "priority_rank": "1", "gene": "B", "position": "20"},
        ]
        chosen = vra.prioritize_unique_rows(rows)
        self.assertEqual([row["event_key"] for row in chosen], ["v2", "v1"])
        v1 = next(row for row in chosen if row["event_key"] == "v1")
        self.assertEqual(v1["stage24_finding_ids"], "f1;f2")
        self.assertIn("MONDO:1", v1["stage24_relation_context_json"])
        self.assertIn("MONDO:2", v1["stage24_relation_context_json"])


if __name__ == "__main__":
    unittest.main()
