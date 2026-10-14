#!/usr/bin/env python3
"""
10_build_astrograss_db.py

Build the AstroGrass database: a catalogue of grass spaceflight and related omics studies,
and the gravitropism candidate-gene list (data/genotypes/gravitropism_candidate_genes.csv,
written by 03_candidate_genes.py) with the candidates' spaceflight DEG calls from Su et al. (2023)
(tables/candidate_spaceflight_degs.csv, written by 16_su2023_degs.py).

No differential-expression analysis of OSD-375 is run here: every expression value is Su et al.'s
(Supplemental Tables S2-S4, CC BY 4.0), and a candidate has values only where Su et al. list one of
its Brachypodium genes as a DEG. The fold changes, adjusted p-values, per-accession DEG calls,
chromosome positions and "GWAS QTL roles" in versions before 2026-10-14 were typed-in values and
were removed.

Outputs:
  - tables/astrograss_master_table.csv (candidate-gene table)
  - docs/js/astrograss_database.json (JSON for the web portal)
  - data/osdr/astrograss_studies_catalog.csv (study catalogue)

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
Affiliation: Collaborative Science Environment Inc.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Studies catalogued in AstroGrass. Each field is taken from the study's own record (OSDR,
# GEO or PRIDE) or from the paper that record links to; fields a record does not state are
# marked "Not stated in record". Checked 2026-10-14.
OSDR_GRASS_STUDIES = [
    {
        "study_id": "OSD-375",
        "title": "Transcriptional profiling of roots and shoots from Brachypodium distachyon seedlings flown on the ISS",
        "mission": "SpaceX-14 (launched 2 April 2018); APEX-06",
        "hardware": "APEX Growth Units in VEGGIE",
        "organism": "Brachypodium distachyon",
        "accessions": "Bd21, Bd21-3, Gaz8",
        "tissues": "Roots, Shoots",
        "assays": "RNA-Seq (Illumina NovaSeq 6000)",
        "doi": "10.26030/2x6b-3v89",
        "reference": "Su et al. (2023) Life 13(3):626, doi:10.3390/life13030626",
        "category": "Spaceflight (ISS)"
    },
    {
        "study_id": "OSD-622",
        "title": "Microgravity effects on leaf morphology, cell structure, carbon metabolism and mRNA expression of dwarf wheat",
        "mission": "STS-110 to the ISS (PESTO; launched 8 April 2002, returned on STS-111)",
        "hardware": "Biomass Production System (BPS), PESTO root modules",
        "organism": "Triticum aestivum",
        "accessions": "USU Apogee",
        "tissues": "Leaves (pooled for the microarray)",
        "assays": "DNA microarray (GenePix)",
        "doi": "10.26030/m2eg-jx23",
        "reference": "Stutte et al. (2006) Planta 224:1038-1049, doi:10.1007/s00425-006-0290-4",
        "category": "Spaceflight (ISS)"
    },
    {
        "study_id": "GSE97940",
        "title": "RNA-Sequencing of 2,4-D treated Brachypodium distachyon roots containing NLS",
        "mission": "Ground study: 2,4-D-induced nodule-like structures (NLS)",
        "hardware": "Not stated in record",
        "organism": "Brachypodium distachyon",
        "accessions": "Not stated in record",
        "tissues": "Roots with nodule-like structures; untreated roots",
        "assays": "RNA-Seq (Illumina HiSeq 2500)",
        "doi": "10.1007/s10142-018-0594-z",
        "reference": "Thomas et al. (2018) Funct Integr Genomics 18:315-326",
        "category": "Ground study (auxin treatment)"
    },
    {
        "study_id": "GSE48040",
        "title": "A global profiling of gene expression under cold stress in Brachypodium",
        "mission": "Ground study: cold treatment (4 degC, 24 h)",
        "hardware": "Not stated in record",
        "organism": "Brachypodium distachyon",
        "accessions": "ABR5",
        "tissues": "Seedlings (12 days old)",
        "assays": "RNA-Seq (Illumina HiSeq 2000)",
        "doi": "10.1186/gb-2013-14-8-r92",
        "reference": "Zhang et al. (2013) Genome Biol 14:R92",
        "category": "Ground study (cold stress)"
    },
    {
        "study_id": "PXD000868",
        "title": "Large-scale phosphoproteome analysis in seedling leaves of Brachypodium distachyon L",
        "mission": "Ground study: phosphoproteome",
        "hardware": "Not stated in record",
        "organism": "Brachypodium distachyon",
        "accessions": "Not stated in record",
        "tissues": "Seedling leaves",
        "assays": "LC-MS/MS phosphoproteomics",
        "doi": "10.1186/1471-2164-15-375",
        "reference": "Lv et al. (2014) BMC Genomics 15:375",
        "category": "Proteomics"
    }
]

CANDIDATES = Path("data/genotypes/gravitropism_candidate_genes.csv")
CANDIDATE_DEGS = Path("tables/candidate_spaceflight_degs.csv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build AstroGrass master database.")
    parser.add_argument("--tables-dir", type=Path, default=Path("tables"), help="Path to tables directory")
    parser.add_argument("--docs-dir", type=Path, default=Path("docs"), help="Path to docs directory")
    parser.add_argument("--osdr-dir", type=Path, default=Path("data/osdr"), help="Path to OSDR directory")
    return parser.parse_args()


def build_database(tables_dir: Path, docs_dir: Path, osdr_dir: Path):
    """Compile and export the AstroGrass database."""
    logger.info("Building AstroGrass Master Knowledgebase...")
    tables_dir.mkdir(parents=True, exist_ok=True)
    docs_dir.mkdir(parents=True, exist_ok=True)
    osdr_dir.mkdir(parents=True, exist_ok=True)

    # 1. Export Studies Catalog
    df_studies = pd.DataFrame(OSDR_GRASS_STUDIES)
    studies_path = osdr_dir / "astrograss_studies_catalog.csv"
    df_studies.to_csv(studies_path, index=False)
    logger.info(f"Saved OSDR Studies Catalog ({len(df_studies)} studies) to {studies_path}")

    # 2. Candidate-gene table, with Su et al.'s DEG calls for each candidate
    df_genes = pd.read_csv(CANDIDATES, keep_default_na=False)
    df_degs = pd.read_csv(CANDIDATE_DEGS, keep_default_na=False)
    deg_records = {c: g.drop(columns=["candidate", "arabidopsis_gene", "arabidopsis_agi", "source"]).to_dict(orient="records")
                   for c, g in df_degs.groupby("candidate", sort=False)}

    def summarise(candidate: str) -> str:
        return "; ".join(
            f"{d['accession']} {d['organ']} {d['deseq_log2fc']:+.2f} ({d['brachypodium_gene_id']}, "
            f"{'Compara' if d['match'].startswith('Compara') else 'Su et al. best hit'})"
            for d in deg_records.get(candidate, []))

    df_genes["su2023_deg_calls"] = df_genes["candidate"].map(lambda c: len(deg_records.get(c, [])))
    df_genes["su2023_degs"] = df_genes["candidate"].map(summarise)

    # Export master CSV
    master_csv_path = tables_dir / "astrograss_master_table.csv"
    df_genes.to_csv(master_csv_path, index=False)
    logger.info(f"Saved AstroGrass Master Table ({len(df_genes)} candidate genes) to {master_csv_path}")

    # 3. Export JSON for client-side search in GitHub Pages
    js_dir = docs_dir / "js"
    js_dir.mkdir(parents=True, exist_ok=True)
    json_path = js_dir / "astrograss_database.json"

    export_payload = {
        "metadata": {
            "name": "AstroGrass Database",
            "version": "2.1.0",
            "curator": "Richard Barker (Collaborative Science Environment Inc.)",
            "primary_study": "NASA OSDR OSD-375 (APEX-06)",
            "total_genes": len(df_genes),
            "verified_arabidopsis_anchors": int((df_genes["arabidopsis_agi"] != "").sum()),
            "brachypodium_gene_ids": "Orthologs of the Arabidopsis anchor in Ensembl Compara (plants); see data/genotypes/compara_orthologs.csv",
            "expression_values": ("DESeq log2 fold change (ISS / ground) and adjusted p-value from Su et al. (2023) Life 13:626, "
                                  "Supplemental Tables S2-S4 (CC BY 4.0), listed only where Su et al. report the gene as a DEG"),
            "total_studies": len(df_studies),
            "organism": "Brachypodium distachyon",
            "date_updated": "2026-10-14"
        },
        "studies": OSDR_GRASS_STUDIES,
        "genes": [{**g, "spaceflight_degs": deg_records.get(g["candidate"], [])}
                  for g in df_genes.to_dict(orient="records")]
    }

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(export_payload, f, indent=2)
    logger.info(f"Saved AstroGrass JSON ({json_path.stat().st_size / 1024:.1f} KB) to {json_path}")

    # Synchronize all tables to docs/tables/ for GitHub Pages downloads
    docs_tables = docs_dir / "tables"
    docs_tables.mkdir(parents=True, exist_ok=True)
    for csv_file in tables_dir.glob("*.csv"):
        (docs_tables / csv_file.name).write_bytes(csv_file.read_bytes())
    if studies_path.exists():
        (docs_tables / studies_path.name).write_bytes(studies_path.read_bytes())
    # Also sync sample_metadata.csv if present
    sample_meta = osdr_dir / "sample_metadata.csv"
    if sample_meta.exists():
        (docs_tables / sample_meta.name).write_bytes(sample_meta.read_bytes())
    logger.info(f"✓ Synchronized all CSV tables to {docs_tables} for GitHub Pages downloads")


def main():
    args = parse_args()
    build_database(args.tables_dir, args.docs_dir, args.osdr_dir)
    print("\n✓ AstroGrass database compilation complete.")


if __name__ == "__main__":
    main()

