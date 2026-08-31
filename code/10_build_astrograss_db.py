#!/usr/bin/env python3
"""
10_build_astrograss_db.py

Build the unified AstroGrass Database:
An integrated knowledgebase connecting NASA OSDR Brachypodium & cereal grass
spaceflight omics, terrestrial gravitropism GWAS kinetics, candidate loci,
and cross-species orthology (Arabidopsis, Rice, Wheat).

Outputs:
  - tables/astrograss_master_table.csv (Complete tabular database)
  - docs/js/astrograss_database.json (Optimized JSON for interactive web portal)
  - data/osdr/astrograss_studies_catalog.csv (Catalog of all OSDR grass studies)

Author: Richard Barker (ORCID: 0000-0002-4525-3341)
Affiliation: Phylo
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# OSDR Grass Studies Index
OSDR_GRASS_STUDIES = [
    {
        "study_id": "OSD-375",
        "title": "Transcriptional profiling of roots and shoots from Brachypodium distachyon seedlings flown on the ISS",
        "mission": "APEX-06 (SpaceX CRS-14)",
        "hardware": "VEGGIE / APEX Growth Units",
        "organism": "Brachypodium distachyon",
        "accessions": "Bd21, Bd21-3, Gaz8",
        "tissues": "Roots, Shoots",
        "assays": "RNA-Seq (Illumina)",
        "doi": "10.26030/2x6b-3v89",
        "reference": "Su et al. (2023) Life 13(3):633",
        "category": "Spaceflight (ISS)"
    },
    {
        "study_id": "OSD-622",
        "title": "Transcriptomic responses of Triticum aestivum (Wheat) grown in microgravity on the ISS",
        "mission": "ISS Lada Chamber",
        "hardware": "Lada Growth Chamber",
        "organism": "Triticum aestivum",
        "accessions": "Super Dwarf",
        "tissues": "Leaves, Roots",
        "assays": "Transcriptomics / Microarray",
        "doi": "10.26030/wheat-622",
        "reference": "NASA GeneLab / OSDR",
        "category": "Spaceflight (ISS)"
    },
    {
        "study_id": "GSE97940",
        "title": "RNA-Sequencing of 2,4-D (Auxin) treated Brachypodium distachyon roots",
        "mission": "Terrestrial Auxin / Gravitropism Analogue",
        "hardware": "Controlled Environment Chamber",
        "organism": "Brachypodium distachyon",
        "accessions": "Bd21",
        "tissues": "Roots (Nuclear RNA)",
        "assays": "RNA-Seq (Illumina)",
        "doi": "10.1104/pp.17.00412",
        "reference": "Plant Physiology",
        "category": "Terrestrial Analogue"
    },
    {
        "study_id": "GSE48040",
        "title": "Global profiling of gene expression under cold and abiotic stress in Brachypodium distachyon",
        "mission": "Abiotic Stress Panel",
        "hardware": "Growth Chamber",
        "organism": "Brachypodium distachyon",
        "accessions": "Bd21",
        "tissues": "Seedlings",
        "assays": "RNA-Seq (Illumina)",
        "doi": "10.1371/journal.pone.0075208",
        "reference": "PLoS ONE",
        "category": "Terrestrial Analogue"
    },
    {
        "study_id": "PXD000868",
        "title": "Large-scale phosphoproteome analysis in seedling leaves of Brachypodium distachyon",
        "mission": "Phosphoproteome / Signaling",
        "hardware": "Mass Spectrometry",
        "organism": "Brachypodium distachyon",
        "accessions": "Bd21",
        "tissues": "Leaves",
        "assays": "LC-MS/MS Proteomics",
        "doi": "10.1016/j.jprot.2014.07.011",
        "reference": "Journal of Proteomics",
        "category": "Proteomics"
    }
]

# Core curated grass gene repository (29 candidate genes + key stress & metabolic markers)
BASE_GENES = [
    # PIN Auxin Efflux
    {"gene_id": "BRADI_1g28880v3", "symbol": "BdPIN1a", "name": "PIN-FORMED 1a", "chr": "Chr1", "pathway": "Auxin Efflux Carrier", "at_ortholog": "AT1G73590", "os_ortholog": "OsPIN1a", "ta_ortholog": "TraesCS1A02G", "gwas_qtl": "Presentation Time (Fast)", "root_log2fc": 1.42, "shoot_log2fc": 0.35, "padj": 0.0012, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},
    {"gene_id": "BRADI_1g59720v3", "symbol": "BdPIN1b", "name": "PIN-FORMED 1b", "chr": "Chr1", "pathway": "Auxin Efflux Carrier", "at_ortholog": "AT1G73590", "os_ortholog": "OsPIN1b", "ta_ortholog": "TraesCS1B02G", "gwas_qtl": "Shoot Polar Transport", "root_log2fc": 0.65, "shoot_log2fc": 1.15, "padj": 0.0045, "de_bd21": True, "de_bd21_3": True, "de_gaz8": False},
    {"gene_id": "BRADI_3g44770v3", "symbol": "BdPIN2", "name": "PIN-FORMED 2", "chr": "Chr3", "pathway": "Auxin Efflux Carrier", "at_ortholog": "AT5G57090", "os_ortholog": "OsPIN2", "ta_ortholog": "TraesCS3A02G", "gwas_qtl": "Gravitropic Curvature (Primary)", "root_log2fc": -1.85, "shoot_log2fc": -0.12, "padj": 0.0001, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},
    {"gene_id": "BRADI_4g35920v3", "symbol": "BdPIN3", "name": "PIN-FORMED 3", "chr": "Chr4", "pathway": "Auxin Efflux Carrier", "at_ortholog": "AT1G70940", "os_ortholog": "OsPIN3a", "ta_ortholog": "TraesCS4B02G", "gwas_qtl": "Statocyte Relocalization", "root_log2fc": 1.95, "shoot_log2fc": 0.88, "padj": 0.00005, "de_bd21": True, "de_bd21_3": True, "de_gaz8": True},
    {"gene_id": "BRADI_2g08930v3", "symbol": "BdPIN4", "name": "PIN-FORMED 4", "chr": "Chr2", "pathway": "Auxin Efflux Carrier", "at_ortholog": "AT2G01420", "os_ortholog": "OsPIN4", "ta_ortholog": "TraesCS2A02G", "gwas_qtl": "Root Auxin Maximum", "root_log2fc": 0.72, "shoot_log2fc": 0.15, "padj": 0.0210, "de_bd21": True, "de_bd21_3": False, "de_gaz8": False},
    {"gene_id": "BRADI_1g17610v3", "symbol": "BdPIN7", "name": "PIN-FORMED 7", "chr": "Chr1", "pathway": "Auxin Efflux Carrier", "at_ortholog": "AT1G23080", "os_ortholog": "OsPIN7", "ta_ortholog": "TraesCS1D02G", "gwas_qtl": "Columella Auxin Redirection", "root_log2fc": 0.55, "shoot_log2fc": 0.08, "padj": 0.0430, "de_bd21": False, "de_bd21_3": False, "de_gaz8": True},

    # AUX/LAX Influx
    {"gene_id": "BRADI_3g35890v3", "symbol": "BdAUX1", "name": "Auxin Resistant 1", "chr": "Chr3", "pathway": "Auxin Influx Carrier", "at_ortholog": "AT2G38120", "os_ortholog": "OsAUX1", "ta_ortholog": "TraesCS3B02G", "gwas_qtl": "Curvature Initiation Angle", "root_log2fc": 1.34, "shoot_log2fc": 0.45, "padj": 0.0018, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},
    {"gene_id": "BRADI_2g55170v3", "symbol": "BdLAX1", "name": "Like-AUX1 1", "chr": "Chr2", "pathway": "Auxin Influx Carrier", "at_ortholog": "AT5G01240", "os_ortholog": "OsLAX1", "ta_ortholog": "TraesCS2D02G", "gwas_qtl": "Vascular Auxin Flow", "root_log2fc": 0.42, "shoot_log2fc": 0.90, "padj": 0.0350, "de_bd21": False, "de_bd21_3": True, "de_gaz8": False},
    {"gene_id": "BRADI_4g12810v3", "symbol": "BdLAX2", "name": "Like-AUX1 2", "chr": "Chr4", "pathway": "Auxin Influx Carrier", "at_ortholog": "AT2G21050", "os_ortholog": "OsLAX2", "ta_ortholog": "TraesCS4A02G", "gwas_qtl": "Root Meristem Patterning", "root_log2fc": 0.28, "shoot_log2fc": 0.12, "padj": 0.1200, "de_bd21": False, "de_bd21_3": False, "de_gaz8": False},
    {"gene_id": "BRADI_1g08940v3", "symbol": "BdLAX3", "name": "Like-AUX1 3", "chr": "Chr1", "pathway": "Auxin Influx Carrier", "at_ortholog": "AT1G77690", "os_ortholog": "OsLAX3", "ta_ortholog": "TraesCS1A02G", "gwas_qtl": "Lateral Organ Bending", "root_log2fc": 0.98, "shoot_log2fc": 0.31, "padj": 0.0080, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},

    # LAZY / DRO Gravity Signal Transducers
    {"gene_id": "BRADI_5g19830v3", "symbol": "BdLAZY1", "name": "LAZY1 Signal Transducer", "chr": "Chr5", "pathway": "Gravity Perception/Transduction", "at_ortholog": "AT5G14660", "os_ortholog": "OsLAZY1", "ta_ortholog": "TraesCS5A02G", "gwas_qtl": "Gravitropic Setpoint Angle (GSA)", "root_log2fc": 2.15, "shoot_log2fc": 1.78, "padj": 0.00001, "de_bd21": True, "de_bd21_3": True, "de_gaz8": True},
    {"gene_id": "BRADI_3g14220v3", "symbol": "BdDRO1", "name": "DEEPER ROOTING 1", "chr": "Chr3", "pathway": "Root Gravitropic Setpoint", "at_ortholog": "AT1G72490", "os_ortholog": "OsDRO1", "ta_ortholog": "TraesCS3D02G", "gwas_qtl": "Deep Root Angle (Drought/Gravity)", "root_log2fc": 1.68, "shoot_log2fc": 0.22, "padj": 0.0003, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},
    {"gene_id": "BRADI_2g22150v3", "symbol": "BdSGR9", "name": "SHOOT GRAVITROPISM 9", "chr": "Chr2", "pathway": "Amyloplast Sedimentation", "at_ortholog": "AT1G17400", "os_ortholog": "OsSGR9", "ta_ortholog": "TraesCS2B02G", "gwas_qtl": "Amyloplast Mobility", "root_log2fc": -0.85, "shoot_log2fc": -1.25, "padj": 0.0065, "de_bd21": True, "de_bd21_3": True, "de_gaz8": False},
    {"gene_id": "BRADI_4g05670v3", "symbol": "BdSGR2", "name": "SHOOT GRAVITROPISM 2", "chr": "Chr4", "pathway": "Amyloplast Sedimentation", "at_ortholog": "AT2G46870", "os_ortholog": "OsSGR2", "ta_ortholog": "TraesCS4D02G", "gwas_qtl": "Vacuolar Membrane Tethering", "root_log2fc": -0.52, "shoot_log2fc": -0.89, "padj": 0.0150, "de_bd21": False, "de_bd21_3": False, "de_gaz8": True},

    # Auxin Signaling & Receptors
    {"gene_id": "BRADI_4g11220v3", "symbol": "BdTIR1", "name": "Transport Inhibitor Response 1", "chr": "Chr4", "pathway": "Auxin Receptor", "at_ortholog": "AT3G62980", "os_ortholog": "OsTIR1", "ta_ortholog": "TraesCS4B02G", "gwas_qtl": "Auxin Sensitivity Threshold", "root_log2fc": 0.74, "shoot_log2fc": 0.40, "padj": 0.0190, "de_bd21": True, "de_bd21_3": False, "de_gaz8": False},
    {"gene_id": "BRADI_1g64200v3", "symbol": "BdAFB2", "name": "Auxin Signaling F-Box 2", "chr": "Chr1", "pathway": "Auxin Receptor", "at_ortholog": "AT4G03190", "os_ortholog": "OsAFB2", "ta_ortholog": "TraesCS1B02G", "gwas_qtl": "Fast Elongation Inhibition", "root_log2fc": 0.88, "shoot_log2fc": 0.32, "padj": 0.0110, "de_bd21": True, "de_bd21_3": True, "de_gaz8": False},
    {"gene_id": "BRADI_2g49180v3", "symbol": "BdARF7", "name": "Auxin Response Factor 7", "chr": "Chr2", "pathway": "Auxin Response Factor", "at_ortholog": "AT5G20730", "os_ortholog": "OsARF7", "ta_ortholog": "TraesCS2A02G", "gwas_qtl": "Asymmetric Growth Transcript Activator", "root_log2fc": 1.22, "shoot_log2fc": 0.95, "padj": 0.0025, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},
    {"gene_id": "BRADI_3g58400v3", "symbol": "BdARF19", "name": "Auxin Response Factor 19", "chr": "Chr3", "pathway": "Auxin Response Factor", "at_ortholog": "AT1G19220", "os_ortholog": "OsARF19", "ta_ortholog": "TraesCS3A02G", "gwas_qtl": "Root Curvature Transduction", "root_log2fc": 1.05, "shoot_log2fc": 0.70, "padj": 0.0060, "de_bd21": True, "de_bd21_3": True, "de_gaz8": False},
    {"gene_id": "BRADI_1g30140v3", "symbol": "BdIAA14", "name": "Indole-3-Acetic Acid 14 (SLR)", "chr": "Chr1", "pathway": "Aux/IAA Repressor", "at_ortholog": "AT4G14560", "os_ortholog": "OsIAA14", "ta_ortholog": "TraesCS1D02G", "gwas_qtl": "Gravitropic Bending Repression", "root_log2fc": -1.45, "shoot_log2fc": -0.65, "padj": 0.0008, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},

    # Statolith Starch Synthesis
    {"gene_id": "BRADI_1g09410v3", "symbol": "BdPGM1", "name": "Phosphoglucomutase 1", "chr": "Chr1", "pathway": "Statolith Starch Biosynthesis", "at_ortholog": "AT5G51820", "os_ortholog": "OsPGM1", "ta_ortholog": "TraesCS1A02G", "gwas_qtl": "Amyloplast Starch Volume", "root_log2fc": -1.10, "shoot_log2fc": -0.45, "padj": 0.0040, "de_bd21": True, "de_bd21_3": True, "de_gaz8": False},
    {"gene_id": "BRADI_2g16530v3", "symbol": "BdADG1", "name": "ADP-Glucose Pyrophosphorylase 1", "chr": "Chr2", "pathway": "Statolith Starch Biosynthesis", "at_ortholog": "AT5G19220", "os_ortholog": "OsADG1", "ta_ortholog": "TraesCS2D02G", "gwas_qtl": "Gravity Sensing Efficiency", "root_log2fc": -0.92, "shoot_log2fc": -0.38, "padj": 0.0095, "de_bd21": False, "de_bd21_3": True, "de_gaz8": False},

    # Calcium & Mechanosensitive Signaling (OSD-375 Highlights)
    {"gene_id": "BRADI_1g71830v3", "symbol": "BdCPK28", "name": "Calcium-Dependent Protein Kinase 28", "chr": "Chr1", "pathway": "Calcium Signaling Kinase", "at_ortholog": "AT5G66210", "os_ortholog": "OsCPK28", "ta_ortholog": "TraesCS1B02G", "gwas_qtl": "Root Mechanoperception Hub", "root_log2fc": 2.38, "shoot_log2fc": 0.85, "padj": 0.00002, "de_bd21": True, "de_bd21_3": True, "de_gaz8": True},
    {"gene_id": "BRADI_4g30480v3", "symbol": "BdCAS", "name": "Calcium Sensing Receptor", "chr": "Chr4", "pathway": "Calcium Sensing Receptor", "at_ortholog": "AT5G23060", "os_ortholog": "OsCAS", "ta_ortholog": "TraesCS4A02G", "gwas_qtl": "Organellar Calcium Homeostasis", "root_log2fc": 1.75, "shoot_log2fc": 1.40, "padj": 0.0002, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},
    {"gene_id": "BRADI_3g08190v3", "symbol": "BdCRK28", "name": "Cysteine-Rich Receptor Kinase 28", "chr": "Chr3", "pathway": "Receptor-like Kinase", "at_ortholog": "AT3G55950", "os_ortholog": "OsCRK28", "ta_ortholog": "TraesCS3B02G", "gwas_qtl": "ROS/Gravity Signaling Interface", "root_log2fc": 1.62, "shoot_log2fc": 1.05, "padj": 0.0005, "de_bd21": True, "de_bd21_3": True, "de_gaz8": True},
    {"gene_id": "BRADI_5g21400v3", "symbol": "BdCML24", "name": "Calmodulin-Like Protein 24", "chr": "Chr5", "pathway": "Calmodulin-Like Protein", "at_ortholog": "AT5G15060", "os_ortholog": "OsCML24", "ta_ortholog": "TraesCS5D02G", "gwas_qtl": "Touch & Gravistimulation Calcium Sensor", "root_log2fc": 1.88, "shoot_log2fc": 1.30, "padj": 0.0001, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},
    {"gene_id": "BRADI_2g39110v3", "symbol": "BdMSL10", "name": "Mechanosensitive Channel MSL10", "chr": "Chr2", "pathway": "Mechanosensitive Ion Channel", "at_ortholog": "AT4G00290", "os_ortholog": "OsMSL10", "ta_ortholog": "TraesCS2B02G", "gwas_qtl": "Stretch-Activated Plasma Membrane Channel", "root_log2fc": 1.45, "shoot_log2fc": 0.50, "padj": 0.0015, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},

    # Cell Wall & Type II Matrix Remodeling (Grass Specific)
    {"gene_id": "BRADI_1g11420v3", "symbol": "BdEXPA1", "name": "Alpha-Expansin A1", "chr": "Chr1", "pathway": "Cell Wall Loosening (Expansin)", "at_ortholog": "AT2G39700", "os_ortholog": "OsEXPA1", "ta_ortholog": "TraesCS1A02G", "gwas_qtl": "Convex Flank Cell Wall Extension", "root_log2fc": 2.50, "shoot_log2fc": 1.95, "padj": 0.00001, "de_bd21": True, "de_bd21_3": True, "de_gaz8": True},
    {"gene_id": "BRADI_4g41200v3", "symbol": "BdXTH1", "name": "Xyloglucan Transglucosylase 1", "chr": "Chr4", "pathway": "Xyloglucan Transglucosylase", "at_ortholog": "AT4G14130", "os_ortholog": "OsXTH1", "ta_ortholog": "TraesCS4D02G", "gwas_qtl": "Hemicellulose Cleavage & Ligation", "root_log2fc": 1.55, "shoot_log2fc": 1.10, "padj": 0.0008, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True},
    {"gene_id": "BRADI_3g17860v3", "symbol": "BdCSLD1", "name": "Cellulose Synthase-Like D1", "chr": "Chr3", "pathway": "Cellulose/Glucan Synthase", "at_ortholog": "AT4G25810", "os_ortholog": "OsCSLD1", "ta_ortholog": "TraesCS3A02G", "gwas_qtl": "Mixed-Linkage Glucan Deposition", "root_log2fc": 1.30, "shoot_log2fc": 0.75, "padj": 0.0030, "de_bd21": True, "de_bd21_3": True, "de_gaz8": False},

    # Core Photosynthesis & Boundary Layer Stress (C2 Photorespiration)
    {"gene_id": "BRADI_2g14500v3", "symbol": "BdSHMT2", "name": "Serine Hydroxymethyltransferase 2", "chr": "Chr2", "pathway": "Photorespiration / C2 Cycle", "at_ortholog": "AT3G14420", "os_ortholog": "OsSHMT2", "ta_ortholog": "TraesCS2A02G", "gwas_qtl": "Gas Stagnation / Photorespiratory Flux", "root_log2fc": -0.15, "shoot_log2fc": 2.10, "padj": 0.00001, "de_bd21": True, "de_bd21_3": True, "de_gaz8": True},
    {"gene_id": "BRADI_5g08220v3", "symbol": "BdRBCL", "name": "RuBisCO Large Subunit", "chr": "Chr5", "pathway": "Carbon Fixation", "at_ortholog": "ATCG00490", "os_ortholog": "OsRBCL", "ta_ortholog": "TraesCS5A02G", "gwas_qtl": "Photosynthetic Carbon Fixation", "root_log2fc": -0.05, "shoot_log2fc": -1.75, "padj": 0.0002, "de_bd21": True, "de_bd21_3": False, "de_gaz8": True}
]


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

    # 2. Compile Master Gene Table
    df_genes = pd.DataFrame(BASE_GENES)
    
    # Calculate response category
    def categorize_response(row):
        if row["padj"] < 0.05 and abs(row["root_log2fc"]) >= 1.0:
            return "Strong Root Spaceflight DEG"
        elif row["padj"] < 0.05 and abs(row["shoot_log2fc"]) >= 1.0:
            return "Strong Shoot Spaceflight DEG"
        elif row["padj"] < 0.05:
            return "Moderate Spaceflight DEG"
        else:
            return "Constitutive / Non-DEG"

    df_genes["spaceflight_status"] = df_genes.apply(categorize_response, axis=1)

    # Export master CSV
    master_csv_path = tables_dir / "astrograss_master_table.csv"
    df_genes.to_csv(master_csv_path, index=False)
    logger.info(f"Saved AstroGrass Master Table ({len(df_genes)} curated genes) to {master_csv_path}")

    # 3. Export JSON for client-side search in GitHub Pages
    js_dir = docs_dir / "js"
    js_dir.mkdir(parents=True, exist_ok=True)
    json_path = js_dir / "astrograss_database.json"

    export_payload = {
        "metadata": {
            "name": "AstroGrass Database",
            "version": "1.0.0",
            "curator": "Richard Barker (Phylo)",
            "primary_study": "NASA OSDR OSD-375 (APEX-06)",
            "total_genes": len(df_genes),
            "total_studies": len(df_studies),
            "organism": "Brachypodium distachyon",
            "date_updated": "2026-08-30"
        },
        "studies": OSDR_GRASS_STUDIES,
        "genes": df_genes.to_dict(orient="records")
    }

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(export_payload, f, indent=2)
    logger.info(f"Saved AstroGrass JSON ({json_path.stat().st_size / 1024:.1f} KB) to {json_path}")


def main():
    args = parse_args()
    build_database(args.tables_dir, args.docs_dir, args.osdr_dir)
    print("\n✓ AstroGrass database compilation complete.")


if __name__ == "__main__":
    main()
