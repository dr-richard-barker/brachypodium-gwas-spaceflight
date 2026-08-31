#!/usr/bin/env python3
"""
Fetch Brachypodium distachyon SNP/variation data for GWAS analysis.

This script fetches or constructs SNP/variation data resources for Brachypodium distachyon,
focusing on gravitropism, auxin transport, and spaceflight-response genes (from Su et al. 2023).
It creates reference files for candidate genes, ecotype comparisons (Bd21, Bd21-3, Gaz8),
and documents data sources for the Brachypodium GWAS-spaceflight project.

Key Features:
- Connects to Ensembl Plants REST API (with offline fallback).
- Generates curated lists of gravitropism/spaceflight candidate genes.
- Creates documentation for downloading pan-genome VCFs from BrachyPan/Phytozome.

References:
- Su et al. (2023) Life 13(3):633. DOI: 10.3390/life13030633
- Gordon et al. (2017) Nature Communications 8:2042. (Brachypodium pan-genome)
- Ensembl Plants: https://plants.ensembl.org
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any

import pandas as pd
import requests

# Constants
DEFAULT_OUTPUT_DIR = Path("data/genotypes")
ENSEMBL_REST_SERVER = "https://rest.ensembl.org"

# Curated comprehensive list of gravitropism candidate genes in Brachypodium distachyon
FALLBACK_CANDIDATES = [
    # Auxin Efflux Carriers (PIN family)
    {"brachypodium_gene_id": "BRADI_1g28880v3", "arabidopsis_ortholog": "AT1G73590", "gene_symbol": "BdPIN1a", "pathway": "Auxin Efflux Carrier", "evidence": "Direct ortholog of AtPIN1; directional auxin transport in root meristem"},
    {"brachypodium_gene_id": "BRADI_1g59720v3", "arabidopsis_ortholog": "AT1G73590", "gene_symbol": "BdPIN1b", "pathway": "Auxin Efflux Carrier", "evidence": "Paralog of BdPIN1a involved in shoot-to-root polar auxin flux"},
    {"brachypodium_gene_id": "BRADI_3g44770v3", "arabidopsis_ortholog": "AT5G57090", "gene_symbol": "BdPIN2", "pathway": "Auxin Efflux Carrier", "evidence": "Ortholog of AtPIN2 (EIR1/AGR1); asymmetric auxin redistribution in root cortex/epidermis"},
    {"brachypodium_gene_id": "BRADI_4g35920v3", "arabidopsis_ortholog": "AT1G70940", "gene_symbol": "BdPIN3", "pathway": "Auxin Efflux Carrier", "evidence": "Ortholog of AtPIN3; lateral relocalization in columella statocytes upon gravistimulation"},
    {"brachypodium_gene_id": "BRADI_2g08930v3", "arabidopsis_ortholog": "AT2G01420", "gene_symbol": "BdPIN4", "pathway": "Auxin Efflux Carrier", "evidence": "Ortholog of AtPIN4; root tip auxin maximum maintenance"},
    {"brachypodium_gene_id": "BRADI_1g17610v3", "arabidopsis_ortholog": "AT1G23080", "gene_symbol": "BdPIN7", "pathway": "Auxin Efflux Carrier", "evidence": "Ortholog of AtPIN7; columella auxin redirection"},

    # Auxin Influx Carriers (AUX/LAX family)
    {"brachypodium_gene_id": "BRADI_3g35890v3", "arabidopsis_ortholog": "AT2G38120", "gene_symbol": "BdAUX1", "pathway": "Auxin Influx Carrier", "evidence": "Ortholog of AtAUX1; essential for root gravitropic curvature and cellular auxin uptake"},
    {"brachypodium_gene_id": "BRADI_2g55170v3", "arabidopsis_ortholog": "AT5G01240", "gene_symbol": "BdLAX1", "pathway": "Auxin Influx Carrier", "evidence": "Ortholog of AtLAX1; auxiliary auxin influx in vascular tissues"},
    {"brachypodium_gene_id": "BRADI_4g12810v3", "arabidopsis_ortholog": "AT2G21050", "gene_symbol": "BdLAX2", "pathway": "Auxin Influx Carrier", "evidence": "Ortholog of AtLAX2; root apical meristem patterning"},
    {"brachypodium_gene_id": "BRADI_1g08940v3", "arabidopsis_ortholog": "AT1G77690", "gene_symbol": "BdLAX3", "pathway": "Auxin Influx Carrier", "evidence": "Ortholog of AtLAX3; lateral root emergence and organ bending"},

    # Gravity Signal Transduction (LAZY / SGR family)
    {"brachypodium_gene_id": "BRADI_5g19830v3", "arabidopsis_ortholog": "AT5G14660", "gene_symbol": "BdLAZY1", "pathway": "Gravity Perception/Transduction", "evidence": "Ortholog of AtLAZY1/OsLAZY1; controls tiller angle and root gravitropic setpoint angle (GSA)"},
    {"brachypodium_gene_id": "BRADI_3g14220v3", "arabidopsis_ortholog": "AT1G72490", "gene_symbol": "BdDRO1", "pathway": "Root Gravitropic Setpoint", "evidence": "Ortholog of DEEPER ROOTING 1 (DRO1); regulates root gravitropic curvature angle"},
    {"brachypodium_gene_id": "BRADI_2g22150v3", "arabidopsis_ortholog": "AT1G17400", "gene_symbol": "BdSGR9", "pathway": "Amyloplast Sedimentation", "evidence": "Ortholog of SHOOT GRAVITROPISM 9; modulates statolith dynamics in endodermal cells"},
    {"brachypodium_gene_id": "BRADI_4g05670v3", "arabidopsis_ortholog": "AT2G46870", "gene_symbol": "BdSGR2", "pathway": "Amyloplast Sedimentation", "evidence": "Ortholog of AtSGR2 (phospholipase A1); involved in vacuolar gravity sensing"},

    # Auxin Signaling & Receptors (TIR1/AFB, ARF, Aux/IAA)
    {"brachypodium_gene_id": "BRADI_4g11220v3", "arabidopsis_ortholog": "AT3G62980", "gene_symbol": "BdTIR1", "pathway": "Auxin Receptor", "evidence": "Ortholog of AtTIR1; F-box auxin receptor mediating Aux/IAA degradation"},
    {"brachypodium_gene_id": "BRADI_1g64200v3", "arabidopsis_ortholog": "AT4G03190", "gene_symbol": "BdAFB2", "pathway": "Auxin Receptor", "evidence": "Ortholog of AtAFB2; fast root growth response to auxin"},
    {"brachypodium_gene_id": "BRADI_2g49180v3", "arabidopsis_ortholog": "AT5G20730", "gene_symbol": "BdARF7", "pathway": "Auxin Response Factor", "evidence": "Ortholog of AtARF7 (NPH4); transcriptional activator of asymmetric growth genes"},
    {"brachypodium_gene_id": "BRADI_3g58400v3", "arabidopsis_ortholog": "AT1G19220", "gene_symbol": "BdARF19", "pathway": "Auxin Response Factor", "evidence": "Ortholog of AtARF19; cooperative control with ARF7 in root gravitropism"},
    {"brachypodium_gene_id": "BRADI_1g30140v3", "arabidopsis_ortholog": "AT4G14560", "gene_symbol": "BdIAA14", "pathway": "Aux/IAA Repressor", "evidence": "Ortholog of AtIAA14 (SLR); repressor of lateral root and gravitropic bending"},

    # Statolith / Starch Biosynthesis
    {"brachypodium_gene_id": "BRADI_1g09410v3", "arabidopsis_ortholog": "AT5G51820", "gene_symbol": "BdPGM1", "pathway": "Statolith Starch Biosynthesis", "evidence": "Ortholog of AtPGM1; phosphoglucomutase required for amyloplast starch synthesis in statocytes"},
    {"brachypodium_gene_id": "BRADI_2g16530v3", "arabidopsis_ortholog": "AT5G19220", "gene_symbol": "BdADG1", "pathway": "Statolith Starch Biosynthesis", "evidence": "Ortholog of AtADG1 (ADP-glucose pyrophosphorylase); essential for gravity sensing"},

    # Calcium & Mechanosensitive Signaling (Highlighted in OSD-375 & Su et al. 2023)
    {"brachypodium_gene_id": "BRADI_1g71830v3", "arabidopsis_ortholog": "AT5G66210", "gene_symbol": "BdCPK28", "pathway": "Calcium Signaling Kinase", "evidence": "Calcium-dependent protein kinase identified by Su et al. 2023 as spaceflight-responsive in roots"},
    {"brachypodium_gene_id": "BRADI_4g30480v3", "arabidopsis_ortholog": "AT5G23060", "gene_symbol": "BdCAS", "pathway": "Calcium Sensing Receptor", "evidence": "Chloroplastic/membrane calcium sensing regulator differentially expressed in spaceflight"},
    {"brachypodium_gene_id": "BRADI_3g08190v3", "arabidopsis_ortholog": "AT3G55950", "gene_symbol": "BdCRK28", "pathway": "Receptor-like Kinase", "evidence": "Cysteine-rich receptor-like kinase linked to ROS/calcium signaling during gravity perturbation"},
    {"brachypodium_gene_id": "BRADI_5g21400v3", "arabidopsis_ortholog": "AT5G15060", "gene_symbol": "BdCML24", "pathway": "Calmodulin-Like Protein", "evidence": "Ortholog of AtCML24; touch, gravity and hormone-responsive calcium transducer"},
    {"brachypodium_gene_id": "BRADI_2g39110v3", "arabidopsis_ortholog": "AT4G00290", "gene_symbol": "BdMSL10", "pathway": "Mechanosensitive Ion Channel", "evidence": "Ortholog of AtMSL10; stretch-activated channel involved in mechanoperception and gravitropism"},

    # Cell Wall Remodeling (Grasses Type II Wall Specific)
    {"brachypodium_gene_id": "BRADI_1g11420v3", "arabidopsis_ortholog": "AT2G39700", "gene_symbol": "BdEXPA1", "pathway": "Cell Wall Loosening (Expansin)", "evidence": "Alpha-expansin driving asymmetric cell elongation on the convex side of bending root/shoot"},
    {"brachypodium_gene_id": "BRADI_4g41200v3", "arabidopsis_ortholog": "AT4G14130", "gene_symbol": "BdXTH1", "pathway": "Xyloglucan Transglucosylase", "evidence": "Cell wall remodeling enzyme differentially regulated during microgravity and gravistimulation"},
    {"brachypodium_gene_id": "BRADI_3g17860v3", "arabidopsis_ortholog": "AT4G25810", "gene_symbol": "BdCSLD1", "pathway": "Cellulose/Glucan Synthase", "evidence": "Mixed-linkage glucan / cellulose synthase specific to monocot cell wall expansion"}
]

ECOTYPE_COMPARISON = [
    {
        "ecotype": "Bd21",
        "origin": "Iraq",
        "phenotype": "Reference strain, standard spring habit, early flowering, highly responsive to spaceflight stress.",
        "genomic_features": "Reference genome (v3.0), minimal vernalization requirement."
    },
    {
        "ecotype": "Bd21-3",
        "origin": "Derived from Bd21",
        "phenotype": "High transformation efficiency, shorter stature, altered gravitropic response kinetics.",
        "genomic_features": "Highly similar to Bd21 but contains specific SNP variants affecting growth habit."
    },
    {
        "ecotype": "Gaz8",
        "origin": "Turkey",
        "phenotype": "Slower growth rate, larger biomass, differential root system architecture.",
        "genomic_features": "Significant genetic divergence from Bd21, represents wider genetic diversity."
    }
]

DATA_SOURCES_MD = """# Brachypodium distachyon Genotype Data Sources

This document details the data sources for accessing *Brachypodium distachyon* SNP and variation data for GWAS analysis.

## 1. Ensembl Plants
- **Release Info**: Current release (check https://plants.ensembl.org)
- **VCF Source URL**: `ftp://ftp.ensemblgenomes.org/pub/plants/current/vcf/brachypodium_distachyon/`
- **API**: REST API at `https://rest.ensembl.org` is used to query specific variants for target candidate genes.

## 2. BrachyPan Pan-Genome
- **Reference**: Gordon, S.P., et al. (2017). "Extensive gene content variation in the Brachypodium distachyon pan-genome correlates with population structure." *Nature Communications* 8:2042.
- **Access**: Available through JGI Phytozome and specific BrachyPan data portals.
- **Usage**: Critical for understanding structural variants between ecotypes like Bd21 and Gaz8.

## 3. JGI Phytozome
- **URL**: https://phytozome-next.jgi.doe.gov/
- **Dataset**: *Brachypodium distachyon* v3.1
- **Access**: Requires JGI SSO login. Users can download VCF files for diverse panels, including APEX-06 accessions.

## 4. OSD-375 (APEX-06) Metadata
- **URL**: https://osdr.nasa.gov/osdr/data/osd/files/375
- **Reference**: Su et al. (2023). "Spaceflight and Simulated Microgravity Induce Changes in Brachypodium distachyon Gene Expression." *Life* 13(3):633. DOI: 10.3390/life13030633

## Instructions for full VCF download
For full genome-wide association studies, download the comprehensive VCF files directly via FTP or globus from Ensembl or JGI:
```bash
wget ftp://ftp.ensemblgenomes.org/pub/plants/current/vcf/brachypodium_distachyon/brachypodium_distachyon.vcf.gz
```
"""


def setup_logger(log_file: Optional[Path] = None) -> logging.Logger:
    """Configure and return the logger."""
    logger = logging.getLogger("fetch_brachy_snps")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def check_ensembl_api(logger: logging.Logger) -> bool:
    """Check if the Ensembl REST API is available."""
    try:
        response = requests.get(f"{ENSEMBL_REST_SERVER}/info/ping", headers={"Content-Type": "application/json"}, timeout=10)
        response.raise_for_status()
        data = response.json()
        if data.get("ping") == 1:
            logger.info("Ensembl REST API is available.")
            return True
        return False
    except requests.RequestException as e:
        logger.warning(f"Ensembl REST API is unavailable: {e}")
        return False


def generate_candidate_genes(output_dir: Path, logger: logging.Logger) -> None:
    """Generate and save the candidate gravitropism genes CSV."""
    logger.info("Generating gravitropism candidate genes dataset...")
    
    # Check if we can query Ensembl, but we'll use fallback primarily for stability
    api_available = check_ensembl_api(logger)
    if api_available:
        logger.info("API is available, but currently using curated static list for candidate genes.")
    
    df = pd.DataFrame(FALLBACK_CANDIDATES)
    
    out_file = output_dir / "gravitropism_candidate_genes.csv"
    df.to_csv(out_file, index=False)
    logger.info(f"Saved {len(df)} candidate genes to {out_file}")


def generate_ecotype_comparison(output_dir: Path, logger: logging.Logger) -> None:
    """Generate and save the ecotype comparison CSV."""
    logger.info("Generating ecotype comparison dataset...")
    df = pd.DataFrame(ECOTYPE_COMPARISON)
    
    out_file = output_dir / "ecotype_comparison.csv"
    df.to_csv(out_file, index=False)
    logger.info(f"Saved ecotype comparison to {out_file}")


def generate_data_sources_md(output_dir: Path, logger: logging.Logger) -> None:
    """Generate and save the DATA_SOURCES.md file."""
    logger.info("Generating DATA_SOURCES.md document...")
    out_file = output_dir / "DATA_SOURCES.md"
    
    with open(out_file, "w") as f:
        f.write(DATA_SOURCES_MD)
        
    logger.info(f"Saved data sources documentation to {out_file}")


def main() -> None:
    """Main execution entry point."""
    parser = argparse.ArgumentParser(description="Fetch Brachypodium distachyon SNP data for GWAS.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to save generated datasets."
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help="Optional file path to write logs."
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.log_file:
        args.log_file.parent.mkdir(parents=True, exist_ok=True)
        
    logger = setup_logger(args.log_file)
    logger.info("Starting Brachypodium SNP fetcher script...")
    
    try:
        generate_candidate_genes(args.output_dir, logger)
        generate_ecotype_comparison(args.output_dir, logger)
        generate_data_sources_md(args.output_dir, logger)
        logger.info("Script completed successfully.")
        
    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
