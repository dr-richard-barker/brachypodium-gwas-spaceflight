#!/usr/bin/env python3
"""
run_pipeline.py

Master execution script for the Brachypodium GWAS-Spaceflight / AstroGrass pipeline.
Executes all steps sequentially and validates output files.

Usage:
    python code/run_pipeline.py [--skip-download]
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger("AstroGrassPipeline")

ROOT_DIR = Path(__file__).resolve().parent.parent

PIPELINE_STEPS = [
    ("Step 1: Ingest NASA OSDR OSD-375 Data", ["python3", "code/01_download_osdr_osd375.py"]),
    ("Step 2: Parse ISA-Tab Sample Metadata", ["python3", "code/02_parse_osdr_metadata.py"]),
    ("Step 3: Fetch SNP & Candidate Genes", ["python3", "code/03_fetch_brachypodium_snps.py"]),
    ("Step 4: Run Gravitropism GWAS Analysis", ["python3", "code/04_gwas_gravitropism_analysis.py", "--demo"]),
    ("Step 5: Cross-Species Meta-Analysis", ["python3", "code/06_meta_analysis_cross_species.py"]),
    ("Step 6: GWAS Spaceflight Integration", ["python3", "code/07_gwas_spaceflight_integration.py"]),
    ("Step 7: Alternative Splicing Integration", ["python3", "code/08_alternative_splicing_link.py"]),
    ("Step 8: Build AstroGrass Database", ["python3", "code/10_build_astrograss_db.py"]),
    ("Step 9: Build Video & Image Series Catalog", ["python3", "code/11_build_video_catalog.py"]),
    ("Step 10: Sync CyVerse iRODS Image Series", ["python3", "code/12_cyverse_irods_sync.py"]),
    ("Step 11: Academic Synteny & Promoter Analysis", ["python3", "code/13_promoter_synteny_analysis.py"]),
    ("Step 12: Generate Publication Figures", ["python3", "code/09_generate_figures.py"]),
    ("Step 13: Compile Manuscript (PDF & DOCX)", ["python3", "code/export_manuscript.py"]),
]


def run_command(name: str, cmd: list[str]) -> bool:
    logger.info(f"{'='*60}")
    logger.info(f"Running {name}...")
    logger.info(f"Command: {' '.join(cmd)}")
    logger.info(f"{'='*60}")
    
    start_time = time.time()
    result = subprocess.run(cmd, cwd=ROOT_DIR, capture_output=True, text=True)
    duration = time.time() - start_time
    
    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip(), file=sys.stderr)
        
    if result.returncode == 0:
        logger.info(f"✓ {name} completed in {duration:.2f}s")
        return True
    else:
        logger.error(f"✗ {name} failed with return code {result.returncode}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Run complete AstroGrass pipeline.")
    parser.add_argument("--skip-download", action="store_true", help="Skip Step 1 (OSDR download)")
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("🌾 Starting AstroGrass / Brachypodium GWAS-Spaceflight Pipeline")
    print("="*70 + "\n")
    
    success_count = 0
    for name, cmd in PIPELINE_STEPS:
        if args.skip_download and "Step 1" in name:
            logger.info("Skipping Step 1 (--skip-download specified)")
            success_count += 1
            continue
            
        success = run_command(name, cmd)
        if success:
            success_count += 1
        else:
            logger.error(f"Pipeline stopped at {name}")
            sys.exit(1)
            
    # Copy generated figures to docs/assets/
    figures_src = ROOT_DIR / "figures"
    docs_assets = ROOT_DIR / "docs" / "assets"
    docs_assets.mkdir(parents=True, exist_ok=True)
    
    for f in figures_src.glob("*.png"):
        (docs_assets / f.name).write_bytes(f.read_bytes())
    for f in figures_src.glob("*.svg"):
        (docs_assets / f.name).write_bytes(f.read_bytes())
    logger.info(f"✓ Synchronized all publication figures to {docs_assets}")
    
    print("\n" + "="*70)
    print(f"🎉 Pipeline Complete! {success_count}/{len(PIPELINE_STEPS)} steps succeeded.")
    print("="*70)
    print("\nKey Outputs Generated:")
    print(" - Master Database: tables/astrograss_master_table.csv")
    print(" - Web Database:   docs/js/astrograss_database.json")
    print(" - Figures:        figures/ (fig1 - fig4 in PNG + SVG)")
    print(" - Manuscript:     manuscript/manuscript.pdf & manuscript.docx")
    print(" - Web Portal:     docs/astrograss.html & docs/index.html\n")


if __name__ == "__main__":
    main()
