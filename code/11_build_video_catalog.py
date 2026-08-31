#!/usr/bin/env python3
"""
11_build_video_catalog.py

Compile and maintain the AstroGrass Video & Time-Lapse Knowledgebase.
Combines:
  1. Curated YouTube scientific videos (Brachypodium 360°, Clinostats, Gravitropism assays)
  2. CyVerse iRODS Brachy_* time-lapse series (Per1, TR7a, Gravatron, FlashLapse)
  3. Interactive canvas kinematic simulators

Outputs:
  - tables/video_catalog.csv
  - docs/tables/video_catalog.csv
  - data/videos/video_catalog.json
  - docs/js/video_database.json

Author: Richard Barker (ORCID: 0000-0002-4525-3341)
Affiliation: Phylo
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

# Master Video Catalog Data
MASTER_VIDEO_CATALOG: List[Dict[str, Any]] = [
    {
        "id": "vid-yt-01",
        "type": "youtube",
        "youtube_id": "DOhR85UD-HE",
        "title": "Brachypodium 360° View",
        "author": "Richard Poiré",
        "author_url": "https://www.youtube.com/@RichardPoir%C3%A9",
        "category": "Brachypodium Morphology",
        "organism": "Brachypodium distachyon",
        "description": "High-resolution 360-degree rotational optical scan of Brachypodium distachyon mature canopy architecture, tillering morphology, and leaf angle distribution.",
        "thumbnail_url": "https://i.ytimg.com/vi/DOhR85UD-HE/hqdefault.jpg",
        "url": "https://www.youtube.com/watch?v=DOhR85UD-HE",
        "tags": ["Brachypodium", "Morphology", "360 View", "Canopy", "Phenomics"]
    },
    {
        "id": "vid-yt-02",
        "type": "youtube",
        "youtube_id": "7ilLAd4-phY",
        "title": "Clinostat Simulated Microgravity Assay",
        "author": "sathish sre",
        "author_url": "https://www.youtube.com/@sathishsre4914",
        "category": "Clinostat / Microgravity",
        "organism": "Plant Biology / Ground Analogue",
        "description": "2D Clinostat experimental hardware operation. Omnilateral continuous rotation at constant RPM disperses unidirectional gravity vectors to simulate spaceflight microgravity in laboratory settings.",
        "thumbnail_url": "https://i.ytimg.com/vi/7ilLAd4-phY/hqdefault.jpg",
        "url": "https://www.youtube.com/watch?v=7ilLAd4-phY",
        "tags": ["Clinostat", "Microgravity", "Hardware", "Ground Analogue", "Rotation"]
    },
    {
        "id": "vid-yt-03",
        "type": "youtube",
        "youtube_id": "YIqCwYKDBqA",
        "title": "WT Arabidopsis Root Gravitropism Time-lapse (15 min interval, 10 fps)",
        "author": "Dr. Richard Barker",
        "author_url": "https://www.youtube.com/@DrRichardBarker",
        "category": "Gravitropism Assays",
        "organism": "Arabidopsis thaliana / Comparative Model",
        "description": "Optical time-lapse capture of root tip gravitropic reorientation following 90-degree reorientation. Frames acquired every 15 minutes at 10 fps, illustrating statolith displacement and asymmetric auxin kinetics.",
        "thumbnail_url": "https://i.ytimg.com/vi/YIqCwYKDBqA/hqdefault.jpg",
        "url": "https://www.youtube.com/watch?v=YIqCwYKDBqA",
        "tags": ["Gravitropism", "Time-lapse", "Auxin", "Root Curvature", "Astrobotany"]
    },
    # CyVerse Brachypodium Time-Lapse Series
    {
        "id": "cyverse-brachy-per1-1",
        "type": "cyverse",
        "series_name": "Brachy_Enhanced_Per1_DB_1_1",
        "title": "Brachypodium Per1 Accession - Series 1.1",
        "author": "Dr. Richard Barker (CyVerse / iRODS)",
        "category": "CyVerse Time-Lapses",
        "organism": "Brachypodium distachyon",
        "accession": "Per1",
        "frames": 84,
        "cadence_seconds": 15.0,
        "irods_path": "/iplant/home/dr_richard_barker/timelapse_extract/Brachy_Enhanced_Per1_DB_1_1",
        "video_url": "https://raw.githubusercontent.com/dr-richard-barker/timelapse-image-series/main/Brachy_Enhanced_Per1_DB_1_1/timelapse.mp4",
        "description": "High-contrast enhanced time-lapse series of Brachypodium distachyon Per1 accession root gravitropic kinematics extracted from CyVerse Data Store.",
        "thumbnail_url": "assets/fig1_gravitropic_kinetics.png",
        "tags": ["CyVerse", "Per1", "Brachypodium", "Time-lapse", "iRODS"]
    },
    {
        "id": "cyverse-brachy-per1-2",
        "type": "cyverse",
        "series_name": "Brachy_Enhanced_Per1_DB_2_3",
        "title": "Brachypodium Per1 Accession - Series 2.3",
        "author": "Dr. Richard Barker (CyVerse / iRODS)",
        "category": "CyVerse Time-Lapses",
        "organism": "Brachypodium distachyon",
        "accession": "Per1",
        "frames": 84,
        "cadence_seconds": 15.0,
        "irods_path": "/iplant/home/dr_richard_barker/timelapse_extract/Brachy_Enhanced_Per1_DB_2_3",
        "video_url": "https://raw.githubusercontent.com/dr-richard-barker/timelapse-image-series/main/Brachy_Enhanced_Per1_DB_2_3/timelapse.mp4",
        "description": "Replicate time-lapse series of Brachypodium distachyon Per1 accession root growth under controlled gravistimulation.",
        "thumbnail_url": "assets/fig1_gravitropic_kinetics.png",
        "tags": ["CyVerse", "Per1", "Brachypodium", "Time-lapse", "iRODS"]
    },
    {
        "id": "cyverse-brachy-per1-3",
        "type": "cyverse",
        "series_name": "Brachy_Enhanced_Per1_DB_3_2",
        "title": "Brachypodium Per1 Accession - Series 3.2",
        "author": "Dr. Richard Barker (CyVerse / iRODS)",
        "category": "CyVerse Time-Lapses",
        "organism": "Brachypodium distachyon",
        "accession": "Per1",
        "frames": 83,
        "cadence_seconds": 15.0,
        "irods_path": "/iplant/home/dr_richard_barker/timelapse_extract/Brachy_Enhanced_Per1_DB_3_2",
        "video_url": "https://raw.githubusercontent.com/dr-richard-barker/timelapse-image-series/main/Brachy_Enhanced_Per1_DB_3_2/timelapse.mp4",
        "description": "Gravitropic reorientation series for Brachypodium distachyon accession Per1.",
        "thumbnail_url": "assets/fig1_gravitropic_kinetics.png",
        "tags": ["CyVerse", "Per1", "Brachypodium", "Time-lapse", "iRODS"]
    },
    {
        "id": "cyverse-brachy-tr7a-1",
        "type": "cyverse",
        "series_name": "Brachy_Enhanced_TR7a_DB_1_3",
        "title": "Brachypodium TR7a Accession - Series 1.3",
        "author": "Dr. Richard Barker (CyVerse / iRODS)",
        "category": "CyVerse Time-Lapses",
        "organism": "Brachypodium distachyon",
        "accession": "BdTR7a",
        "frames": 83,
        "cadence_seconds": 15.0,
        "irods_path": "/iplant/home/dr_richard_barker/timelapse_extract/Brachy_Enhanced_TR7a_DB_1_3",
        "video_url": "https://raw.githubusercontent.com/dr-richard-barker/timelapse-image-series/main/Brachy_Enhanced_TR7a_DB_1_3/timelapse.mp4",
        "description": "High-contrast time-lapse capture of Brachypodium distachyon Turkish accession BdTR7a during gravitropic curvature.",
        "thumbnail_url": "assets/fig1_gravitropic_kinetics.png",
        "tags": ["CyVerse", "BdTR7a", "Brachypodium", "Turkish Diversity", "iRODS"]
    },
    {
        "id": "cyverse-brachy-tr7a-2",
        "type": "cyverse",
        "series_name": "Brachy_Enhanced_TR7a_DB_2_2",
        "title": "Brachypodium TR7a Accession - Series 2.2",
        "author": "Dr. Richard Barker (CyVerse / iRODS)",
        "category": "CyVerse Time-Lapses",
        "organism": "Brachypodium distachyon",
        "accession": "BdTR7a",
        "frames": 84,
        "cadence_seconds": 15.0,
        "irods_path": "/iplant/home/dr_richard_barker/timelapse_extract/Brachy_Enhanced_TR7a_DB_2_2",
        "video_url": "https://raw.githubusercontent.com/dr-richard-barker/timelapse-image-series/main/Brachy_Enhanced_TR7a_DB_2_2/timelapse.mp4",
        "description": "Replicate kinematic series of Brachypodium distachyon accession BdTR7a.",
        "thumbnail_url": "assets/fig1_gravitropic_kinetics.png",
        "tags": ["CyVerse", "BdTR7a", "Brachypodium", "Turkish Diversity", "iRODS"]
    },
    {
        "id": "cyverse-gravatron-demo",
        "type": "cyverse",
        "series_name": "Gravatron_FlashLapse_DEMO",
        "title": "Gravatron FlashLapse Microgravity Demonstration",
        "author": "Dr. Richard Barker (CyVerse / iRODS)",
        "category": "Clinostat / Microgravity",
        "organism": "Plant Biology / Multi-species",
        "frames": 200,
        "cadence_seconds": 30.0,
        "irods_path": "/iplant/home/dr_richard_barker/timelapse_extract/Gravatron_FlashLapse_DEMO",
        "video_url": "https://raw.githubusercontent.com/dr-richard-barker/timelapse-image-series/main/Gravatron_FlashLapse_DEMO/timelapse.mp4",
        "description": "Demonstration series of seedling reorientation inside the Gravatron FlashLapse simulated microgravity rotational apparatus.",
        "thumbnail_url": "assets/fig1_gravitropic_kinetics.png",
        "tags": ["Gravatron", "FlashLapse", "Microgravity", "Simulation"]
    },
    {
        "id": "cyverse-flashlapse-straight",
        "type": "cyverse",
        "series_name": "FlashLapse_Straight_growth",
        "title": "FlashLapse Straight Growth Baseline Control",
        "author": "Dr. Richard Barker (CyVerse / iRODS)",
        "category": "Gravitropism Assays",
        "organism": "Pisum sativum / Model Seedling",
        "frames": 612,
        "cadence_seconds": 15.0,
        "irods_path": "/iplant/home/dr_richard_barker/timelapse_extract/FlashLapse_Straight_growth",
        "video_url": "https://raw.githubusercontent.com/dr-richard-barker/timelapse-image-series/main/FlashLapse_Straight_growth/timelapse.mp4",
        "description": "High-frequency 612-frame straight vertical growth baseline capture for kinematic calibration.",
        "thumbnail_url": "assets/fig1_gravitropic_kinetics.png",
        "tags": ["FlashLapse", "Baseline", "Control", "Calibration"]
    },
    # Interactive Simulators
    {
        "id": "sim-90deg-reorientation",
        "type": "simulator",
        "title": "90° Gravistimulation Kinetic Reorientation Simulator",
        "author": "AstroGrass Modeling Engine",
        "category": "Interactive Simulators",
        "organism": "Brachypodium distachyon",
        "description": "Real-time HTML5 Canvas biophysical simulation of monocot root columella bending from 0° to 42° following a 90° gravistimulation vector shift.",
        "thumbnail_url": "assets/fig5_mechanistic_model.png",
        "tags": ["Simulator", "Kinematics", "1g", "Interactive", "Canvas"]
    },
    {
        "id": "sim-2d-clinostat",
        "type": "simulator",
        "title": "2D Clinostat Omnilateral Rotation Simulator (1.0 RPM)",
        "author": "AstroGrass Modeling Engine",
        "category": "Interactive Simulators",
        "organism": "Plant Biology / Ground Analogue",
        "description": "Interactive rotational simulation modeling the dispersion of gravitational sedimentation vectors during 1.0 RPM clinorotation.",
        "thumbnail_url": "assets/fig5_mechanistic_model.png",
        "tags": ["Simulator", "Clinostat", "RPM", "Microgravity", "Interactive"]
    }
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build AstroGrass video catalog.")
    parser.add_argument("--tables-dir", type=Path, default=Path("tables"), help="Tables directory")
    parser.add_argument("--docs-dir", type=Path, default=Path("docs"), help="Docs directory")
    parser.add_argument("--data-dir", type=Path, default=Path("data/videos"), help="Data videos directory")
    return parser.parse_args()


def build_catalog(tables_dir: Path, docs_dir: Path, data_dir: Path):
    logger.info("Building AstroGrass Video & Image Series Catalog...")
    tables_dir.mkdir(parents=True, exist_ok=True)
    docs_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    # 1. Export JSON catalog for Data Store & Web App
    json_path = data_dir / "video_catalog.json"
    web_json_path = docs_dir / "js" / "video_database.json"
    web_json_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "metadata": {
            "name": "AstroGrass Video Knowledgebase",
            "version": "1.1.0",
            "curator": "Richard Barker (Phylo)",
            "total_items": len(MASTER_VIDEO_CATALOG),
            "youtube_videos": sum(1 for v in MASTER_VIDEO_CATALOG if v["type"] == "youtube"),
            "cyverse_series": sum(1 for v in MASTER_VIDEO_CATALOG if v["type"] == "cyverse"),
            "simulators": sum(1 for v in MASTER_VIDEO_CATALOG if v["type"] == "simulator"),
            "updated_at": "2026-08-31"
        },
        "videos": MASTER_VIDEO_CATALOG
    }

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    with open(web_json_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    logger.info(f"Saved Video Catalog JSON to {json_path} and {web_json_path}")

    # 2. Export CSV Table
    df = pd.DataFrame(MASTER_VIDEO_CATALOG)
    # Convert list tags to string for CSV
    if "tags" in df.columns:
        df["tags"] = df["tags"].apply(lambda t: "; ".join(t) if isinstance(t, list) else str(t))

    csv_path = tables_dir / "video_catalog.csv"
    df.to_csv(csv_path, index=False)
    logger.info(f"Saved Video Catalog CSV ({len(df)} records) to {csv_path}")

    # Sync to docs/tables/
    docs_tables = docs_dir / "tables"
    docs_tables.mkdir(parents=True, exist_ok=True)
    (docs_tables / "video_catalog.csv").write_bytes(csv_path.read_bytes())
    logger.info(f"✓ Synchronized video catalog to {docs_tables / 'video_catalog.csv'}")


def main():
    args = parse_args()
    build_catalog(args.tables_dir, args.docs_dir, args.data_dir)
    print("\n✓ AstroGrass video catalog compilation complete.")


if __name__ == "__main__":
    main()
