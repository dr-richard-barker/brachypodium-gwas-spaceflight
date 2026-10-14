#!/usr/bin/env python3
"""
04_gravitropism_kinetics.py

Summarise root gravitropic curvature measurements per Brachypodium accession and test for
differences between accessions (one-way ANOVA at each stimulation time).

No measurements are in the repository yet. The gravitropism_data.csv committed before
2026-10-14 was produced by this script's former --demo mode (seeded random numbers), and the
"GWAS" step only copied the candidate-gene list; both were removed. This script is not part of
run_pipeline.py until measured data, with a documented source, are added.

Input columns (CSV):
    accession, replicate, stimulation_min, curvature_angle, root_length_mm
    (optional: root_angle_pre, root_angle_post, measurement_tool, notes)

Usage:
    python code/04_gravitropism_kinetics.py --phenotype-file <measurements.csv>

Outputs (in --output-dir, default tables/):
    gravitropism_kinetics.csv, gravitropism_summary.csv, gravitropism_anova.csv

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = ["accession", "stimulation_min", "curvature_angle", "root_length_mm"]

# Accessions flown in NASA OSDR OSD-375, flagged in the summary table
OSD375_ACCESSIONS = ["Bd21", "Bd21-3", "Gaz8"]


def load_phenotype_data(filepath: Path) -> pd.DataFrame:
    """Load and validate gravitropic curvature measurements."""
    if not filepath.exists():
        logger.error(f"Phenotype file not found: {filepath}")
        logger.info("Expected columns: " + ", ".join(REQUIRED_COLUMNS))
        sys.exit(1)

    df = pd.read_csv(filepath)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        logger.error(f"Missing required columns: {missing}")
        sys.exit(1)

    logger.info(f"Loaded {len(df)} observations from {df['accession'].nunique()} accessions")
    return df


def compute_reorientation_kinetics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mean, SD and SE of curvature per accession and stimulation time, plus a per-accession summary."""
    kinetics = (
        df.groupby(["accession", "stimulation_min"])
        .agg(
            mean_curvature=("curvature_angle", "mean"),
            sd_curvature=("curvature_angle", "std"),
            se_curvature=("curvature_angle", lambda x: x.std() / np.sqrt(len(x))),
            n_reps=("curvature_angle", "count"),
            mean_root_length=("root_length_mm", "mean"),
        )
        .reset_index()
    )

    last_time = df["stimulation_min"].max()
    at_last = df[df["stimulation_min"] == last_time].groupby("accession")["curvature_angle"].mean()
    summary = (
        df.groupby("accession")
        .agg(
            overall_mean_curvature=("curvature_angle", "mean"),
            overall_sd=("curvature_angle", "std"),
            n_total=("curvature_angle", "count"),
        )
        .reset_index()
    )
    summary[f"mean_curvature_{last_time}min"] = summary["accession"].map(at_last)
    summary["in_osd375"] = summary["accession"].isin(OSD375_ACCESSIONS)
    return kinetics, summary.sort_values("overall_mean_curvature", ascending=False)


def run_anova(df: pd.DataFrame) -> pd.DataFrame:
    """One-way ANOVA of curvature across accessions at each stimulation time."""
    results = {}
    for stim_time in sorted(df["stimulation_min"].unique()):
        subset = df[df["stimulation_min"] == stim_time]
        groups = [g["curvature_angle"].values for _, g in subset.groupby("accession")]
        if len(groups) >= 2:
            f_stat, p_val = stats.f_oneway(*groups)
            results[stim_time] = {"f_statistic": f_stat, "p_value": p_val, "n_accessions": len(groups)}
            logger.info(f"  {stim_time} min: F={f_stat:.3f}, p={p_val:.2e}, n_groups={len(groups)}")
    out = pd.DataFrame.from_dict(results, orient="index")
    out.index.name = "stimulation_min"
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarise gravitropic curvature measurements.")
    parser.add_argument("--phenotype-file", type=Path, required=True,
                        help="CSV of measured curvature angles (see module docstring)")
    parser.add_argument("--output-dir", type=Path, default=Path("tables"), help="Output directory")
    args = parser.parse_args()

    df = load_phenotype_data(args.phenotype_file)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    kinetics, summary = compute_reorientation_kinetics(df)
    kinetics.to_csv(args.output_dir / "gravitropism_kinetics.csv", index=False)
    summary.to_csv(args.output_dir / "gravitropism_summary.csv", index=False)
    run_anova(df).to_csv(args.output_dir / "gravitropism_anova.csv")
    logger.info(f"Saved kinetics, summary and ANOVA tables to {args.output_dir}/")


if __name__ == "__main__":
    main()
