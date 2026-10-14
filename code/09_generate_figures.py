#!/usr/bin/env python3
"""
09_generate_figures.py

Generate the figures for the Brachypodium spaceflight / AstroGrass resource.

Produces:
  - Fig 1: OSD-375 sample design (counts from the OSDR ISA sample table)
  - Fig 2: AstroGrass resource architecture (schematic)
  - Fig 3: spaceflight DEGs per accession, organ and direction, from Su et al. 2023
           Supplemental Tables S2-S4 (tables/su2023_deg_counts.csv)
  - Fig 4: candidate genes listed as DEGs by Su et al. 2023, with their DESeq log2 fold changes
           (tables/candidate_spaceflight_degs.csv)
  - Fig S1: APEX-06 hardware and seedlings (Su et al. 2023, CC BY 4.0) with mission facts
           taken from Su et al. 2023 and the OSD-375 ISA metadata

Figures removed on 2026-10-14 because their content was random or typed in: the gravitropism
kinetics/GWAS figure, the volcano, DEG-count and enrichment panels, the chromosome ideogram and
fold-change heatmap, the cross-species figure, the mechanistic model's spaceflight panel, and
Figs S2-S4 (synteny, promoter motifs, 46-accession correlation).

Outputs saved to figures/ (300 DPI PNG + SVG); run_pipeline.py copies them to docs/assets/.

Author: Richard Barker (ORCID: 0000-0001-5681-9857)
Affiliation: Collaborative Science Environment Inc.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

# Style configuration
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.9
plt.rcParams["font.size"] = 9.5
plt.rcParams["figure.titlesize"] = 12
plt.rcParams["svg.hashsalt"] = "brachypodium-gwas-spaceflight"  # reproducible SVG ids

SVG_META = {"Date": None}

# Up/down encoding: the validated blue/red diverging pair, with a neutral gray midpoint.
UP, DOWN, MID = "#e34948", "#2a78d6", "#f0efec"
INK, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
ACCESSIONS = ["Bd21", "Bd21-3", "Gaz8"]

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate figures.")
    parser.add_argument("--osdr-dir", type=Path, default=Path("data/osdr"), help="Path to OSDR data")
    parser.add_argument("--tables-dir", type=Path, default=Path("tables"), help="Path to tables")
    parser.add_argument("--figures-dir", type=Path, default=Path("figures"), help="Output figures directory")
    return parser.parse_args()


def save(fig, figures_dir: Path, stem: str) -> None:
    fig.savefig(figures_dir / f"{stem}.png", dpi=300)
    fig.savefig(figures_dir / f"{stem}.svg", metadata=SVG_META)
    plt.close(fig)
    logger.info(f"Saved {stem}.")


# =============================================================================
# FIGURE 1: OSD-375 sample design
# =============================================================================
def plot_figure_1(osdr_dir: Path, figures_dir: Path):
    logger.info("Generating Figure 1: OSD-375 sample design...")
    summary_file = osdr_dir / "experimental_design_summary.csv"
    if not summary_file.exists():
        raise FileNotFoundError(f"{summary_file} not found: run code/02_parse_osdr_metadata.py")

    df = pd.read_csv(summary_file)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    pivot = df.pivot_table(index=["Ecotype", "Tissue"], columns="Condition", values="Count", fill_value=0)
    pivot.index = [f"{e}\n{t.replace('Plant ', '')}" for e, t in pivot.index]
    pivot.plot(kind="bar", color=["#1565c0", "#e65100"], edgecolor="#222222", ax=ax, width=0.7)
    ax.set_ylabel("Samples (n)", fontweight="bold")
    ax.set_xlabel("Accession and organ", fontweight="bold")
    ax.set_title(f"OSD-375 sample design (n = {int(df['Count'].sum())} samples)", loc="left", fontweight="bold", fontsize=11)
    ax.tick_params(axis="x", rotation=0)
    ax.grid(True, axis="y", linestyle=":", alpha=0.6)
    ax.set_ylim(0, df["Count"].max() + 1.5)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.legend(title="Condition", frameon=True, facecolor="white", loc="upper center", ncol=2)
    fig.text(0.01, 0.01, "Source: OSDR OSD-375 ISA sample table (s_GLDS-375.txt). Each flight or ground-control "
             "APEX Growth Unit gave one root and one shoot sample.", fontsize=7.5, color="#555555")
    plt.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, figures_dir, "fig1_osd375_sample_design")


# =============================================================================
# FIGURE 2: AstroGrass architecture (schematic)
# =============================================================================
def plot_figure_2(figures_dir: Path):
    logger.info("Generating Figure 2: AstroGrass architecture...")
    fig, ax = plt.subplots(figsize=(12, 6.5), dpi=300)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.set_title("AstroGrass: grass spaceflight study catalogue and gravitropism candidate genes",
                 loc="center", fontweight="bold", fontsize=13, color="#0b1d3a")

    inputs = [
        ("NASA OSDR OSD-375", "Brachypodium, ISS (SpaceX-14)\nDEG tables: Su et al. 2023"),
        ("NASA OSDR OSD-622", "Wheat, ISS (STS-110, PESTO)\ncv. USU Apogee"),
        ("Ground studies (GEO)", "GSE97940 (2,4-D, NLS)\nGSE48040 (cold, ABR5)"),
        ("Proteomics (PRIDE)", "PXD000868\n(seedling leaf phosphoproteome)"),
    ]
    for i, (head, desc) in enumerate(inputs):
        rect = patches.FancyBboxPatch((0.5 + i*2.8, 5.2), 2.5, 1.8, boxstyle="round,pad=0.04", ec="#1565c0", fc="#e3f2fd", lw=1.5)
        ax.add_patch(rect)
        ax.text(1.75 + i*2.8, 6.4, head, ha="center", fontweight="bold", fontsize=9.5, color="#0b1d3a")
        ax.text(1.75 + i*2.8, 5.6, desc, ha="center", fontsize=8.5, color="#4a5568")
        ax.annotate("", xy=(1.75 + i*2.8, 4.3), xytext=(1.75 + i*2.8, 5.1),
                    arrowprops=dict(arrowstyle="->", lw=2, color="#0b1d3a"))

    core_box = patches.FancyBboxPatch((2.0, 2.6), 8.0, 1.6, boxstyle="round,pad=0.05", ec="#2e7d32", fc="#e8f5e9", lw=2)
    ax.add_patch(core_box)
    ax.text(6.0, 3.7, "AstroGrass catalogue and candidate-gene table", ha="center", fontweight="bold", fontsize=11, color="#1b5e20")
    ax.text(6.0, 3.2, "Study records with checked DOIs • 29 gravitropism candidates: UniProt-checked Arabidopsis anchors,",
            ha="center", fontsize=9, color="#2e7d32")
    ax.text(6.0, 2.85, "Ensembl Compara orthologs, and spaceflight DEG calls from Su et al. 2023 (Tables S2-S4) • CSV and JSON",
            ha="center", fontsize=9, color="#2e7d32")

    ax.annotate("", xy=(3.5, 1.6), xytext=(4.0, 2.5), arrowprops=dict(arrowstyle="->", lw=2, color="#2e7d32"))
    ax.annotate("", xy=(8.5, 1.6), xytext=(8.0, 2.5), arrowprops=dict(arrowstyle="->", lw=2, color="#2e7d32"))

    deliv_1 = patches.FancyBboxPatch((1.0, 0.4), 4.5, 1.1, boxstyle="round,pad=0.04", ec="#e65100", fc="#fff3e0", lw=1.5)
    ax.add_patch(deliv_1)
    ax.text(3.25, 1.0, "Web explorer", ha="center", fontweight="bold", fontsize=10, color="#e65100")
    ax.text(3.25, 0.65, "docs/astrograss.html (search and filter)", ha="center", fontsize=8.5, color="#555")

    deliv_2 = patches.FancyBboxPatch((6.5, 0.4), 4.5, 1.1, boxstyle="round,pad=0.04", ec="#4a148c", fc="#f3e5f5", lw=1.5)
    ax.add_patch(deliv_2)
    ax.text(8.75, 1.0, "Open repository", ha="center", fontweight="bold", fontsize=10, color="#4a148c")
    ax.text(8.75, 0.65, "GitHub: code, tables and figures", ha="center", fontsize=8.5, color="#555")

    plt.tight_layout()
    save(fig, figures_dir, "fig2_astrograss_architecture")


# =============================================================================
# FIGURE 3: Su et al. 2023 DEG counts
# =============================================================================
def plot_figure_3(tables_dir: Path, figures_dir: Path):
    logger.info("Generating Figure 3: Su et al. 2023 DEG counts...")
    counts = pd.read_csv(tables_dir / "su2023_deg_counts.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), dpi=300)
    x = np.arange(len(ACCESSIONS))
    for ax, organ in zip(axes, ["root", "shoot"]):
        sub = counts[counts["organ"] == organ].set_index("accession").loc[ACCESSIONS]
        ax.bar(x, sub["up"], width=0.55, color=UP, label="Up in spaceflight")
        ax.bar(x, -sub["down"], width=0.55, color=DOWN, label="Down in spaceflight")
        top = max(sub["up"].max(), sub["down"].max())
        for xi, (u, d) in enumerate(zip(sub["up"], sub["down"])):
            ax.text(xi, u + top * 0.02, f"{u:,}", ha="center", va="bottom", fontsize=9, color=INK)
            ax.text(xi, -d - top * 0.02, f"{d:,}", ha="center", va="top", fontsize=9, color=INK)
        ax.axhline(0, color=INK_MUTED, lw=0.8)
        ax.set_ylim(-top * 1.2, top * 1.2)
        ax.set_xticks(x)
        ax.set_xticklabels(ACCESSIONS)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(int(v)):,}"))
        ax.set_ylabel("Genes (up above, down below zero)", color=INK_MUTED)
        ax.set_title(f"{organ.capitalize()}s", loc="left", fontweight="bold", fontsize=11)
        ax.grid(True, axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[1].legend(frameon=False, loc="upper left")
    fig.suptitle("Spaceflight DEGs reported by Su et al. (2023), ISS vs ground control",
                 x=0.01, ha="left", fontweight="bold", fontsize=12)
    fig.text(0.01, 0.01, "Source: Su et al. 2023, Life 13:626, Supplemental Tables S2-S4 (CC BY 4.0). Genes listed as DEGs "
             "(DESeq padj and edgeR q < 0.05); direction from the DESeq log2 fold change. Panels have separate y scales.\n"
             "Bd21 was sequenced at a different facility and read length from Bd21-3 and Gaz8, so differences between "
             "accessions are confounded with sequencing batch.", fontsize=7.5, color=INK_MUTED)
    plt.tight_layout(rect=(0, 0.07, 1, 0.95))
    save(fig, figures_dir, "fig3_su2023_deg_counts")


# =============================================================================
# FIGURE 4: candidate genes in Su et al. 2023's DEG lists
# =============================================================================
def plot_figure_4(tables_dir: Path, figures_dir: Path):
    logger.info("Generating Figure 4: candidate genes in Su et al. 2023's DEG lists...")
    hits = pd.read_csv(tables_dir / "candidate_spaceflight_degs.csv")
    columns = [(a, o) for a in ACCESSIONS for o in ("root", "shoot")]
    genes = list(dict.fromkeys(hits.sort_values(["arabidopsis_gene", "brachypodium_gene_id"])["brachypodium_gene_id"]))

    labels = []
    for gid in genes:
        g = hits[hits["brachypodium_gene_id"] == gid]
        names = []
        for anchor, m in g[["arabidopsis_gene", "match"]].drop_duplicates().itertuples(index=False):
            names.append(anchor + ("" if m.startswith("Compara") else "*"))
        labels.append(f"{gid}  ({', '.join(dict.fromkeys(names))})")

    values = np.full((len(genes), len(columns)), np.nan)
    for i, gid in enumerate(genes):
        g = hits[hits["brachypodium_gene_id"] == gid].drop_duplicates(["accession", "organ"])
        for _, r in g.iterrows():
            values[i, columns.index((r["accession"], r["organ"]))] = r["deseq_log2fc"]

    cmap = LinearSegmentedColormap.from_list("updown", [DOWN, MID, UP])
    lim = float(np.ceil(np.nanmax(np.abs(values))))
    fig, ax = plt.subplots(figsize=(10, 0.45 * len(genes) + 2.2), dpi=300)
    ax.set_facecolor("#fcfcfb")
    im = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=-lim, vmax=lim, aspect="auto")
    for i in range(len(genes)):
        for j in range(len(columns)):
            v = values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8.5,
                        color="white" if abs(v) > lim * 0.55 else INK)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([f"{a}\n{o}" for a, o in columns])
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xticks(np.arange(-0.5, len(columns)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(genes)), minor=True)
    ax.grid(which="minor", color="white", lw=2)
    ax.tick_params(which="minor", length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("DESeq log2 fold change (ISS / ground)", color=INK_MUTED)
    ax.set_title("Candidate genes listed as spaceflight DEGs by Su et al. (2023)", loc="left", fontweight="bold", fontsize=11)
    fig.text(0.01, 0.01, "Values from Su et al. 2023, Life 13:626, Supplemental Tables S2-S4 (CC BY 4.0). Blank: not listed as a DEG "
             "(not a measured zero). Genes are Ensembl Compara orthologs of the\ncandidate's Arabidopsis anchor; * marks a match "
             "through Su et al.'s best-hit Arabidopsis ortholog only. Accession differences are confounded with sequencing batch.",
             fontsize=7.5, color=INK_MUTED)
    plt.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, figures_dir, "fig4_candidate_spaceflight_degs")


# =============================================================================
# SUPPLEMENTARY FIGURE S1: APEX-06 hardware, seedlings and mission facts
# =============================================================================
def plot_figure_s1(figures_dir: Path):
    logger.info("Generating Supplementary Figure S1: APEX-06 hardware and mission...")
    fig = plt.figure(figsize=(15, 9.5), dpi=300)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1], height_ratios=[1.6, 1])

    for pos, img_name, title in [
        (gs[0, 0], "apex06_hardware_growth_unit.png", "A. APEX Growth Unit (Su et al. 2023, Fig. 1)"),
        (gs[0, 1], "apex06_brachypodium_seedlings_iss.png", "B. ISS and ground-control seedlings (Su et al. 2023, Fig. 2)"),
    ]:
        ax = fig.add_subplot(pos)
        ax.axis("off")
        img_path = figures_dir / img_name
        if not img_path.exists():
            raise FileNotFoundError(img_path)
        ax.imshow(mpimg.imread(str(img_path)))
        ax.set_title(title, loc="left", fontweight="bold", fontsize=11)

    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.axis("off")
    ax_c.set_title("C. APEX-06 timeline (Su et al. 2023)", loc="left", fontweight="bold", fontsize=11)
    ax_c.add_patch(patches.FancyBboxPatch((0.02, 0.05), 0.96, 0.88, boxstyle="round,pad=0.03",
                                          ec="#1565c0", fc="#e3f2fd", lw=1.5))
    timeline_text = (
        "• Launch: SpaceX-14, 2 April 2018; docked with the ISS about 2 days later\n"
        "• Activation: 12 April 2018; 90 mL half-strength MS medium per APEX Growth Unit\n"
        "• 12 flight units (4 per accession), 24 seeds each, in VEGGIE\n"
        "• Germination: 24 h under red LEDs, then 4 days under red, blue and green LEDs\n"
        "• Fixation: RNAlater in Kennedy Fixation Tubes, 24 h at room temperature\n"
        "• Storage: MELFI (-80 °C) for 20 days; returned to Earth on 5 May 2018\n"
        "• Sequencing: Illumina NovaSeq 6000; Bd21 2×150 bp (UW–Madison),\n"
        "   Bd21-3 and Gaz8 2×100 bp (NASA GeneLab)"
    )
    ax_c.text(0.06, 0.5, timeline_text, va="center", fontsize=9.5, color="#0b1d3a", linespacing=1.5)

    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis("off")
    ax_d.set_title("D. Growth conditions (OSD-375 ISA metadata; Su et al. 2023)", loc="left", fontweight="bold", fontsize=11)
    ax_d.add_patch(patches.FancyBboxPatch((0.02, 0.05), 0.96, 0.88, boxstyle="round,pad=0.03",
                                          ec="#2e7d32", fc="#e8f5e9", lw=1.5))
    env_text = (
        "• Growth temperature: 22 °C\n"
        "• Light: continuous; red during germination, then red/blue/green\n"
        "• Light intensity: 40–55 µmol m⁻² s⁻¹\n"
        "• Growth time: 4 days after germination\n"
        "• Absorbed radiation dose (flight): 1.397 mGy (0.349 mGy/day)\n"
        "• Ground control: VEGGIE in the ISS Environmental Simulator at\n"
        "   Kennedy Space Center, run with a 2-day delay\n"
        "• Samples: 48 (3 accessions × 2 organs × flight/ground × 4 units)"
    )
    ax_d.text(0.06, 0.5, env_text, va="center", fontsize=9.5, color="#1b5e20", linespacing=1.5)

    fig.text(0.01, 0.005, "Panels A and B reproduced from Su, Levine & Masson (2023) Life 13:626, "
             "doi:10.3390/life13030626, under CC BY 4.0.", fontsize=8, color="#555555")
    plt.tight_layout(rect=(0, 0.02, 1, 1))
    save(fig, figures_dir, "figS1_mission_experimental_design")


def main():
    args = parse_args()
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    plot_figure_1(args.osdr_dir, args.figures_dir)
    plot_figure_2(args.figures_dir)
    plot_figure_3(args.tables_dir, args.figures_dir)
    plot_figure_4(args.tables_dir, args.figures_dir)
    plot_figure_s1(args.figures_dir)
    logger.info("✓ Figures 1-4 and S1 generated in figures/")


if __name__ == "__main__":
    main()
