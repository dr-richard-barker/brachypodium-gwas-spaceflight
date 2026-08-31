# Brachypodium GWAS-Spaceflight Integration

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data_License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-green.svg)](environment.yml)
[![R 4.3](https://img.shields.io/badge/R-4.3-blue.svg)](environment.yml)

> **Brachypodium GWAS-Spaceflight Integration: Connecting Gravitropic Reorientation Natural Variation with ISS Microgravity Transcriptomics (OSD-375)**

This repository provides an integrated computational genomics framework linking terrestrial natural genetic variation (*Brachypodium distachyon* GWAS, pan-genomics, and gravitropic reorientation architecture) with spaceflight transcriptomic responses from the NASA Open Science Data Repository (**OSDR OSD-375 / APEX-06**).

---

## 🌟 Key Findings

1. **Ecotype-Specific Microgravity Responses (APEX-06 / OSD-375):** RNA-Seq profiling of three distinct *Brachypodium distachyon* ecotypes (**Bd21**, **Bd21-3**, and **Gaz8**) across two tissues (**root** and **shoot**) under spaceflight vs ground control conditions reveals pronounced genotype-by-environment ($G \times E$) transcriptomic divergence in response to microgravity.
2. **Gravitropism GWAS & Spaceflight Overlap:** Natural genetic variation underlying root and shoot gravitropic reorientation kinetics co-localizes significantly with differentially expressed gene (DEG) networks triggered by ISS microgravity.
3. **Monocot vs Dicot Response Divergence:** Cross-species comparative meta-analysis against *Arabidopsis thaliana* spaceflight datasets highlights conserved core stress-response pathways alongside distinct monocot-specific cell wall remodeling and auxin transport machinery.
4. **Tissue-Partitioned Transcriptional Reprogramming:** Organ-specific partitioning demonstrates distinct primary root sensing modules versus shoot photomorphogenic/circadian compensation under orbital microgravity.
5. **Predictive Machine Learning Biomarkers:** Multimodal classification and knowledge graph modeling identify key candidate quantitative trait loci (QTLs) for bioengineering cereal crop resilience in closed-loop space life support systems (BLSS).

---

## 📁 Repository Structure

```
brachypodium-gwas-spaceflight/
├── README.md                          # Project overview, quickstart & citation
├── LICENSE                            # MIT License
├── CITATION.cff                       # Citation metadata (CFF v1.2.0)
├── .zenodo.json                       # Zenodo archival metadata
├── environment.yml                    # Conda environment specification
├── Dockerfile                         # Container definition for reproducible execution
├── pyproject.toml                     # Python packaging & dependency metadata
│
├── manuscript/                        # Publication manuscript, figures & references
│   ├── manuscript.tex
│   └── references.bib
│
├── code/                              # Analytical pipeline scripts (Python & R)
│   ├── 01_download_osdr.py            # Download OSD-375 transcriptomics from NASA OSDR
│   ├── 02_download_genotypes.py       # Download Ensembl/BrachyPan/Phytozome SNP & GWAS data
│   ├── 03_differential_expression.R   # Tissue- and ecotype-level DESeq2 differential expression
│   ├── 04_gwas_spaceflight_overlap.py # Gravitropism GWAS vs spaceflight DEG overlap & enrichments
│   ├── 05_cross_species_meta.py       # Comparative analysis vs Arabidopsis spaceflight response
│   ├── 06_ml_classifier.py           # Machine learning gene & ecotype responsiveness classifiers
│   └── 07_knowledge_graph.py          # Multimodal knowledge graph & network analysis
│
├── figures/                           # High-resolution publication figures (PNG & SVG)
├── tables/                            # Processed data tables, DEGs, GWAS overlaps (CSV & JSON)
├── data/
│   ├── osdr/                          # NASA OSDR OSD-375 RNA-Seq processed counts & metadata
│   ├── genotypes/                     # VCFs, variant tables & pan-genome annotations
│   ├── gwas_phenotypes/               # Gravitropic reorientation & root architecture GWAS phenotypes
│   └── knowledge_graph/               # Cytoscape.js & GraphML network files
│
└── docs/                              # Interactive GitHub Pages web platform
    ├── assets/                        # Static media & figures
    ├── css/                           # Styling sheets
    └── js/                            # Interactive visualization scripts
```

---

## 🚀 Interactive GitHub Pages Website

The project includes an interactive web report hosted via GitHub Pages in the `docs/` directory:
- **Interactive Volcano & PCA Explorers:** Explore tissue- and ecotype-specific gene expression shifts across Bd21, Bd21-3, and Gaz8.
- **GWAS-Spaceflight Hotspot Browser:** Dynamic Manhattan & locus zoom plots intersecting gravitropism loci with spaceflight candidate genes.
- **Cross-Species Alignment:** Direct orthology comparisons between *Arabidopsis* and *Brachypodium* spaceflight transcriptomes.
- **Data Tables & Downloads:** Full access to processed DEG tables, GWAS summary statistics, and network graphs.

---

## 💻 Quickstart & Reproducibility

### Conda Environment

```bash
# Clone repository
git clone https://github.com/dr-richard-barker/brachypodium-gwas-spaceflight.git
cd brachypodium-gwas-spaceflight

# Create and activate conda environment
conda env create -f environment.yml
conda activate brachypodium-gwas-spaceflight

# Install editable package
pip install -e .
```

### Docker Container

```bash
# Build Docker image
docker build -t brachypodium-gwas-spaceflight .

# Run interactive container session
docker run -it -v $(pwd):/workspace/brachypodium-gwas-spaceflight brachypodium-gwas-spaceflight bash
```

---

## 📜 Citation & Attribution

If you use this repository, pipeline, or processed datasets in your research, please cite:

1. **Project Publication:**
   > Barker, R. (2026). *Brachypodium GWAS-Spaceflight Integration: Connecting Gravitropic Reorientation Natural Variation with ISS Microgravity Transcriptomics (OSD-375)*. GitHub: [dr-richard-barker/brachypodium-gwas-spaceflight](https://github.com/dr-richard-barker/brachypodium-gwas-spaceflight).

2. **Primary Spaceflight Dataset (APEX-06 / OSD-375):**
   > Su, S.-H., et al. (2023). *Molecular and Cellular Adaptations of Brachypodium distachyon to Microgravity in the APEX-06 Spaceflight Experiment*. **Life**, 13(3), 633. DOI: [10.3390/life13030633](https://doi.org/10.3390/life13030633).
   > NASA OSDR: [OSD-375](https://osdr.nasa.gov/osdr/data/osd/files/375) / DOI: [10.26030/2x6b-3v89](https://doi.org/10.26030/2x6b-3v89).
