# Provenance of results files

Which code produced each committed table, data file and figure, and from what. Written 14 Oct 2026. Scripts are in
`code/` and run from the repository root; `code/run_pipeline.py --skip-download` reruns every offline step, and CI
fails if a rebuilt table or JSON differs from the committed one. Scripts marked *network* are run by hand and their
outputs committed.

## Inputs (not generated)

| File | What it is |
|---|---|
| `data/osdr/Study_Metadata_Files/` | OSD-375 ISA-Tab metadata, file listing and raw MultiQC report from NASA OSDR (`01_download_osdr_osd375.py`, network) |
| `data/genotypes/candidate_genes_curated.csv` | Hand-curated candidate list; anchor corrections of 2026-10-14 chosen by hand from UniProt gene-name searches |
| `data/su2023/su2023_supplementary_tables.xlsx` | Su et al. (2023) Supplementary Tables, CC BY 4.0; source and SHA-256 in `data/su2023/SOURCE.md` |
| `figures/apex06_*.png`, `figures/life-13-00626-g00*.png` | Su et al. (2023) Figures 1 and 2, CC BY 4.0 (byte-identical pairs) |
| `data/phenotypes/` | "Barker_2016" orientation and displacement plots (PNG only; source and numbers TODO) |

## Tables and data

| File | Produced by | From |
|---|---|---|
| `data/osdr/sample_metadata.csv`, `experimental_design_summary.csv`, `sequencing_parameters.csv` | `02_parse_osdr_metadata.py` | OSD-375 ISA sample and assay tables |
| `data/genotypes/arabidopsis_anchor_check.csv` | `14_check_arabidopsis_anchors.py` (network) | UniProt REST API |
| `data/genotypes/compara_orthologs.csv` | `15_map_orthologs.py` (network) | Ensembl REST API, Compara (plants) |
| `data/genotypes/withdrawn_bradi_id_check.csv` | `15_map_orthologs.py` (network) | Ensembl REST API |
| `data/genotypes/gravitropism_candidate_genes.csv` | `03_candidate_genes.py` | curated list, anchor check, Compara orthologs |
| `data/su2023/su2023_degs.csv` | `16_su2023_degs.py` | Su et al. Supplemental Tables S2–S4 |
| `data/su2023/id_mapping_check.csv` | `16_su2023_degs.py --check-ids` (network) | Ensembl REST API |
| `tables/su2023_deg_counts.csv` | `16_su2023_degs.py` | `su2023_degs.csv` |
| `tables/candidate_spaceflight_degs.csv` | `16_su2023_degs.py` | `su2023_degs.csv`, curated list, Compara orthologs |
| `tables/astrograss_master_table.csv`, `docs/js/astrograss_database.json` | `10_build_astrograss_db.py` | candidate table, candidate DEG calls |
| `data/osdr/astrograss_studies_catalog.csv` | `10_build_astrograss_db.py` | study records typed into the script after checking each against OSDR, GEO, PRIDE, CrossRef and DataCite (14 Oct 2026) |
| `tables/video_catalog.csv`, `data/videos/video_catalog.json`, `docs/js/video_database.json` | `11_build_video_catalog.py` | video entries typed into the script (not re-checked against YouTube or CyVerse) |
| `tables/cyverse_brachypodium_series.csv`, `data/osdr/cyverse_image_series.json` | `12_cyverse_irods_sync.py` | series metadata typed into the script (frame counts not re-checked against the CyVerse collection) |
| `docs/tables/*.csv` | copies made by `10`, `11` and `12` | `tables/`, `data/osdr/` |

## Figures

| Figure | Produced by | From |
|---|---|---|
| `fig1_osd375_sample_design` | `09_generate_figures.py` | `data/osdr/experimental_design_summary.csv` |
| `fig2_astrograss_architecture` | `09_generate_figures.py` | schematic; labels typed into the script |
| `fig3_su2023_deg_counts` | `09_generate_figures.py` | `tables/su2023_deg_counts.csv` |
| `fig4_candidate_spaceflight_degs` | `09_generate_figures.py` | `tables/candidate_spaceflight_degs.csv` |
| `figS1_mission_experimental_design` | `09_generate_figures.py` | Su et al. Figures 1–2; timeline and conditions typed into the script from Su et al. (2023) and the OSD-375 ISA table |

`run_pipeline.py` copies the figures to `docs/assets/`.

## Manuscript

| File | Produced by |
|---|---|
| `manuscript/manuscript.pdf`, `manuscript/manuscript.docx` | `export_manuscript.py`, from `manuscript.tex` and `references.bib` |
| `manuscript/references.bib` | generated from CrossRef records on 14 Oct 2026; all 15 entries re-verified against CrossRef/OpenAlex |
