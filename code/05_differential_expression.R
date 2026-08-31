#!/usr/bin/env Rscript
# ==============================================================================
# 05_differential_expression.R
#
# Differential expression analysis of OSD-375 Brachypodium distachyon
# spaceflight RNA-Seq data using DESeq2.
#
# Study: APEX-06 (SpaceX CRS-14 -> ISS, April 2018)
# Ecotypes: Bd21 (reference), Bd21-3 (transformable), Gaz8 (Turkish)
# Tissues: Root, Shoot
# Conditions: Spaceflight (ISS VEGGIE) vs Ground Control
#
# Reference: Su et al. (2023) Life 13(3):633. DOI: 10.3390/life13030633
# Data source: NASA OSDR OSD-375 (DOI: 10.26030/2x6b-3v89)
#
# Usage:
#   Rscript code/05_differential_expression.R [--input-dir data/osdr] [--output-dir tables]
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(EnhancedVolcano)
  library(optparse)
})

# --- CLI Arguments -----------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input-dir"), type = "character",
              default = "data/osdr",
              help = "Directory containing OSD-375 processed count matrices"),
  make_option(c("-o", "--output-dir"), type = "character",
              default = "tables",
              help = "Output directory for DEG tables"),
  make_option(c("-f", "--figures-dir"), type = "character",
              default = "figures",
              help = "Output directory for volcano/heatmap plots"),
  make_option(c("-p", "--padj-cutoff"), type = "double",
              default = 0.05,
              help = "Adjusted p-value cutoff for significance"),
  make_option(c("-l", "--lfc-cutoff"), type = "double",
              default = 0.5,
              help = "Log2 fold-change cutoff for significance")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Configuration -----------------------------------------------------------
ECOTYPES <- c("Bd21", "Bd21-3", "Gaz8")
TISSUES  <- c("root", "shoot")

cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
cat("Brachypodium OSD-375 Differential Expression Analysis\n")
cat("APEX-06 | ISS VEGGIE | Spaceflight vs Ground Control\n")
cat("=" %>% rep(70) %>% paste(collapse = ""), "\n\n")

# --- Create output directories -----------------------------------------------
dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`figures-dir`, recursive = TRUE, showWarnings = FALSE)

# --- Helper Functions ---------------------------------------------------------

#' Run DESeq2 on a count matrix for one ecotype x tissue combination
#'
#' @param counts_df Data frame with gene IDs as rownames, samples as columns
#' @param sample_info Data frame with sample metadata (condition column)
#' @param ecotype Character: ecotype name (Bd21, Bd21-3, Gaz8)
#' @param tissue Character: tissue type (root, shoot)
#' @param padj_cutoff Numeric: adjusted p-value cutoff
#' @param lfc_cutoff Numeric: log2 fold-change cutoff
#' @return Data frame of DESeq2 results
run_deseq2 <- function(counts_df, sample_info, ecotype, tissue,
                       padj_cutoff = 0.05, lfc_cutoff = 0.5) {
  cat(sprintf("\n--- Processing: %s %s ---\n", ecotype, tissue))

  # Subset to this ecotype x tissue
  samples <- sample_info %>%
    filter(accession == ecotype, tissue_type == tolower(tissue)) %>%
    pull(sample_id)

  if (length(samples) == 0) {
    cat(sprintf("  WARNING: No samples found for %s %s. Skipping.\n",
                ecotype, tissue))
    return(NULL)
  }

  # Subset count matrix
  counts_sub <- counts_df[, colnames(counts_df) %in% samples, drop = FALSE]
  meta_sub <- sample_info %>%
    filter(sample_id %in% colnames(counts_sub)) %>%
    mutate(condition = factor(condition, levels = c("ground_control", "spaceflight")))

  # Ensure matching order
  counts_sub <- counts_sub[, match(meta_sub$sample_id, colnames(counts_sub))]

  cat(sprintf("  Samples: %d spaceflight, %d ground control\n",
              sum(meta_sub$condition == "spaceflight"),
              sum(meta_sub$condition == "ground_control")))

  # Filter low-count genes
  keep <- rowSums(counts_sub >= 10) >= 2
  counts_sub <- counts_sub[keep, ]
  cat(sprintf("  Genes after filtering: %d\n", nrow(counts_sub)))

  # DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData  = round(counts_sub),
    colData    = meta_sub,
    design     = ~ condition
  )

  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "spaceflight", "ground_control"))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj) %>%
    mutate(
      ecotype = ecotype,
      tissue  = tissue,
      significant = padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff
    )

  n_up   <- sum(res_df$significant & res_df$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res_df$significant & res_df$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("  DEGs (padj < %.2f, |LFC| > %.1f): %d up, %d down\n",
              padj_cutoff, lfc_cutoff, n_up, n_down))

  return(res_df)
}

#' Generate volcano plot for DEG results
#'
#' @param res_df DESeq2 results data frame
#' @param ecotype Character: ecotype name
#' @param tissue Character: tissue type
#' @param figures_dir Character: output directory
make_volcano <- function(res_df, ecotype, tissue, figures_dir,
                         padj_cutoff = 0.05, lfc_cutoff = 0.5) {
  outfile <- file.path(figures_dir,
                       sprintf("volcano_%s_%s.png", tolower(ecotype), tissue))

  p <- EnhancedVolcano(
    res_df,
    lab         = res_df$gene_id,
    x           = "log2FoldChange",
    y           = "padj",
    title       = sprintf("OSD-375: %s %s (Spaceflight vs Ground)", ecotype, tissue),
    subtitle    = "APEX-06 / ISS VEGGIE / DESeq2",
    pCutoff     = padj_cutoff,
    FCcutoff    = lfc_cutoff,
    pointSize   = 1.5,
    labSize     = 3.0,
    col         = c("grey30", "forestgreen", "royalblue", "red2"),
    legendLabels = c("NS", "|LFC|", "padj", "padj & |LFC|")
  )

  ggsave(outfile, p, width = 10, height = 8, dpi = 300)
  cat(sprintf("  Saved: %s\n", outfile))
}

# --- Main Analysis -----------------------------------------------------------

# Check for count matrix file
count_files <- list.files(opt$`input-dir`,
                          pattern = "counts?\\.csv$|unnormalized.*counts",
                          full.names = TRUE, recursive = TRUE)

if (length(count_files) == 0) {
  cat("\n")
  cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
  cat("NOTE: No count matrix files found in ", opt$`input-dir`, "\n")
  cat("Expected files: *counts.csv or *unnormalized*counts*\n\n")
  cat("To obtain the data, run:\n")
  cat("  python code/01_download_osdr_osd375.py\n")
  cat("  python code/02_parse_osdr_metadata.py\n\n")
  cat("Generating example output structure for downstream scripts...\n")
  cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")

  # Generate expected output schema for downstream scripts
  schema_df <- tibble(
    gene_id          = character(),
    baseMean         = double(),
    log2FoldChange   = double(),
    lfcSE            = double(),
    stat             = double(),
    pvalue           = double(),
    padj             = double(),
    ecotype          = character(),
    tissue           = character(),
    significant      = logical()
  )

  for (eco in ECOTYPES) {
    for (tis in TISSUES) {
      outfile <- file.path(opt$`output-dir`,
                           sprintf("deg_%s_%s.csv", tolower(gsub("-", "", eco)), tis))
      write_csv(schema_df, outfile)
      cat(sprintf("  Created schema: %s\n", outfile))
    }
  }

  # Summary file
  summary_df <- tibble(
    ecotype   = rep(ECOTYPES, each = 2),
    tissue    = rep(TISSUES, 3),
    n_deg_up  = NA_integer_,
    n_deg_down = NA_integer_,
    n_total   = NA_integer_,
    note      = "Awaiting count matrix data from OSDR download"
  )
  write_csv(summary_df, file.path(opt$`output-dir`, "deg_summary.csv"))
  cat(sprintf("  Created summary: %s/deg_summary.csv\n", opt$`output-dir`))

} else {
  cat(sprintf("Found count matrix: %s\n", count_files[1]))

  # Load count matrix
  counts_df <- read_csv(count_files[1], show_col_types = FALSE)

  # Try to set gene IDs as rownames
  if ("gene_id" %in% colnames(counts_df)) {
    counts_df <- counts_df %>% column_to_rownames("gene_id")
  } else if ("Gene" %in% colnames(counts_df)) {
    counts_df <- counts_df %>% column_to_rownames("Gene")
  }

  # Load sample metadata
  meta_file <- file.path(opt$`input-dir`, "sample_metadata.csv")
  if (!file.exists(meta_file)) {
    cat("ERROR: Sample metadata not found at ", meta_file, "\n")
    cat("Run: python code/02_parse_osdr_metadata.py\n")
    quit(status = 1)
  }

  sample_info <- read_csv(meta_file, show_col_types = FALSE)

  # Run DESeq2 for each ecotype x tissue
  all_results <- list()
  summary_rows <- list()

  for (eco in ECOTYPES) {
    for (tis in TISSUES) {
      res <- run_deseq2(counts_df, sample_info, eco, tis,
                        opt$`padj-cutoff`, opt$`lfc-cutoff`)

      if (!is.null(res)) {
        # Save individual results
        outfile <- file.path(opt$`output-dir`,
                             sprintf("deg_%s_%s.csv",
                                     tolower(gsub("-", "", eco)), tis))
        write_csv(res, outfile)
        cat(sprintf("  Saved: %s\n", outfile))

        # Volcano plot
        tryCatch(
          make_volcano(res, eco, tis, opt$`figures-dir`,
                       opt$`padj-cutoff`, opt$`lfc-cutoff`),
          error = function(e) cat(sprintf("  Warning: volcano plot failed: %s\n", e$message))
        )

        all_results[[paste(eco, tis, sep = "_")]] <- res

        # Summary stats
        summary_rows[[length(summary_rows) + 1]] <- tibble(
          ecotype    = eco,
          tissue     = tis,
          n_deg_up   = sum(res$significant & res$log2FoldChange > 0, na.rm = TRUE),
          n_deg_down = sum(res$significant & res$log2FoldChange < 0, na.rm = TRUE),
          n_total    = sum(res$significant, na.rm = TRUE)
        )
      }
    }
  }

  # Combined summary
  summary_df <- bind_rows(summary_rows)
  write_csv(summary_df, file.path(opt$`output-dir`, "deg_summary.csv"))

  cat("\n")
  cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
  cat("DEG Summary (Spaceflight vs Ground Control):\n")
  cat("=" %>% rep(70) %>% paste(collapse = ""), "\n")
  print(summary_df, n = Inf)

  # Combined results for all comparisons
  all_degs <- bind_rows(all_results)
  write_csv(all_degs, file.path(opt$`output-dir`, "deg_all_comparisons.csv"))

  # Heatmap of top DEGs across all comparisons
  cat("\nGenerating cross-ecotype heatmap...\n")
  tryCatch({
    top_genes <- all_degs %>%
      filter(significant == TRUE) %>%
      group_by(gene_id) %>%
      summarize(mean_lfc = mean(abs(log2FoldChange), na.rm = TRUE),
                n_comparisons = n()) %>%
      arrange(desc(n_comparisons), desc(mean_lfc)) %>%
      head(50)

    heatmap_data <- all_degs %>%
      filter(gene_id %in% top_genes$gene_id) %>%
      select(gene_id, ecotype, tissue, log2FoldChange) %>%
      mutate(label = paste(ecotype, tissue, sep = "_")) %>%
      select(gene_id, label, log2FoldChange) %>%
      pivot_wider(names_from = label, values_from = log2FoldChange) %>%
      column_to_rownames("gene_id")

    png(file.path(opt$`figures-dir`, "heatmap_top_degs.png"),
        width = 12, height = 16, units = "in", res = 300)
    pheatmap(as.matrix(heatmap_data),
             cluster_rows = TRUE, cluster_cols = TRUE,
             main = "OSD-375: Top 50 DEGs across Ecotypes & Tissues",
             fontsize_row = 6)
    dev.off()
    cat("  Saved: figures/heatmap_top_degs.png\n")
  }, error = function(e) {
    cat(sprintf("  Warning: heatmap generation failed: %s\n", e$message))
  })
}

cat("\n✓ Differential expression analysis complete.\n")
