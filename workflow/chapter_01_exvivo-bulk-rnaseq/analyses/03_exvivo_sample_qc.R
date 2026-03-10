# ==============================================================================
# Script:   03_exvivo_sample_qc.R
# Chapter:  [1] — exvivo bulk rnaseq
# Step:     3 of [Y] in the chapter_01 workflow
#           Preceded by: 02_exvivo_sample_exploratory.R
#           Followed by: 04_[Y]
#
# Description:
#   Takes the pre-filter snapshot RDS produced by 02_exvivo_sample_exploratory.R
#   and applies sample- and gene-level QC to produce a clean object ready for
#   normalisation and downstream analysis.
#   Steps are:
#     2. Apply sample-level inclusion criteria (library size, age, tissue
#        type, gene detection rate).
#     3. Filter lowly expressed and non-autosomal genes.
#     4. Remove outlier samples using a multi-study QC approach
#        (inter-study connectivity, within-study Z-score, WGCNA Z.k) —
#        all metrics computed simultaneously on the full sample set before
#        any removal is applied.
#        → RDS object saved here for post-QC figures and downstream analysis.
#
# Inputs:
#   - rds/chapter_01/exvivo_sample_exploratory.rds
#                                 List: $counts, $metadata — pre-filter
#                                 snapshot produced by 02_exvivo_sample_exploratory.R
#
# Outputs:
#   - rds/chapter_01/exvivo_sample_qc.rds
#                                 List: $counts (filtered, autosomal),
#                                 $metadata (post-QC), $SE_norm (normalised
#                                 SummarizedExperiment), $qc_summary.
#                                 Loaded by [FIGURE SCRIPT FOR CHECKPOINT B].
#
# Environment:
#   All scripts in this repository are intended to run inside a Docker
#   container. The image used for the RStudio session is available at:
#     https://hub.docker.com/r/stodo3569/rstudio-server_microglial-identity
#   See docker/ in this repository for the Dockerfile and build instructions.
#   Common utility functions shared across scripts are sourced from
#   common/R/00_functions.R (see repository root).
# ==============================================================================

source("microglial-identity/common/R/00_functions.R")

# ==============================================================================
# RDS management
# ==============================================================================
# This script reads one input RDS and produces one output RDS object
# (Checkpoint B — post-QC snapshot for downstream analysis and figures):
#
#   Input:
#     rds/chapter_01/exvivo_sample_exploratory.rds  List: $counts, $metadata
#
#   Output:
#     rds/chapter_01/exvivo_sample_qc.rds           RDS checkpoint: post-QC snapshot
#
# Two ways to proceed — set the flags below accordingly:
#
#   1) Skip this script entirely (DOWNLOAD_OUTPUT_RDS = TRUE):
#      Downloads the output RDS object from Zenodo and stops.
#      No intermediate objects needed.
#
#   2) Load the input RDS from Zenodo (LOAD_INPUT_RDS = TRUE, default):
#      Downloads exvivo_sample_exploratory.rds from Zenodo if not present
#      locally, then runs this script in full.
#
#   3) Run fully from local files (LOAD_INPUT_RDS = FALSE):
#      Expects exvivo_sample_exploratory.rds to be present locally at
#      rds/chapter_01/exvivo_sample_exploratory.rds (produced by
#      running 02_exvivo_sample_exploratory.R).
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE <- "https://zenodo.org/records/XXXXXXX/files/chapter_01/rds"

# ── Input RDS ────────────────────────────────────────────────────────────────
input_rds <- list(
  pre_filter = "rds/chapter_01/exvivo_sample_exploratory.rds"  # List: $counts, $metadata
)

LOAD_INPUT_RDS <- TRUE

if (LOAD_INPUT_RDS) {
  for (rds_name in names(input_rds)) {
    local_path <- input_rds[[rds_name]]
    local_file <- basename(local_path)
    if (file.exists(local_path)) {
      message(local_file, " already exists locally — skipping download.")
    } else if (file.exists(local_file)) {
      message(local_file, " found in working directory — skipping download.")
      local_path <- local_file
    } else {
      url <- paste0(ZENODO_BASE, "/", local_file)
      message("Downloading ", local_file, " from Zenodo...")
      download.file(url, destfile = local_path, mode = "wb")
    }
  }
  checkpoint_a    <- readRDS(local_path)
  exvivo_counts   <- checkpoint_a$counts
  exvivo_metadata <- checkpoint_a$metadata
  remove(checkpoint_a)
}

# ── Output RDS ───────────────────────────────────────────────────────────────
output_rds <- list(
  post_qc = "rds/chapter_01/exvivo_sample_qc.rds"
)

# DOWNLOAD_OUTPUT_RDS: set to TRUE to download the output RDS from Zenodo
# instead of running this script. The script will stop after downloading.
DOWNLOAD_OUTPUT_RDS <- FALSE

if (DOWNLOAD_OUTPUT_RDS) {
  for (rds_name in names(output_rds)) {
    local_path <- output_rds[[rds_name]]
    if (!file.exists(local_path)) {
      url <- paste0(ZENODO_BASE, "/", basename(local_path))
      message("Downloading ", basename(local_path), " from Zenodo...")
      download.file(url, destfile = local_path, mode = "wb")
    } else {
      message(basename(local_path), " already exists locally — skipping download.")
    }
  }
  message("Output RDS downloaded. Proceed to the next script.")
  stop("Stopping early — set DOWNLOAD_OUTPUT_RDS <- FALSE to run the full script.", call. = FALSE)
}

# SAVE_OUTPUT_RDS: set to TRUE (default) to save the output RDS to disk.
# Set to FALSE to skip saving — useful when exploring interactively.
SAVE_OUTPUT_RDS <- TRUE

# ==============================================================================
# Set up working environment
# ==============================================================================

set.seed(123)

library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DaMiRseq)

# ==============================================================================
# 2. Apply sample-level inclusion criteria
# ==============================================================================

exvivo_metadata_filtered <- exvivo_metadata %>%
  dplyr::filter(total_reads >= 5000000) %>%
  dplyr::filter(!is.na(age)) %>%
  dplyr::filter(Control_tissue == TRUE) %>%
  dplyr::filter(sample_gene_count > 7999)
exvivo_counts <- exvivo_counts[, exvivo_metadata_filtered$name_id]

# ==============================================================================
# 3. Filter lowly expressed genes (autosomal only)
# ==============================================================================

rows_chrom     <- mapIds(org.Hs.eg.db, column = "CHR", keytype = "ENSEMBL",
                         keys = gsub("\\..*", "", rownames(exvivo_counts)), unique = TRUE)
rows_autosomal <- rows_chrom %in% 1:22
exvivo_counts  <- exvivo_counts[rows_autosomal, ]

cpm_exvivo_counts <- cpm(exvivo_counts)
keep              <- filterByExpr(y = cpm_exvivo_counts, min.count = 1, min.prop = 0.7)
exvivo_counts     <- exvivo_counts[keep, ]

# Prepare SummarizedExperiment for normalisation and QC
exvivo_metadata_filtered$class <- "class"
exvivo_metadata_filtered        <- column_to_rownames(exvivo_metadata_filtered, "name_id")
exvivo_SE <- DaMiR.makeSE(round(exvivo_counts), exvivo_metadata_filtered)

# ==============================================================================
# 4. Normalise and remove outlier samples
# ==============================================================================
# Multi-study QC using custom_DaMiR.sampleFilt() in simultaneous mode:
# all four metrics are computed on the full sample set before any removal,
# so every sample is evaluated against the same reference distribution.
#
#   1. Study-size-weighted inter-study connectivity (Z-scored):
#      Mean correlation to each OTHER study computed separately, then averaged
#      — giving every study an equal vote. Result is Z-scored globally.
#      Threshold: th.inter_study_z
#
#   2. Within-study Z-score:
#      Mean correlation to own study-mates, Z-scored WITHIN each study.
#      Threshold: th.within_study_z
#
#   3. Small-study exemption:
#      Studies with <= small_study_threshold samples are exempt from the
#      within-study criterion (too few samples for a stable Z-score).
#
#   4. WGCNA standardised connectivity (Z.k):
#      Computed on the full sample set alongside the correlation-based metrics.
#      Catches samples with globally low connectivity regardless of whether
#      they pass the study-level criteria.
#      Threshold: th.wgcna_z
#
#   All four pass/fail decisions are applied simultaneously in a single
#   removal pass. $qc_summary contains complete metrics for every sample,
#   including those that are removed.

exvivo_SE_norm <- DaMiR.normalization(exvivo_SE, minCounts = 0, hyper = "no")

exvivo_samplefilt_result <- custom_DaMiR.sampleFilt(
  data                  = exvivo_SE_norm,
  study_df              = exvivo_metadata_filtered,
  th.inter_study_z      = -1,
  th.within_study_z     = -1,
  th.wgcna_z            = -1,
  th.duplicate_cor      = 0.99,
  small_study_threshold = 3,
  type                  = "pearson"
)

exvivo_qc_summary <- exvivo_samplefilt_result$sample_qc_summary
keep_sample       <- exvivo_samplefilt_result$keep_samples

exvivo_metadata_filtered <- rownames_to_column(exvivo_metadata_filtered, "name_id") %>%
  filter(name_id %in% keep_sample)
exvivo_counts <- exvivo_counts[, exvivo_metadata_filtered$name_id]

 # ── RDS checkpoint: post-QC snapshot ──────────────────────────────────────────
# Saved here so that fig_exvivo_sample_qc.R can load this object and generate
# the multi-study QC panel (inter-study Z, within-study Z, WGCNA Z.k).
# $qc_summary contains complete metrics for all samples including removed ones,
# so the figure script can show the full distribution with outliers highlighted.
# Also serves as the primary input for downstream analysis scripts
# (SVA, differential expression).
if (SAVE_OUTPUT_RDS) {
  dir.create("rds/chapter_01", recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      counts     = exvivo_counts,
      metadata   = exvivo_metadata_filtered,
      SE_norm    = exvivo_SE_norm,
      qc_summary = exvivo_qc_summary
    ),
    file = output_rds$post_qc
  )
  message("Saved RDS object: ", output_rds$post_qc)
}
