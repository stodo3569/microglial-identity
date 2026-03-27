# ==============================================================================
# Script:   02_exvivo_sample_exploratory.R
# Chapter:  [1] — exvivo bulk rnaseq
# Step:     2 of [Y] in the chapter_01 workflow
#           Preceded by: 01_tximport.R
#           Followed by: 03_exvivo_sample_qc.R
#           Used by: fig_exvivo_exploratory_prefilter.R; tab_exvivo_exploratory_prefilter.R
#
# Description:
#   Takes the merged gene-level count matrix produced by 01_tximport.R
#   and performs sample-level QC to produce a pre-filter snapshot.
#   Steps are:
#     1. Align count matrix to metadata; impute sex; compute per-sample
#        QC metrics (library size, gene detection rate at CPM > 5).
#        → RDS object saved here for pre-filter exploratory figures.
#
# Inputs:
#   - rds/chapter_01/exvivo_counts.rds
#                                 Merged tximport object from 01_tximport.R
#   - common/exvivo_sample_metadata.tsv
#                                 Sample-level metadata; version-controlled in
#                                 the repository (not an RDS — loaded directly
#                                 with read_tsv())
#
# Outputs:
#   - rds/chapter_01/exvivo_sample_exploratory.rds
#                                 List: $counts (raw), $metadata — after QC
#                                 metric computation, before any sample removal.
#                                 Loaded by [FIGURE SCRIPT FOR CHECKPOINT A].
#
# Environment:
#   All scripts in this repository are intended to run inside a Docker
#   container. The image used for the RStudio session is available at:
#     https://hub.docker.com/r/stodo3569/rstudio-server_microglial-identity
#   See docker/ in this repository for the Dockerfile and build instructions.
#   Common utility functions shared across scripts are sourced from
#   common/R/00_functions.R (see repository root).
# ==============================================================================

source("common/R/00_functions.R")

# ==============================================================================
# RDS management
# ==============================================================================
# This script reads one input RDS and produces one output RDS object
# (Checkpoint A — pre-filter snapshot for exploratory figures):
#
#   Input:
#      rds/chapter_01/exvivo_counts.rds     [DESCRIPTION OF OBJECT]
#
#   Output:
#     rds/chapter_01/exvivo_sample_exploratory.rds  RDS checkpoint: pre-filter snapshot
#
#   Note: common/exvivo_sample_metadata.tsv is also read in this script but
#   is version-controlled in the repository and does not need to be downloaded.
#
# Two ways to proceed — set the flags below accordingly:
#
#   1) Skip this script entirely (DOWNLOAD_OUTPUT_RDS = TRUE):
#      Downloads the output RDS object from Zenodo and stops.
#      No raw quant.sf files or intermediate objects needed.
#
#   2) Load the input RDS from Zenodo (LOAD_INPUT_RDS = TRUE, default):
#      Downloads exvivo_counts.rds from Zenodo if not present
#      locally, then runs this script in full.
#
#   3) Run fully from local files (LOAD_INPUT_RDS = FALSE):
#      Expects exvivo_counts.rds to be present locally at
#      rds/chapter_01/exvivo_counts.rds (either
#      produced by running 01_tximport.R, or placed manually).
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE <- "https://zenodo.org/records/XXXXXXX/files/chapter_01/rds"

# ── Input RDS ────────────────────────────────────────────────────────────────
input_rds <- list(
  txi = "rds/chapter_01/exvivo_counts.rds"   # [DESCRIPTION OF OBJECT]
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
  exvivo_counts <- readRDS(local_path)
}

# ── Output RDS ───────────────────────────────────────────────────────────────
output_rds <- list(
  pre_filter = "rds/chapter_01/exvivo_sample_exploratory.rds"
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
# 1. Load metadata and compute per-sample QC metrics
# ==============================================================================

exvivo_metadata <- read_tsv("common/exvivo_sample_metadata.tsv")

# Align with the counts matrix
exvivo_metadata <- filter(exvivo_metadata, name_id %in% colnames(exvivo_counts))
exvivo_counts   <- exvivo_counts[, exvivo_metadata$name_id]

# Impute sex
exvivo_metadata <- SexDetect(exvivo_counts, exvivo_metadata)

# Add total reads mapped variable
exvivo_metadata <- add_column(exvivo_metadata, total_reads = colSums(exvivo_counts))

# Calculate gene count per sample (CPM > 5)
tmp_cpm_count <- cpm(exvivo_counts)
for (index in 1:ncol(tmp_cpm_count)) {
  sample_name       <- colnames(tmp_cpm_count)[index]
  sample_gene_count <- sum(tmp_cpm_count[, index] > 5)
  message("Gene count for ", sample_name, " is ", sample_gene_count)
  exvivo_metadata[
    which(exvivo_metadata$geo_accession == sample_name |
            exvivo_metadata$sample_title  == sample_name),
    "sample_gene_count"
  ] <- sample_gene_count
}
remove(tmp_cpm_count)

# ── RDS checkpoint: pre-filter snapshot ────────────────────────────────────────
# Saved here so that fig_exvivo_exploratory_prefilter.R and
# tab_exvivo_exploratory_prefilter.R can load this object
# and generate figures showing the full, unfiltered sample distribution.
# or 03_exvivo_sample_qc.R to proceed with the quality control

if (SAVE_OUTPUT_RDS) {
  dir.create("rds/chapter_01", recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(counts = exvivo_counts, metadata = exvivo_metadata),
    file = output_rds$pre_filter
  )
  message("Saved RDS object: ", output_rds$pre_filter)
}
