# ==============================================================================
# Script:   04_exvivo_correction.R
# Chapter:  [1] — exvivo bulk rnaseq
# Step:     4 of [Y] in the chapter_01 workflow
#           Preceded by: 03_exvivo_sample_qc.R
#           Followed by: 05_[Y]
#
# Description:
#   Takes the post-QC snapshot RDS produced by 03_exvivo_sample_qc.R and
#   generates three normalised/batch-corrected expression matrices for
#   downstream analysis and visualisation:
#
#     1. VST — variance-stabilising transformation (DESeq2). Retains all post-QC
#        samples; used for unsupervised exploration and visualisation.
#
#     2. Study batch correction (limma) — removeBatchEffect() applied to
#        vst-normalised counts. Regresses out discrete study labels
#        while preserving age and sex effects.
#
#     3. SVA correction — surrogate variable analysis applied to  vst-normalised counts.
#        SVs are estimated via DaMiR.SV_modified()
#        and regressed out with removeBatchEffect(), preserving
#        the same design covariates (age, sex).
#
# Inputs:
#   - rds/chapter_01/exvivo_sample_qc.rds
#                                 List: $counts, $metadata, $SE_norm,
#                                 $qc_summary — post-QC snapshot produced
#                                 by 03_exvivo_sample_qc.R
#
# Outputs:
#   - rds/chapter_01/exvivo_counts_corrected.rds
#                                 List: $vst (VST + quantile-normalised +
#                                 scaled matrix), $sva_adjusted (SVA-corrected
#                                 log-CPM matrix), $study_adjusted (limma
#                                 study-batch-corrected log-CPM matrix),
#                                 $metadata_full (all post-QC samples),
#                                 $metadata_filtered (study-filtered samples
#                                 used for SVA and study correction).
#                                 Loaded by [FIGURE SCRIPT FOR CHECKPOINT C].
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
# (Checkpoint C — three normalised expression matrices for downstream use):
#
#   Input:
#     rds/chapter_01/exvivo_sample_qc.rds         List: $counts, $metadata,
#                                                  $SE_norm, $qc_summary
#
#   Output:
#     rds/chapter_01/exvivo_counts_corrected.rds      RDS checkpoint: normalised matrices
#
# Two ways to proceed — set the flags below accordingly:
#
#   1) Skip this script entirely (DOWNLOAD_OUTPUT_RDS = TRUE):
#      Downloads the output RDS object from Zenodo and stops.
#      No intermediate objects needed.
#
#   2) Load the input RDS from Zenodo (LOAD_INPUT_RDS = TRUE, default):
#      Downloads exvivo_sample_qc.rds from Zenodo if not present locally,
#      then runs this script in full.
#
#   3) Run fully from local files (LOAD_INPUT_RDS = FALSE):
#      Expects exvivo_sample_qc.rds to be present locally at
#      rds/chapter_01/exvivo_sample_qc.rds (produced by running
#      03_exvivo_sample_qc.R).
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE <- "https://zenodo.org/records/XXXXXXX/files/chapter_01/rds"

# ── Input RDS ────────────────────────────────────────────────────────────────
input_rds <- list(
  post_qc = "rds/chapter_01/exvivo_sample_qc.rds"  # List: $counts, $metadata, $SE_norm, $qc_summary
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
  checkpoint_b         <- readRDS(local_path)
  exvivo_counts        <- checkpoint_b$counts
  exvivo_metadata      <- checkpoint_b$metadata
  remove(checkpoint_b)
}

# ── Output RDS ───────────────────────────────────────────────────────────────
output_rds <- list(
  normalisation = "rds/chapter_01/exvivo_normalisation.rds"
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
library(DESeq2)
library(edgeR)
library(limma)
library(sva)
library(DaMiRseq)
library(SummarizedExperiment)

# ==============================================================================
# Studies to exclude, if needed
# ==============================================================================
# Identify studies that are too small or otherwise unsuitable for batch
# correction modelling. These are excluded from `exvivo_counts_filtered`
# (used for both SVA and limma study correction)

exclude_studies <- c()  # e.g. c("study_X", "study_Y") — populate as needed

# ==============================================================================
# Second lowly-expressed-gene filter (shared across all three methods)
# ==============================================================================
# Applied once here on all post-QC samples before any normalisation or
# study filtering. More stringent than the QC-step filter (min.prop = 0.95),
# requiring detectable expression in at least 95 % of samples. The resulting
# gene set is used as-is by all three downstream branches so that VST, limma
# study correction, and SVA operate on an identical set of genes.

cpm_exvivo_counts <- cpm(exvivo_counts)
keep              <- filterByExpr(y = cpm_exvivo_counts, min.count = 1, min.prop = 0.95)
exvivo_counts     <- exvivo_counts[keep, ]

# ── Filter samples by study ───────────────────────────────────────────────────
samples_to_keep <- exvivo_metadata %>%
  filter(!study %in% exclude_studies) %>%
  pull(name_id)

exvivo_counts <- exvivo_counts[, colnames(exvivo_counts) %in% samples_to_keep]
exvivo_metadata <- exvivo_metadata %>%
  filter(name_id %in% samples_to_keep)
rownames(exvivo_metadata) <- exvivo_metadata$name_id

# ==============================================================================
# 1. VST: variance-stabilising transformation
# ==============================================================================

exvivo_counts_vst <- vst(round(exvivo_counts))

# ==============================================================================
# 2 & 3. Study-filtered normalisation for limma and SVA correction
# ==============================================================================

# ── Shared design matrix ──────────────────────────────────────────────────────
exvivo_mod  <- model.matrix(~age + sex_imputed, data = exvivo_metadata)
exvivo_mod0 <- exvivo_mod[, 1]

# ==============================================================================
# 2. Study batch correction (limma removeBatchEffect)
# ==============================================================================

exvivo_study_adjusted <- removeBatchEffect(
  exvivo_counts_vst,
  batch  = exvivo_metadata$study,
  design = exvivo_mod
)

# ==============================================================================
# 3. SVA correction
# ==============================================================================
# SVA is run on the VST matrix subsetted to
# study-filtered samples. sva::num.sv() estimates the number of surrogate
# variables directly from the VST matrix. SVs are regressed out alongside the design
# covariates using removeBatchEffect().

exvivo_n.sv <- num.sv(exvivo_counts_vst, exvivo_mod, method = "be")

exvivo_sva <- sva(
  exvivo_counts_vst,
  exvivo_mod,
  exvivo_mod0,
  n.sv = exvivo_n.sv
)$sv

exvivo_sva_adjusted <- removeBatchEffect(
  exvivo_counts_vst,
  covariates = exvivo_sva,
  design     = exvivo_mod
)

# ── RDS checkpoint: normalised matrices ────────────────────────────────────────
# Saved here so that downstream analysis scripts and figure scripts can load
# a single object containing all three expression matrices:
#   $vst            — VST (all post-QC samples)
#   $study_adjusted — limma study-batch-corrected log-CPM 
#   $sva_adjusted   — SVA-corrected log-CPM 
#   $metadata_full  — metadata for all post-QC samples (matches $vst columns)

if (SAVE_OUTPUT_RDS) {
  dir.create("rds/chapter_01", recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      vst                = exvivo_counts_vst,
      study_adjusted     = exvivo_study_adjusted,
      sva_adjusted       = exvivo_sva_adjusted,
      metadata_full      = exvivo_metadata
    ),
    file = output_rds$normalisation
  )
  message("Saved RDS object: ", output_rds$normalisation)
}
