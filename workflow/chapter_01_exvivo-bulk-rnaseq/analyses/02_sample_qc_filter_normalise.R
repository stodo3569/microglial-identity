# ==============================================================================
# Script:   02_sample_qc_filter_normalise.R
# Chapter:  [1] — exvivo bulk rnaseq
# Step:     2 of [Y] in the CHAPTER 01 workflow
#           Preceded by: 01_tximport.R
#           Followed by: 03_[Y]
#
# Description:
#   Takes the merged gene-level count matrix produced by 01_tximport.R
#   and performs sample-level QC, multi-step filtering, and normalisation.
#   Steps are:
#     1. Align count matrix to metadata; impute sex; compute per-sample
#        QC metrics (library size, gene detection rate at CPM > 5).
#        → Checkpoint A saved here for pre-filter exploratory figures.
#     2. Apply sample-level inclusion criteria (minimum library size,
#        age annotation, neurological confounders, gene detection floor).
#     3. Filter lowly expressed genes (autosomal only; filterByExpr).
#        → Checkpoint B saved here for post-filter QC figures.
#     4. Multi-study outlier removal using a custom extension of
#        DaMiR.sampleFilt() with study-size-weighted inter-study
#        connectivity Z-score, within-study Z-score, small-study exemption,
#        and a final WGCNA standardised connectivity (Z.k) pass.
#     5. Second gene-expression filter on the retained sample set.
#     6. Study batch correction with limma::removeBatchEffect.
#     7. Study batch correction wtih unwanted and hidden variable identification (SVA).
#
# Inputs:
#   - CHAPTER_01/rds/txi_all_datasets.rds
#                                 Merged tximport object from 01_tximport.R
#   - common/exvivo_sample_metadata.tsv
#                                 Sample-level metadata; version-controlled in
#                                 the repository (not an RDS — loaded directly
#                                 with read_tsv())
#
# Outputs:
#   - CHAPTER_01/rds/exvivo_pre_filter.rds
#                                 List: $counts (raw), $metadata — after QC
#                                 metric computation, before any sample removal.
#                                 Loaded by [FIGURE SCRIPT FOR CHECKPOINT A].
#   - CHAPTER_01/rds/exvivo_post_sample_filter.rds
#                                 List: $counts, $metadata, $SE — after sample
#                                 and gene filters, before outlier removal.
#                                 Loaded by [FIGURE SCRIPT FOR CHECKPOINT B].
#   - CHAPTER_01/rds/exvivo_normalised.rds
#                                 List: $counts_filtered, $metadata_filtered,
#                                 $SE_norm, $qc_summary, $wgcna_connectivity —
#                                 final output; loaded by [FOLLOWING SCRIPT NAME].
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
# This script reads one input RDS and produces three output RDS objects
# (two intermediate checkpoints + one final output):
#
#   Input:
#     exvivo_counts.rds     [DESCRIPTION OF OBJECT]
#
#   Outputs:
#     [XX]a_exvivo_pre_filter.rds               Checkpoint A: pre-filter snapshot
#     [XX]b_exvivo_post_sample_filter.rds        Checkpoint B: post-filter snapshot
#     [XX]_exvivo_normalised.rds                 Final output for downstream scripts
#
#   Note: common/exvivo_sample_metadata.tsv is also read in this script but
#   is version-controlled in the repository and does not need to be downloaded.
#
# Three ways to proceed — set the flags below accordingly:
#
#   1) Skip this script entirely (DOWNLOAD_OUTPUT_RDS = TRUE):
#      Downloads all three output RDS objects from Zenodo and stops.
#      No raw quant.sf files or intermediate objects needed.
#
#   2) Load the input RDS from Zenodo (LOAD_INPUT_RDS = TRUE, default):
#      Downloads exvivo_counts.rds from Zenodo if not present
#      locally, then runs this script in full.
#
#   3) Run fully from local files (LOAD_INPUT_RDS = FALSE):
#      Expects exvivo_counts.rds to be present locally (either
#      produced by running 01_tximport.R, or placed manually).
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE <- "https://zenodo.org/records/XXXXXXX/files"

# ── Input RDS ────────────────────────────────────────────────────────────────
input_rds <- list(
  txi = "CHAPTER_01/rds/exvivo_counts.rds"   # [DESCRIPTION OF OBJECT]
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
  pre_filter  = "[CHAPTER_XX]/rds/[XX]a_exvivo_pre_filter.rds",
  post_filter = "[CHAPTER_XX]/rds/[XX]b_exvivo_post_sample_filter.rds",
  normalised  = "[CHAPTER_XX]/rds/[XX]_exvivo_normalised.rds"
)

# DOWNLOAD_OUTPUT_RDS: set to TRUE to download all output RDS from Zenodo
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
  message("All output RDS downloaded. Proceed to the next script.")
  stop("Stopping early — set DOWNLOAD_OUTPUT_RDS <- FALSE to run the full script.", call. = FALSE)
}

# SAVE_OUTPUT_RDS: set to TRUE (default) to save all output RDS to disk.
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
library(DESeq2)
library(sva)
library(WGCNA)

# ==============================================================================
# 1. Load metadata and compute per-sample QC metrics
# ==============================================================================

exvivo_metadata <- read_tsv("microglial-identity/common/exvivo_sample_metadata.tsv")

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

# ── Checkpoint A: pre-filter snapshot ────────────────────────────────────────
# Saved here so that [FIGURE SCRIPT FOR CHECKPOINT A] can load this object
# and generate figures showing the full, unfiltered sample distribution.
if (SAVE_OUTPUT_RDS) {
  saveRDS(
    list(counts = exvivo_counts, metadata = exvivo_metadata),
    file = output_rds$pre_filter
  )
  message("Saved checkpoint A: ", output_rds$pre_filter)
}

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

# ── Checkpoint B: post-sample-filter snapshot ─────────────────────────────────
# Saved here so that [FIGURE SCRIPT FOR CHECKPOINT B] can load this object
# and generate QC figures on the filtered-but-pre-outlier-removal dataset.
if (SAVE_OUTPUT_RDS) {
  saveRDS(
    list(counts = exvivo_counts, metadata = exvivo_metadata_filtered, SE = exvivo_SE),
    file = output_rds$post_filter
  )
  message("Saved checkpoint B: ", output_rds$post_filter)
}

# ==============================================================================
# 4. Normalise and remove outlier samples
# ==============================================================================
# Multi-study QC using the improved custom_DaMiR.sampleFilt():
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
#   4. WGCNA standardised connectivity (Z.k) — final pass:
#      Recomputed on surviving samples. Catches residual outliers not captured
#      by correlation-based filters.
#      Threshold: th.wgcna_z

exvivo_SE_norm <- DaMiR.normalization(exvivo_SE, minCounts = 0, hyper = "no")

exvivo_samplefilt_result <- custom_DaMiR.sampleFilt(
  data                  = exvivo_SE_norm,
  study_df              = exvivo_metadata_filtered,
  th.inter_study_z      = -1,
  th.within_study_z     = -1,
  th.wgcna_z            = -1,
  th.duplicate_cor      = 0.99,
  small_study_threshold = 3,
  type                  = "spearman"
)

exvivo_qc_summary  <- exvivo_samplefilt_result$sample_qc_summary
exvivo_wgcna_conn  <- exvivo_samplefilt_result$wgcna_connectivity
keep_sample        <- exvivo_samplefilt_result$keep_samples

exvivo_metadata_filtered <- rownames_to_column(exvivo_metadata_filtered, "name_id") %>%
  filter(name_id %in% keep_sample)
exvivo_counts <- exvivo_counts[, exvivo_metadata_filtered$name_id]

# ==============================================================================
# 5. Second gene filter on the retained sample set
# ==============================================================================

cpm_exvivo_counts <- cpm(exvivo_counts)
keep              <- filterByExpr(y = cpm_exvivo_counts, min.count = 1, min.prop = 0.95)
exvivo_counts     <- exvivo_counts[keep, ]

# ==============================================================================
# 6. Final normalisation
# ==============================================================================

exvivo_metadata_filtered <- column_to_rownames(exvivo_metadata_filtered, "name_id")
exvivo_SE_fill <- DaMiR.makeSE(round(exvivo_counts), exvivo_metadata_filtered)
exvivo_SE_norm <- DaMiR.normalization(exvivo_SE_fill, minCounts = 0, hyper = "no")

# ── Final output ──────────────────────────────────────────────────────────────
# Loaded by: [FOLLOWING SCRIPT NAME]
if (SAVE_OUTPUT_RDS) {
  saveRDS(
    list(
      counts_filtered    = exvivo_counts,
      metadata_filtered  = exvivo_metadata_filtered,
      SE_norm            = exvivo_SE_norm,
      qc_summary         = exvivo_qc_summary,
      wgcna_connectivity = exvivo_wgcna_conn
    ),
    file = output_rds$normalised
  )
  message("Saved final output: ", output_rds$normalised)
}