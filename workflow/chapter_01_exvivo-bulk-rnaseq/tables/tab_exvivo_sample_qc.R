# ==============================================================================
# Script:   tab_exvivo_sample_qc.R
# Chapter:  [1] — exvivo bulk rnaseq
# Step:     Table script — loads RDS checkpoint from 03_exvivo_sample_qc.R
#           Preceded by: 03_exvivo_sample_qc.R
#           Followed by: n/a (terminal output script)
#
# Description:
#   Generates a per-study QC filtering summary table from the post-QC
#   snapshot of the ex vivo sample set. The table shows, for each study, how
#   many samples passed or failed each of the simultaneous multi-study QC
#   criteria applied by custom_DaMiR.sampleFilt() in 03_exvivo_sample_qc.R.
#   A summary row gives totals across all studies.
#
#   QC criteria summarised:
#     - Inter-study connectivity Z-score  >= -1
#     - Within-study connectivity Z-score >= -1 (studies with n <= 3 exempt)
#     - WGCNA standardised connectivity   >= -1
#     - Not flagged as a near-duplicate   (pairwise r < 0.99)
#
#   Inputs:
#     - rds/chapter_01/exvivo_sample_qc.rds
#                                 List: $counts, $metadata, $SE_norm,
#                                       $qc_summary —
#                                 post-QC checkpoint from 03_exvivo_sample_qc.R
#
#   Outputs:
#     - outputs/chapter_01/tables/tab_exvivo_sample_qc.html
#                                 gt table summarising per-study pass / fail
#                                 counts for each QC criterion; exported as a
#                                 self-contained HTML file.
#
# Environment:
#   All scripts in this repository are intended to run inside a Docker
#   container. The image used for the RStudio session is available at:
#     https://hub.docker.com/r/stodo3569/rstudio-server_microglial-identity
#   See docker/ in this repository for the Dockerfile and build instructions.
#   Common utility functions are sourced from
#   common/R/00_functions.R (see repository root).
# ==============================================================================

source("microglial-identity/common/R/00_functions.R")

# ==============================================================================
# RDS management
# ==============================================================================
# This script loads one input RDS and produces no new RDS objects.
# SAVE_OUTPUT_RDS is always FALSE in table scripts.
#
#   Input:
#     rds/chapter_01/exvivo_sample_qc.rds     List: $counts, $metadata,
#                                                    $SE_norm, $qc_summary
#
# Two ways to proceed — set the flag below accordingly:
#
#   1) Load from Zenodo (LOAD_INPUT_RDS = TRUE, default):
#      Downloads exvivo_sample_qc.rds from Zenodo if not present
#      locally, then runs this script in full.
#
#   2) Run fully from local files (LOAD_INPUT_RDS = FALSE):
#      Expects exvivo_sample_qc.rds to be present locally at
#      rds/chapter_01/exvivo_sample_qc.rds (produced by running
#      03_exvivo_sample_qc.R).
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE     <- "https://zenodo.org/records/XXXXXXX/files/chapter_01/rds"
LOAD_INPUT_RDS  <- TRUE
SAVE_OUTPUT_RDS <- FALSE   # Table scripts never save new RDS objects

input_rds <- list(
  post_qc = "rds/chapter_01/exvivo_sample_qc.rds"
)

if (LOAD_INPUT_RDS) {
  for (rds_name in names(input_rds)) {
    local_path <- input_rds[[rds_name]]
    if (file.exists(local_path)) {
      message(basename(local_path), " already exists locally — skipping download.")
    } else {
      url <- paste0(ZENODO_BASE, "/", basename(local_path))
      message("Downloading ", basename(local_path), " from Zenodo...")
      download.file(url, destfile = local_path, mode = "wb")
    }
  }
  dat <- readRDS(input_rds$post_qc)
}

exvivo_qc_summary <- dat$qc_summary

# ==============================================================================
# Set up working environment
# ==============================================================================

library(tidyverse)
library(gt)

# ==============================================================================
# Table 1.02 — Per-study sample counts: multi-study QC filter outcomes
# ==============================================================================
# Summarises per study how many samples passed or failed each QC criterion
# individually, and how many passed all criteria simultaneously. Because all
# four metrics are computed on the full sample set before any removal, every
# sample contributes to every column. Studies with n <= 3 are exempt from the
# within-study Z-score criterion; these are marked with an asterisk.
# A grand-total row is appended and rendered in bold.
# ==============================================================================

strip_accession <- function(x) sub("_(GSE|SRP|EGAD)\\d+$", "", x)

tab_data <- exvivo_qc_summary %>%
  mutate(study = strip_accession(study)) %>%
  group_by(study) %>%
  summarise(
    total          = n(),
    fail_inter     = sum(!pass_correlation_filters & pass_dedup_filter & pass_wgcna_filter, na.rm = TRUE),
    fail_within    = sum(pass_correlation_filters  & !pass_dedup_filter & pass_wgcna_filter & !is.na(within_study_z), na.rm = TRUE),
    fail_duplicate = sum(removed_as_duplicate, na.rm = TRUE),
    fail_wgcna     = sum(!pass_wgcna_filter & pass_correlation_filters & pass_dedup_filter, na.rm = TRUE),
    pass           = sum(pass_all_filters, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  add_row(
    study          = "Total",
    total          = sum(.$total),
    fail_inter     = sum(.$fail_inter),
    fail_within    = sum(.$fail_within),
    fail_duplicate = sum(.$fail_duplicate),
    fail_wgcna     = sum(.$fail_wgcna),
    pass           = sum(.$pass)
  )

tab_gt <- tab_data %>%
  gt() %>%
  cols_label(
    study          = "Study",
    total          = "Total",
    fail_inter     = "Failed inter-study Z",
    fail_within    = "Failed within-study Z",
    fail_duplicate = "Removed as duplicate",
    fail_wgcna     = "Failed WGCNA Z.k",
    pass           = "Passed all criteria"
  ) %>%
  tab_spanner(
    label   = "Removed by criterion",
    columns = c(fail_inter, fail_within, fail_duplicate, fail_wgcna)
  ) %>%
  fmt_number(
    columns  = c(total, fail_inter, fail_within, fail_duplicate, fail_wgcna, pass),
    decimals = 0
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = nrow(tab_data))
  ) %>%
  tab_footnote(
    footnote  = "Studies with n \u2264 3 samples are exempt from the within-study Z-score criterion.",
    locations = cells_column_labels(columns = fail_within)
  )

# ==============================================================================
# Save output
# ==============================================================================

dir.create(
  "outputs/chapter_01/tables",
  recursive    = TRUE,
  showWarnings = FALSE
)

gtsave(
  tab_gt,
  filename = "outputs/chapter_01/tables/tab_exvivo_sample_qc.html"
)

message("Saved: outputs/chapter_01/tables/tab_exvivo_sample_qc.html")
