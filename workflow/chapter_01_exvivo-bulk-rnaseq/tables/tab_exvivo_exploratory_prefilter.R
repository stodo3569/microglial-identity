# ==============================================================================
# Script:   tab_exvivo_exploratory_prefilter.R
# Chapter:  [1] — exvivo bulk rnaseq
# Step:     Table script — loads RDS checkpoint from 02_exvivo_sample_exploratory.R
#           Preceded by: 02_exvivo_sample_exploratory.R
#           Followed by: n/a (terminal output script)
#
# Description:
#   Generates a per-study QC filtering summary table from the pre-filter
#   snapshot of the ex vivo sample set. The table shows, for each study, how
#   many samples passed or failed the full set of inclusion criteria applied in
#   02_exvivo_sample_exploratory.R. A summary row gives totals across all studies.
#
#   Inclusion criteria applied:
#     - Mapping rate        > 0 %
#     - Total mapped reads  > 5,000,000
#     - Gene count (CPM>5)  > 8,000
#     - Donor age           > 0 (non-missing)
#     - Control tissue      == TRUE
#
#   Inputs:
#     - rds/chapter_01/exvivo_sample_exploratory.rds
#                                 List: $counts (raw count matrix),
#                                       $metadata (sample-level metadata) —
#                                 pre-filter rds object from 02_exvivo_sample_exploratory.R
#
#   Outputs:
#     - outputs/chapter_01/tables/tab_exvivo_exploratory_prefilter.html
#                                 gt table summarising per-study pass / fail
#                                 counts; exported as a self-contained HTML file.
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
#     rds/chapter_01/exvivo_sample_exploratory.rds     List: $counts, $metadata
#
# Two ways to proceed — set the flag below accordingly:
#
#   1) Load from Zenodo (LOAD_INPUT_RDS = TRUE, default):
#      Downloads exvivo_sample_exploratory.rds from Zenodo if not present
#      locally, then runs this script in full.
#
#   2) Run fully from local files (LOAD_INPUT_RDS = FALSE):
#      Expects exvivo_sample_exploratory.rds to be present locally at
#      rds/chapter_01/exvivo_sample_exploratory.rds (produced
#      by running 02_exvivo_sample_exploratory.R).
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE     <- "https://zenodo.org/records/XXXXXXX/files/chapter_01/rds"
LOAD_INPUT_RDS  <- TRUE
SAVE_OUTPUT_RDS <- FALSE   # Table scripts never save new RDS objects

input_rds <- list(
  pre_filter = "rds/chapter_01/exvivo_sample_exploratory.rds"
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
  dat <- readRDS(input_rds$pre_filter)
}

exvivo_counts   <- dat$counts
exvivo_metadata <- dat$metadata

# ==============================================================================
# Set up working environment
# ==============================================================================

library(tidyverse)
library(gt)

# ==============================================================================
# Table 1.01 — Per-study sample counts: pass / fail QC inclusion criteria
# ==============================================================================
# Summarises, per study, how many samples passed all five inclusion criteria
# simultaneously. A sample fails the table if it does not meet any single
# criterion. A grand-total row is appended and rendered in bold.
#
# Inclusion criteria mirror those applied in 02_exvivo_sample_exploratory.R:
#   mapping_rate     > 60
#   total_reads      > 4,999,999   (i.e. >= 5,000,000)
#   sample_gene_count > 7,999      (i.e. >= 8,000)
#   age              > 19 and non-missing
#   Control_tissue   == TRUE
# ==============================================================================

tab_data <- exvivo_metadata %>%
  mutate(
    passed_all_criteria =
      (mapping_rate      > 0)       &
      (total_reads       > 4999999)  &
      (sample_gene_count > 7999)     &
      (age               > 0)       &
      (!is.na(age))                  &
      (Control_tissue    == TRUE)
  ) %>%
  group_by(study) %>%
  summarise(
    total = n(),
    pass  = sum(passed_all_criteria),
    fail  = total - pass
  ) %>%
  add_row(
    study = "Total",
    total = sum(.$total),
    pass  = sum(.$pass),
    fail  = sum(.$fail)
  )

tab_gt <- tab_data %>%
  gt() %>%
  cols_label(
    study = "Study",
    total = "Total Samples",
    pass  = "Passed All Criteria",
    fail  = "Failed At Least One Criterion"
  ) %>%
  fmt_number(
    columns  = c(total, pass, fail),
    decimals = 0
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = nrow(tab_data))
  )

# ==============================================================================
# Save output
# ==============================================================================

dir.create(
  "outputs/chapter_01/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

gtsave(
  tab_gt,
  filename = "outputs/chapter_01/tables/tab_exvivo_exploratory_prefilter.html"
)

message("Saved: outputs/chapter_01/tables/tab_exvivo_exploratory_prefilter.html")
