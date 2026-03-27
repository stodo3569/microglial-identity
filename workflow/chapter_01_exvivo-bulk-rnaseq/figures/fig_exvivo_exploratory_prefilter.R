# ==============================================================================
# Script:   fig_exvivo_exploratory_prefilter.R
# Chapter:  [1] — exvivo bulk rnaseq
# Step:     Figure script — loads RDS checkpoint from 02_exvivo_sample_exploratory.R
#           Preceded by: 02_exvivo_sample_exploratory.R
#           Followed by: n/a (terminal output script)
#
# Description:
#   Generates exploratory QC figures from the pre-filter snapshot of the
#   ex vivo sample set. These figures show the full sample distribution
#   before any inclusion criteria are applied and are used to justify the
#   filtering thresholds chosen in 02_exvivo_sample_exploratory.R.
#
#   Outputs:
#     Fig A — Mapping rate per study (threshold: none)
#     Fig B — Total mapped reads per study (threshold: 5,000,000)
#     Fig C — Gene detection rate per study (CPM > 5; threshold: 8,000 genes)
#     Fig D — Donor age per study (threshold: none)
#     Fig E — Control-tissue status per study (bar chart)
#
#   Inputs:
#     - rds/chapter_01/exvivo_sample_exploratory.rds
#                                 List: $counts (raw count matrix),
#                                       $metadata (sample-level metadata) —
#                                 pre-filter rds object from 02_exvivo_sample_exploratory.R
#
#   Outputs:
#     - outputs/chapter_01/figures/fig_exvivo_exploratory_prefilter.svg
#                                 Patchwork panel (A–E) of per-study QC boxplots
#                                 and bar chart; A4 portrait dimensions.
#
# Environment:
#   All scripts in this repository are intended to run inside a Docker
#   container. The image used for the RStudio session is available at:
#     https://hub.docker.com/r/stodo3569/rstudio-server_microglial-identity
#   See docker/ in this repository for the Dockerfile and build instructions.
#   Common utility functions (plot_qc_boxplot, qc_bar) are sourced from
#   common/R/00_functions.R (see repository root).
# ==============================================================================

source("microglial-identity/common/R/00_functions.R")

# ==============================================================================
# RDS management
# ==============================================================================
# This script loads one input RDS and produces no new RDS objects.
# SAVE_OUTPUT_RDS is always FALSE in figure scripts.
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

ZENODO_BASE <-  "https://zenodo.org/records/XXXXXXX/files/chapter_01/rds"
LOAD_INPUT_RDS  <- TRUE
SAVE_OUTPUT_RDS <- FALSE   # Figure scripts never save new RDS objects

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
library(patchwork)

# ==============================================================================
# Figure 2.01 — Per-study QC metric distributions (pre-filter)
# ==============================================================================
# Five panels showing the distribution of each sample-level QC metric across
# studies. Dashed threshold lines correspond to the inclusion criteria applied
# in 02_exvivo_sample_exploratory.R. All panels use exvivo_metadata, which at this
# checkpoint includes every sample before any filtering has been applied.
#
# plot_qc_boxplot() and qc_bar() are defined in common/R/00_functions.R.
# ==============================================================================

p1 <- plot_qc_boxplot(
  data_frame = exvivo_metadata,
  metric     = "mapping_rate",
  group_by   = "study",
  threshold  = 0
) +
  xlab("") +
  ylab("Mapping rate (%)")

p2 <- plot_qc_boxplot(
  data_frame      = exvivo_metadata,
  metric    = "total_reads",
  group_by  = "study",
  threshold = 5000000
) +
  xlab("") +
  ylab("Mapped reads")

p3 <- plot_qc_boxplot(
  data_frame      = exvivo_metadata,
  metric    = "sample_gene_count",
  group_by  = "study",
  threshold = 8000
) +
  xlab("") +
  ylab("Gene count per sample (CPM > 5)")

p4 <- plot_qc_boxplot(
  data_frame      = exvivo_metadata,
  metric    = "age",
  group_by  = "study",
  threshold = 0
) +
  xlab("") +
  ylab("Age (years)")

p5 <- qc_bar(
  data_frame = exvivo_metadata,
  metric     = "Control_tissue",
  group_by   = "study"
) +
  xlab("") +
  ylab("No. of samples")

fig <- (p1 | p2) / (p3 | p4) / (p5 | plot_spacer()) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14))

dir.create("outputs/chapter_01/figures", recursive = TRUE, showWarnings = FALSE)

svg(
  "outputs/chapter_01/figures/fig_exvivo_exploratory_prefilter.svg",
  width  = 8.27,
  height = 11.69
)
print(fig)
dev.off()

message("Saved: outputs/chapter_01/figures/fig_exvivo_exploratory_prefilter.svg")
