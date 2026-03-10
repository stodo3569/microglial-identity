# ==============================================================================
# Script:   fig_exvivo_sample_qc.R
# Chapter:  [1] — exvivo bulk rnaseq
# Step:     Figure script — loads RDS checkpoint from 03_exvivo_sample_qc.R
#           Preceded by: 03_exvivo_sample_qc.R
#           Followed by: n/a (terminal output script)
#
# Description:
#   Generates multi-study QC figures from the post-QC snapshot of the
#   ex vivo sample set. The three panels visualise each of the simultaneous
#   outlier-removal criteria applied by custom_DaMiR.sampleFilt() in
#   03_exvivo_sample_qc.R. Because all metrics are computed on the full
#   sample set before any removal, every sample has a value in every panel.
#
#   Outputs:
#     Fig A — Inter-study connectivity Z-score per study (threshold: -1)
#              All samples have a value; no exemptions apply.
#
#     Fig B — Within-study connectivity Z-score per study (threshold: -1)
#              Studies with <= 3 samples cannot produce a stable Z-score and
#              are shown as open triangles at y = 0 to indicate exemption.
#
#     Fig C — WGCNA standardised connectivity (Z.k) per study (threshold: -1)
#              All samples have a value; no placeholder markers needed.
#
#   In all panels: green points = retained; black points = removed by any
#   filter. Open triangles (panel B only) = exempt from within-study filter.
#
#   Inputs:
#     - rds/chapter_01/exvivo_sample_qc.rds
#                                 List: $counts, $metadata, $SE_norm,
#                                       $qc_summary —
#                                 post-QC checkpoint from 03_exvivo_sample_qc.R
#
#   Outputs:
#     - outputs/chapter_01/figures/fig_exvivo_sample_qc.svg
#                                 Patchwork panel (A–C) of per-study
#                                 connectivity and WGCNA Z.k dot plots;
#                                 A4 landscape dimensions.
#
# Environment:
#   All scripts in this repository are intended to run inside a Docker
#   container. The image used for the RStudio session is available at:
#     https://hub.docker.com/r/stodo3569/rstudio-server_microglial-identity
#   See docker/ in this repository for the Dockerfile and build instructions.
#   Common utility functions (plot_qc_boxplot) are sourced from
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
SAVE_OUTPUT_RDS <- FALSE   # Figure scripts never save new RDS objects

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
library(patchwork)

# ==============================================================================
# Figure 2.02 — Multi-study QC: per-study outlier detection
# ==============================================================================
# Three panels corresponding to the simultaneous outlier-removal criteria in
# custom_DaMiR.sampleFilt(). All metrics were computed on the full sample set
# before any removal, so every sample appears in every panel.
#
# Colouring convention (applied by plot_qc_boxplot):
#   Green filled point  — sample retained (pass_all_filters == TRUE)
#   Black filled point  — sample removed  (pass_all_filters == FALSE)
#   Open triangle       — panel B only: study exempt from within-study Z
#                         (study_n <= 3; within_study_z is NA)
#
# Label formatting:
#   Database accession suffixes (GSE/SRP/EGAD followed by digits) are stripped
#   once here before any data reaches plot_qc_boxplot, leaving study names in
#   the form Author_Year (e.g. Kana_2019). Since plot_qc_boxplot no longer
#   modifies labels internally, boxplot layers and all external overlays
#   (triangle markers) use identical labels and align correctly on the x-axis.
#
# Panel B:
#   sort_by_mean = FALSE preserves the factor level ordering defined here,
#   ensuring all studies — including exempt ones with entirely NA
#   within_study_z — appear on the x-axis. NA values are left as-is so
#   geom_jitter produces no point for exempt studies; only the open triangle
#   overlay marks them at y = 0.
#
# plot_qc_boxplot() is defined in common/R/00_functions.R.
# ==============================================================================

# Strip database accession suffix (e.g. _GSE133357, _SRP115940, _EGAD00001005736)
# leaving Author_Year as the display label.
strip_accession <- function(x) sub("_(GSE|SRP|EGAD)\\d+$", "", x)

study_levels <- unique(exvivo_qc_summary$study)

qc_plot_df <- mutate(exvivo_qc_summary,
                     study = strip_accession(study),
                     study = factor(study,
                                    levels = unique(strip_accession(study_levels))))

# Panel B data: flag exempt studies but keep within_study_z as NA so that
# geom_jitter produces no point for them — only the triangle overlay appears.
qc_plot_B <- mutate(qc_plot_df, within_study_exempt = is.na(within_study_z))

# Panel A — inter-study connectivity Z-score
pA <- plot_qc_boxplot(
  data_frame = qc_plot_df,
  metric     = "inter_study_z",
  group_by   = "study",
  threshold  = -1,
  sort_by_mean = FALSE
) +
  xlab("") +
  ylab("Inter-study connectivity\n(Z-score)")

# Panel B — within-study connectivity Z-score
pB <- plot_qc_boxplot(
  data_frame   = qc_plot_B,
  metric       = "within_study_z",
  group_by     = "study",
  threshold    = -1,
  sort_by_mean = FALSE
) +
  geom_point(
    data  = filter(qc_plot_B, within_study_exempt == TRUE),
    aes(x = study, y = 0),
    shape = 2, size = 2.5, colour = "grey50"
  ) +
  xlab("") +
  ylab("Within-study connectivity\n(Z-score)")

# Panel C — WGCNA standardised connectivity (Z.k)
qc_plot_C <- mutate(qc_plot_df, zk_exempt = is.na(Z.k))

zk_exempt_studies <- qc_plot_df %>%
  group_by(study) %>%
  filter(all(is.na(Z.k))) %>%
  distinct(study) %>%
  mutate(y = 0)

pC <- plot_qc_boxplot(
  data_frame   = qc_plot_df,
  metric       = "Z.k",
  group_by     = "study",
  threshold    = -1,
  sort_by_mean = FALSE
) +
  geom_point(
    data  = zk_exempt_studies,
    aes(x = study, y = y),
    shape = 4, size = 2.5, colour = "grey50"
  ) +
  xlab("") +
  ylab("WGCNA standardised\nconnectivity (Z.k)")

fig <- (pA / pB / pC) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14))

dir.create("outputs/chapter_01/figures", recursive = TRUE, showWarnings = FALSE)

svg(
  "outputs/chapter_01/figures/fig_exvivo_sample_qc.svg",
  width  = 11.69,
  height = 8.27
)
print(fig)
dev.off()

message("Saved: outputs/chapter_01/figures/fig_exvivo_sample_qc.svg")
