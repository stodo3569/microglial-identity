# ==============================================================================
# Script:   tab_geneset_database.R
# Chapter:  [common] — shared across all chapters
# Step:     Table script — loads RDS from 00_geneset_database.R
#           Preceded by: 00_geneset_database.R
#           Followed by: n/a (terminal output script)
#
# Description:
#   Generates a summary table of the compiled geneset database, showing for
#   each gene set library: the date it was accessed, the number of unique
#   gene set terms, and the average number of genes per term. A summary row
#   gives column sums and means across all libraries.
#
#   Inputs:
#     - rds/common/geneset_database.rds
#                                 Named list of data.frames, one per source.
#                                 Each element: data.frame(Term, ensembl_gene_id)
#                                 Names follow the pattern "{SourceName}_{YYYY-MM-DD}"
#
#   Outputs:
#     - outputs/common/tables/tab_geneset_database.html
#                                 gt table summarising per-source term counts
#                                 and mean genes per term; exported as a
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
#     rds/common/geneset_database.rds     Named list, one data.frame(Term,
#                                         ensembl_gene_id) per source.
#
# Two ways to proceed — set the flag below accordingly:
#
#   1) Load from Zenodo (LOAD_INPUT_RDS = TRUE, default):
#      Downloads geneset_database.rds from Zenodo if not present locally,
#      then runs this script in full.
#
#   2) Run fully from local files (LOAD_INPUT_RDS = FALSE):
#      Expects geneset_database.rds to be present locally at
#      rds/common/geneset_database.rds (produced by running
#      00_geneset_database.R).
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE     <- "https://zenodo.org/records/XXXXXXX/files/common/rds"
LOAD_INPUT_RDS  <- TRUE
SAVE_OUTPUT_RDS <- FALSE   # Table scripts never save new RDS objects

input_rds <- list(
  geneset_database = "rds/common/geneset_database.rds"
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
  geneset_db <- readRDS(input_rds$geneset_database)
}

# ==============================================================================
# Set up working environment
# ==============================================================================

library(tidyverse)
library(gt)

# ==============================================================================
# Table — Geneset database summary
# ==============================================================================
# One row per gene set library. Source name and access date are parsed from
# the list key (format: "{SourceName}_{YYYY-MM-DD}"). Unique terms are counted
# with n_distinct(Term); mean genes per term is the mean of per-term
# ensembl_gene_id counts. A grand-total row is appended and rendered in bold:
# summing unique terms across libraries and averaging mean genes per library.
# ==============================================================================

tab_rows <- lapply(seq_along(geneset_db), function(i) {
  key         <- names(geneset_db)[i]
  df          <- geneset_db[[i]]
  terms       <- as.character(unlist(df$Term, use.names = FALSE))
  date_match  <- regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}$", key)
  date_str    <- if (date_match > 0) regmatches(key, date_match) else NA_character_
  source_name <- sub("_[0-9]{4}-[0-9]{2}-[0-9]{2}$", "", key)
  counts      <- as.integer(table(terms))   # one count per unique term

  data.frame(
    source        = source_name,
    date_accessed = date_str,
    n_terms       = length(counts),
    mean_genes    = if (length(counts) > 0) mean(counts) else NA_real_,
    stringsAsFactors = FALSE
  )
})

tab_data <- do.call(rbind, tab_rows)

tab_data <- rbind(
  tab_data,
  data.frame(
    source        = "Total",
    date_accessed = NA_character_,
    n_terms       = sum(tab_data$n_terms),
    mean_genes    = mean(tab_data$mean_genes),
    stringsAsFactors = FALSE
  )
)

tab_gt <- tab_data %>%
  gt() %>%
  cols_label(
    source        = "Gene set library",
    date_accessed = "Date accessed",
    n_terms       = "Unique terms",
    mean_genes    = "Mean genes per term"
  ) %>%
  fmt_number(
    columns  = n_terms,
    decimals = 0
  ) %>%
  fmt_number(
    columns  = mean_genes,
    decimals = 1
  ) %>%
  sub_missing(
    columns      = date_accessed,
    missing_text = "—"
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_body(rows = nrow(tab_data))
  )

# ==============================================================================
# Save output
# ==============================================================================

dir.create(
  "outputs/common/tables",
  recursive    = TRUE,
  showWarnings = FALSE
)

gtsave(
  tab_gt,
  filename = "outputs/common/tables/tab_geneset_database.html"
)

message("Saved: outputs/common/tables/tab_geneset_database.html")
