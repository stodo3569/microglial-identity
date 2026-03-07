# ==============================================================================
# Script:   [FILENAME].R
# Chapter:  [CHAPTER NUMBER] — [CHAPTER TITLE]
# Step:     [X] of [Y] in the [CHAPTER NUMBER] workflow
#           Preceded by: [PRECEDING SCRIPT NAME OR "none — first step"]
#           Followed by: [FOLLOWING SCRIPT NAME]
#
# Description:
#   Aggregates Salmon transcript-level quantification outputs (quant.sf files)
#   from all RNA-seq datasets used in the microglial identity analysis and
#   imports them into a single gene-level count matrix using tximport with the
#   'lengthScaledTPM' method. This produces count-scale values corrected for
#   per-sample transcript-length differences, making them suitable for DESeq2
#   input. All datasets were quantified against GENCODE v49 (GRCh38), except
#   Gosselin_2017, which used an older GENCODE release; this dataset is
#   imported separately and merged at the gene level (zero-filling genes absent
#   from either set) to produce a unified matrix across all samples. Basic
#   sanity checks on library sizes and gene detection rates are performed.
#   The script expects quant.sf files to be organised as:
#     data/[DATASET_NAME]/Aligned_data/[SAMPLE_ID]/quant.sf
#   where [SAMPLE_ID] is a GSM or SRX accession used as the column name in
#   the output matrix. The working directory is set to /home/rstudio, so
#   the data/ directory should be present there inside the container.
#
# Inputs:
#   - data/**/Aligned_data/[SAMPLE_ID]/quant.sf
#                                 Salmon quantification output, one per sample
#   - gencode.v49.annotation.gtf.gz
#                                 Downloaded automatically from GENCODE FTP if
#                                 not already present at /home/rstudio/
#
# Outputs:
#   - txi_all_datasets.rds        Merged tximport object containing $counts
#                                 (lengthScaledTPM-derived), $abundance (TPM),
#                                 and $length (effective transcript lengths)
#                                 for all samples; used by downstream scripts
#
# Environment:
#   All scripts in this repository are intended to run inside a Docker
#   container. The image used for the RStudio session is available at:
#     https://hub.docker.com/r/stodo3569/rstudio-server_microglial-identity
#   See docker/ in this repository for the Dockerfile and build instructions.
#   Common utility functions shared across scripts are sourced from
#   common/R/00_functions.R (see repository root).
# ==============================================================================

# Common utility functions shared across all scripts in this repository.
# The path is relative to this script's location within the local clone of the
# repository — R resolves it on disk, not from GitHub. It will work as long as
# you run this script from inside a cloned copy of the repository with its
# folder structure intact. The file can also be browsed online at:
# https://github.com/[YOUR-GITHUB-USERNAME]/[YOUR-REPO-NAME]/blob/main/common/R/00_functions.R
source("../../common/R/00_functions.R")

# ==============================================================================
# RDS management
# ==============================================================================
# This script reads raw Salmon quant.sf files and produces one RDS object:
#   txi_all_datasets.rds
#
# Three ways to proceed — set the flags below accordingly:
#
#   1) Skip this script entirely (DOWNLOAD_OUTPUT_RDS = TRUE):
#      Downloads the output RDS from Zenodo and stops. No quant.sf files
#      needed. Use this to jump straight to the next script.
#
#   2) Run this script from scratch (DOWNLOAD_OUTPUT_RDS = FALSE, default):
#      Requires raw Salmon quant.sf files on disk, organised as described in
#      the header above. These are available from GEO: [INSERT GEO ACCESSION].
#      Note: this script has no input RDS objects — its inputs are raw files,
#      so LOAD_INPUT_RDS does not apply here (see note below).
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE <- "https://zenodo.org/records/XXXXXXX/files"

# Input RDS objects required by this script:
#   None — this is the first script in the pipeline. Inputs are raw quant.sf
#   files read directly from disk (see step 2 below).
#
# LOAD_INPUT_RDS is included here for interface consistency with downstream
# scripts but has no effect in this script.
LOAD_INPUT_RDS <- TRUE

# Output RDS produced by this script (also available on Zenodo)
output_rds <- list(
  txi = "txi_all_datasets.rds"   # Merged tximport object: $counts, $abundance (TPM), $length
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
  message("Output RDS downloaded. Load with readRDS() and proceed to the next script.")
  stop("Stopping early — set DOWNLOAD_OUTPUT_RDS <- FALSE to run the full script.", call. = FALSE)
}

# SAVE_OUTPUT_RDS: set to TRUE (default) to save the output RDS to disk.
# Set to FALSE to skip saving — useful when exploring interactively.
SAVE_OUTPUT_RDS <- TRUE

# ==============================================================================
# Import Salmon quant.sf files into a unified gene expression matrix
# ==============================================================================
# Uses tximport with lengthScaledTPM to produce count-scale values corrected
# for per-sample transcript length differences — suitable for DESeq2 input.
#
# All datasets quantified with Salmon against GENCODE v49 (GRCh38), except
# Gosselin_2017 which used an older GENCODE release. tximport handles the
# transcript set differences transparently via ignoreTxVersion = TRUE.
# ==============================================================================

# ── 1. Load packages ───────────────────────────────────────────────
library(tximport)
library(readr)

# ── 2. Discover all quant.sf files ───────────────────────────────────────────
setwd("/home/rstudio")
base_dir <- "data"

all_quant_files <- list.files(
  path       = base_dir,
  pattern    = "quant\\.sf$",
  recursive  = TRUE,
  full.names = TRUE
)

# Keep only those inside Aligned_data directories
all_quant_files <- all_quant_files[grepl("/Aligned_data/", all_quant_files)]

# Name each file by its sample ID (GSM or SRX)
names(all_quant_files) <- basename(dirname(all_quant_files))

cat("Found", length(all_quant_files), "quant files across",
    length(unique(dirname(dirname(dirname(all_quant_files))))), "datasets\n")

# ── 3. Build tx2gene mapping from GENCODE v49 GTF ───────────────────────────
gtf_url  <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
gtf_file <- "/home/rstudio/gencode.v49.annotation.gtf.gz"

if (!file.exists(gtf_file)) {
  download.file(gtf_url, gtf_file)
}

# Parse transcript_id → gene_id from GTF (version suffixes stripped)
gtf_lines <- read_lines(gtf_file)
gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
tx_lines  <- gtf_lines[grepl("\ttranscript\t", gtf_lines)]

extract_attr <- function(lines, attr_name) {
  pattern <- paste0(attr_name, ' "([^"]+)"')
  m <- regmatches(lines, regexpr(pattern, lines))
  gsub(paste0(attr_name, ' "|"'), "", m)
}

tx_id   <- extract_attr(tx_lines, "transcript_id")
gene_id <- extract_attr(tx_lines, "gene_id")

# Strip version suffixes to match ignoreTxVersion = TRUE
tx_id   <- sub("\\.\\d+$", "", tx_id)
gene_id <- sub("\\.\\d+$", "", gene_id)

tx2gene <- data.frame(TXNAME = tx_id, GENEID = gene_id, stringsAsFactors = FALSE)
tx2gene <- unique(tx2gene)

cat("tx2gene mapping:", nrow(tx2gene), "transcripts →",
    length(unique(tx2gene$GENEID)), "genes\n")

# Free memory — GTF lines can be large
rm(gtf_lines, tx_lines, tx_id, gene_id)
gc()

# ── 4. Import with tximport ──────────────────────────────────────────────────
# Split files
gosselin_files <- all_quant_files[grepl("Gosselin_2017", all_quant_files)]
other_files    <- all_quant_files[!grepl("Gosselin_2017", all_quant_files)]

cat("Non-Gosselin:", length(other_files), "files\n")
cat("Gosselin:", length(gosselin_files), "files\n")

# Import separately
txi_other <- tximport(
  files               = other_files,
  type                = "salmon",
  tx2gene             = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion     = TRUE
)

txi_gosselin <- tximport(
  files               = gosselin_files,
  type                = "salmon",
  tx2gene             = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion     = TRUE
)

# Merge at gene level — union of all genes, zero-fill where absent
all_genes <- union(rownames(txi_other$counts), rownames(txi_gosselin$counts))

merge_matrices <- function(mat1, mat2, all_rows) {
  out <- matrix(0, nrow = length(all_rows),
                ncol = ncol(mat1) + ncol(mat2),
                dimnames = list(all_rows, c(colnames(mat1), colnames(mat2))))
  out[rownames(mat1), colnames(mat1)] <- mat1
  out[rownames(mat2), colnames(mat2)] <- mat2
  out
}

# Build merged txi object
txi <- list(
  counts           = merge_matrices(txi_other$counts, txi_gosselin$counts, all_genes),
  abundance        = merge_matrices(txi_other$abundance, txi_gosselin$abundance, all_genes),
  length           = merge_matrices(txi_other$length, txi_gosselin$length, all_genes),
  countsFromAbundance = "lengthScaledTPM"
)

cat("Merged gene expression matrix:", nrow(txi$counts), "genes ×",
    ncol(txi$counts), "samples\n")
cat("Genes in both:", length(intersect(rownames(txi_other$counts), rownames(txi_gosselin$counts))), "\n")
cat("Genes only in main:", length(setdiff(rownames(txi_other$counts), rownames(txi_gosselin$counts))), "\n")
cat("Genes only in Gosselin:", length(setdiff(rownames(txi_gosselin$counts), rownames(txi_other$counts))), "\n")

# ── 5. Quick sanity checks ──────────────────────────────────────────────────
# Library sizes per sample
lib_sizes <- colSums(txi$counts)
cat("\nLibrary size range:",
    round(min(lib_sizes) / 1e6, 1), "M –",
    round(max(lib_sizes) / 1e6, 1), "M\n")

# Genes with non-zero expression
genes_detected <- rowSums(txi$counts > 0)
cat("Genes detected in at least 1 sample:", sum(genes_detected >= 1), "\n")
cat("Genes detected in all samples:", sum(genes_detected == ncol(txi$counts)), "\n")

# ── 6. Save outputs ─────────────────────────────────────────────────────────
if (SAVE_OUTPUT_RDS) {
  saveRDS(txi, file = "/home/rstudio/txi_all_datasets.rds")
  cat("\nSaved txi object to /home/rstudio/txi_all_datasets.rds\n")
  cat("Contents: txi$counts, txi$abundance (TPM), txi$length\n")
}

exvivo_counts <- txi$counts
