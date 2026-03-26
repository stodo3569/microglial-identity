# ==============================================================================
# Script:   00_geneset_database.R
# Chapter:  [common] — shared across all chapters
# Step:     Standalone — run once to build the geneset database
#           Preceded by: n/a
#           Followed by: any script performing enrichment analysis
#
# Description:
#   Downloads 13+ curated pathway and gene-set databases (Reactome, KEGG,
#   WikiPathways, GO, MSigDB collections, PathBank, Jensen Compartments,
#   Panther, NCI-PID, ACSN, BioCarta, HumanCyc, NetPath), maps all gene
#   identifiers to current HGNC symbols and then to Ensembl gene IDs, and
#   compiles the results into a single named-list RDS for use by all
#   downstream enrichment analysis scripts.
#
#   Gene identifier handling follows the conventions in CLAUDE.md:
#   update_gene_symbols() first, then update_gene_to_ensembl().
#   HGNC-ID sources use update_hgnc_to_symbols() before that chain.
#
# Inputs:
#   - HGNC complete set (downloaded automatically from storage.googleapis.com)
#   - Per-source pathway/GMT files (downloaded automatically; see each
#     process_*() function for the specific URL)
#
# Outputs:
#   - rds/common/geneset_database.rds
#                                 Named list of data.frames, one per source.
#                                 Each element: data.frame(Term, ensembl_gene_id)
#                                 Names: "Reactome", "KEGG", "WikiPathways",
#                                 "GO-Biological-Processes", etc.
#                                 Loaded by enrichment analysis scripts via
#                                 readRDS("rds/common/geneset_database.rds")
#
# Downstream usage:
#   geneset_db <- readRDS("rds/common/geneset_database.rds")
#
#   # Single source:
#   custom_enricher(gene = my_genes, universe = background,
#                   geneSets = geneset_db[["Reactome"]])
#
#   # Combined sources:
#   combined <- do.call(rbind, geneset_db[c("Reactome", "KEGG", "WikiPathways")])
#   custom_enricher(gene = my_genes, universe = background, geneSets = combined)
#
# Environment:
#   All scripts in this repository are intended to run inside a Docker
#   container. The image used for the RStudio session is available at:
#     https://hub.docker.com/r/stodo3569/rstudio-server_microglial-identity
#   See docker/ in this repository for the Dockerfile and build instructions.
#   Common utility functions shared across scripts are sourced from
#   common/R/00_functions.R (see repository root).
# ==============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(readr)
library(data.table)
library(httr)
library(jsonlite)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(KEGGREST)
library(rWikiPathways)
library(foreach)
library(doParallel)
library(iterators)
library(rvest)
library(GO.db)
library(DBI)
library(RSQLite)
library(xml2)
library(ndexr)
library(tidyr)

# ==============================================================================
# RDS management
# ==============================================================================
# This script produces one output RDS:
#   rds/common/geneset_database.rds
#
# Two ways to proceed — set the flags below accordingly:
#
#   1) Skip this script entirely (DOWNLOAD_OUTPUT_RDS = TRUE):
#      Downloads the output RDS from Zenodo and stops. No external database
#      connections needed. Use this to jump straight to enrichment scripts.
#
#   2) Run this script from scratch (DOWNLOAD_OUTPUT_RDS = FALSE, default):
#      Downloads all pathway databases from their original sources, processes
#      gene identifiers, and compiles the RDS. Requires internet access.
#
# Replace XXXXXXX in ZENODO_BASE with the actual Zenodo record ID once the
# dataset has been deposited.
# ==============================================================================

ZENODO_BASE <- "https://zenodo.org/records/XXXXXXX/files/common/rds"

output_rds <- list(
  geneset_database = "rds/common/geneset_database.rds"
)

DOWNLOAD_OUTPUT_RDS <- FALSE

if (DOWNLOAD_OUTPUT_RDS) {
  dir.create("rds/common", recursive = TRUE, showWarnings = FALSE)
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
  message("Output RDS downloaded. Load with readRDS() and proceed to enrichment scripts.")
  stop("Stopping early — set DOWNLOAD_OUTPUT_RDS <- FALSE to run the full script.", call. = FALSE)
}

SAVE_OUTPUT_RDS <- TRUE

# Set to TRUE to also write per-source TSV files to term_lib/ (for inspection/archiving).
# Default FALSE — the RDS is the primary output.
SAVE_TSVS <- FALSE

term_dir <- "term_lib/"

# ==============================================================================
# Reference data
# ==============================================================================

# HGNC database; general
hgnc_database <- fread("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt")

# Subset databases for converting ENTREZ and ENSEMBL ids to symbols
symbol_entrez_db  <- hgnc_database %>% dplyr::select(symbol, entrez_id)
symbol_ensembl_db <- hgnc_database %>% dplyr::select(symbol, ensembl_gene_id)
hgnc_ensembl_db   <- hgnc_database %>% dplyr::select(hgnc_id, ensembl_gene_id) %>% mutate(hgnc_id = gsub("HGNC:", "", hgnc_id))
hgnc_symbol_db    <- hgnc_database %>% dplyr::select(hgnc_id, symbol) %>% mutate(hgnc_id = gsub("HGNC:", "", hgnc_id))

# HGNC database containing current and previous symbols
genes <- hgnc_database %>%
  dplyr::select(symbol, prev_symbol) %>%
  separate_rows(prev_symbol, sep = "\\|") %>%
  pivot_longer(cols = -symbol, names_to = "origin", values_to = "alt_symbol") %>%
  dplyr::select(symbol, alt_symbol) %>%
  as.data.table()

# ==============================================================================
# Gene identifier helpers
# ==============================================================================

# Function to calculate Jaccard similarity index between two sets
jaccard_similarity <- function(set1, set2) {
  length(intersect(set1, set2)) / length(union(set1, set2))
}

# Function to update gene symbols using the provided genes data.table
update_gene_symbols <- function(lib_data, genes) {
  tmp <- tempfile(fileext = ".sqlite")
  con <- dbConnect(RSQLite::SQLite(), dbname = tmp)
  on.exit({
    dbDisconnect(con)
    file.remove(tmp)
  })

  dbWriteTable(con, "lib_data", lib_data)
  dbWriteTable(con, "genes", genes)

  query <- '
  SELECT lib_data.Term, lib_data.Gene, genes.symbol
  FROM lib_data
  LEFT JOIN genes
  ON lib_data.Gene = genes.alt_symbol
  '
  lib_data_joined <- dbGetQuery(con, query) %>% as.data.table()
  lib_data_joined[, symbol := ifelse(Gene %in% genes$symbol, Gene, symbol)]
  lib_data_updated <- lib_data_joined[!is.na(symbol), .(Term, Gene = symbol)]
  distinct(as_tibble(lib_data_updated))
}

# Function to update gene symbols to Ensembl IDs
update_gene_to_ensembl <- function(lib_data, symbol_ensembl_db) {
  tmp <- tempfile(fileext = ".sqlite")
  con <- dbConnect(RSQLite::SQLite(), dbname = tmp)
  on.exit({
    dbDisconnect(con)
    file.remove(tmp)
  })

  dbWriteTable(con, "lib_data", lib_data)
  dbWriteTable(con, "symbol_ensembl_db", symbol_ensembl_db)

  query <- '
  SELECT lib_data.Term, lib_data.Gene, symbol_ensembl_db.ensembl_gene_id
  FROM lib_data
  LEFT JOIN symbol_ensembl_db
  ON lib_data.Gene = symbol_ensembl_db.symbol
  '
  lib_data_joined <- dbGetQuery(con, query) %>% as.data.table()
  lib_data_updated <- lib_data_joined[!is.na(ensembl_gene_id) & ensembl_gene_id != "", .(Term, ensembl_gene_id)]
  distinct(as_tibble(lib_data_updated))
}

# Function to update hgnc ids to gene symbols
update_hgnc_to_symbols <- function(lib_data, hgnc_symbol_db) {
  tmp <- tempfile(fileext = ".sqlite")
  con <- dbConnect(RSQLite::SQLite(), dbname = tmp)
  on.exit({
    dbDisconnect(con)
    file.remove(tmp)
  })

  dbWriteTable(con, "lib_data", lib_data)
  dbWriteTable(con, "hgnc_symbol_db", hgnc_symbol_db)

  query <- '
  SELECT lib_data.Term, hgnc_symbol_db.symbol AS Gene
  FROM lib_data
  LEFT JOIN hgnc_symbol_db
  ON lib_data.Gene = hgnc_symbol_db.hgnc_id
  '

  lib_data_joined <- dbGetQuery(con, query) %>% as.data.table()
  distinct(as_tibble(lib_data_joined[!is.na(Gene)]))
}

# ==============================================================================
# Source-specific processing functions
# ==============================================================================

# process_pathbank() -----------------------------------------------------------
# Downloads PathBank all-proteins CSV (zipped), filters to Homo sapiens,
# resolves gene symbols, maps to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_pathbank <- function(pathbank_url, term_dir, genes, symbol_ensembl_db,
                             save_tsv = FALSE) {
  if (!httr::http_error(pathbank_url)) {
    temp_zip <- tempfile()
    download.file(pathbank_url, temp_zip, mode = "wb")
    pathbank_data <- read_csv(unz(temp_zip, "pathbank_all_proteins.csv"),
                              col_types = cols(.default = col_character()),
                              show_col_types = FALSE)
    unlink(temp_zip)

    pathbank_processed <- pathbank_data %>%
      filter(Species == "Homo sapiens") %>%
      mutate(Term = `Pathway Name`, Gene = `Gene Name`) %>%
      dplyr::select(Term, Gene) %>%
      filter(!is.na(Gene) & Gene != "" & !grepl("^[0-9]", Gene))

    pathbank_processed <- update_gene_symbols(pathbank_processed, genes)
    pathbank_processed <- update_gene_to_ensembl(pathbank_processed, symbol_ensembl_db)

    if (save_tsv) {
      write_tsv(pathbank_processed, paste0(term_dir, "PathBank_", Sys.Date(), ".tsv"))
    }

    return(pathbank_processed)
  } else {
    cat("Error: Invalid PathBank source URL. Please check the URL and try again.")
    return(invisible(NULL))
  }
}

# process_jensen_compartment() -------------------------------------------------
# Downloads Jensen Lab compartment associations TSV (gene symbol, GO term, score).
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_jensen_compartment <- function(jencom_url, term_dir, genes, symbol_ensembl_db,
                                       save_tsv = FALSE) {
  if (!httr::http_error(jencom_url)) {
    jencom_data <- fread(jencom_url, header = FALSE)

    # Columns: gene_symbol, gene_alias, GO_id, term_name, score
    colnames(jencom_data) <- c("Gene", "Gene_alias", "GO_id", "Term", "Score")

    jencom_processed <- jencom_data %>%
      dplyr::select(Term, Gene)

    jencom_processed <- update_gene_symbols(jencom_processed, genes)
    jencom_processed <- update_gene_to_ensembl(jencom_processed, symbol_ensembl_db)

    if (save_tsv) {
      write_tsv(jencom_processed, paste0(term_dir, "Jensen-Compartments_", Sys.Date(), ".tsv"))
    }

    return(jencom_processed)
  } else {
    cat("Error: Invalid Jensen-Compartment source URL. Please check the URL and try again.")
    return(invisible(NULL))
  }
}

# process_go() -----------------------------------------------------------------
# Downloads the human GOA annotation file (GAF format, gzip) and splits by
# GO aspect into three data.frames.
#
# Returns: named list with elements "GO-Biological-Processes",
#          "GO-Molecular-Function", "GO-Cellular-Component",
#          each a data.frame(Term, ensembl_gene_id); or invisible(NULL) on error
process_go <- function(url, term_dir, genes, symbol_ensembl_db, save_tsv = FALSE) {

  if (!httr::http_error(url)) {

    # Read the GO data, skipping comment lines starting with '!'
    GO_db <- data.table::fread(
      cmd = paste("curl -sL", url, "| zcat | grep -v '^!'"),
      header = FALSE
    )

    # GAF columns: V1=DB, V2=DB_Object_ID, V3=Gene_Symbol, V4=Qualifier,
    #              V5=GO_ID, V6=DB_Ref, V7=Evidence, V8=With, V9=Aspect,
    #              V10=Name, V11=Synonym, V12=Type, V13=Taxon, V14=Date, V15=Assigned_By

    go_term_names <- AnnotationDbi::Term(GOTERM)

    aspects <- list(
      P = "GO-Biological-Processes",
      F = "GO-Molecular-Function",
      C = "GO-Cellular-Component"
    )

    result <- list()

    for (code in names(aspects)) {
      aspect_name <- aspects[[code]]

      go_subset <- GO_db %>%
        filter(V9 == code & !grepl("NOT", V4)) %>%
        dplyr::select(Term = V5, Gene = V3) %>%
        mutate(Term = go_term_names[Term]) %>%
        filter(!is.na(Term) & Term != "")

      go_subset <- update_gene_symbols(go_subset, genes)
      go_subset <- update_gene_to_ensembl(go_subset, symbol_ensembl_db)

      if (save_tsv) {
        write_tsv(go_subset, paste0(term_dir, "Gene-Ontology_", aspect_name, "_", Sys.Date(), ".tsv"))
      }

      result[[aspect_name]] <- go_subset
    }

    return(result)

  } else {
    cat("Error: Invalid GO source URL. Please check the URL and try again.")
    return(invisible(NULL))
  }
}

# process_kegg() ---------------------------------------------------------------
# Fetches KEGG pathway-gene links via KEGGREST, converts Entrez IDs to symbols,
# and maps to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_kegg <- function(term_dir, genes, symbol_entrez_db, symbol_ensembl_db,
                         save_tsv = FALSE) {

  res <- keggLink("hsa", "pathway")

  kegg_lib <- data.frame(Term = names(res), Gene = unname(res)) %>%
    mutate(Term = gsub("path:", "", Term),
           Gene = gsub("hsa:", "", Gene))

  # Replace pathway IDs with human-readable names
  all_pathways <- keggList("pathway", "hsa")
  pathway_id_name_mapping <- data.frame(
    Term = names(all_pathways),
    Name = all_pathways,
    stringsAsFactors = FALSE
  )
  kegg_lib <- kegg_lib %>%
    inner_join(pathway_id_name_mapping, by = "Term") %>%
    mutate(Term = Name) %>%
    dplyr::select(-Name)

  # Convert Entrez IDs to gene symbols
  kegg_lib <- kegg_lib %>%
    mutate(Gene = as.numeric(Gene)) %>%
    left_join(symbol_entrez_db, by = c("Gene" = "entrez_id")) %>%
    filter(!is.na(symbol) & symbol != "" & !grepl("^[0-9]", symbol)) %>%
    mutate(Gene = symbol) %>%
    dplyr::select(Term, Gene)

  kegg_lib <- update_gene_symbols(kegg_lib, genes)
  kegg_lib <- update_gene_to_ensembl(kegg_lib, symbol_ensembl_db)

  if (save_tsv) {
    write_tsv(kegg_lib, paste0(term_dir, "KEGG_", as.character(Sys.Date()), ".tsv"))
  }

  return(kegg_lib)
}

# process_wiki() ---------------------------------------------------------------
# Scrapes the WikiPathways GMT download page to find the latest Homo sapiens
# GMT, converts Entrez IDs to symbols, and maps to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_wiki <- function(url, term_dir, genes, symbol_entrez_db, symbol_ensembl_db,
                         save_tsv = FALSE) {

  page <- read_html(url)

  # Find all links on the page
  links <- html_nodes(page, "a") %>% html_attr("href")

  homo_sapiens_link <- links[grepl("Homo_sapiens.gmt", links)] %>% paste0("https://data.wikipathways.org/current/gmt/", .)

  # Read WikiPathways data
  hs_WikiPathways_list <- readLines(homo_sapiens_link) %>%
    lapply(function(x) strsplit(x, "\t", fixed = TRUE)[[1]])

  hs_WikiPathways <- do.call(rbind, lapply(hs_WikiPathways_list, function(x) {
    data.frame(Term = x[1], Gene = x[-(1:2)], stringsAsFactors = FALSE)
  }))

  # Convert Entrez IDs to gene symbols
  hs_WikiPathways_processed <- hs_WikiPathways %>%
    filter(Gene != "", !grepl("http", Gene)) %>%
    mutate(Term = sub("%WikiPathways.*", "", Term),
           Gene = as.numeric(Gene)) %>%
    left_join(symbol_entrez_db, by = c("Gene" = "entrez_id")) %>%
    filter(!is.na(symbol) & symbol != "" & !grepl("^[0-9]", symbol)) %>%
    mutate(Gene = symbol) %>%
    dplyr::select(Term, Gene)

  hs_WikiPathways_processed <- update_gene_symbols(hs_WikiPathways_processed, genes)
  hs_WikiPathways_processed <- update_gene_to_ensembl(hs_WikiPathways_processed, symbol_ensembl_db)

  if (save_tsv) {
    write_tsv(hs_WikiPathways_processed, paste0(term_dir, "WikiPathways_", as.character(Sys.Date()), ".tsv"))
  }

  return(hs_WikiPathways_processed)
}

# process_panther() ------------------------------------------------------------
# Scrapes the PANTHER FTP current release page to find the SequenceAssociation
# pathway file, extracts HGNC IDs, converts to symbols, and maps to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_panther <- function(url, term_dir, hgnc_symbol_db, genes, symbol_ensembl_db,
                            save_tsv = FALSE) {

  response <- GET(url)
  webpage <- rvest::read_html(response, encoding = "UTF-8")
  links <- webpage %>% html_nodes("a") %>% html_attr("href")

  matching_links <- grep("SequenceAssociationPathway.*\\.txt$", links, value = TRUE)
  full_url <- paste0(url, matching_links)

  if (!http_error(full_url)) {
    file_content <- readLines(full_url)

    lines_with_hgnc <- file_content[grepl("HGNC", file_content)]

    panther_db <- tibble(
      Term = str_extract(lines_with_hgnc, "(?<=\\t)[^\\t]+"),
      Gene = str_extract(lines_with_hgnc, "(?<=HGNC=)\\d+")
    ) %>%
      filter(!is.na(Term) & !is.na(Gene) & Gene != "")

    panther_db <- update_hgnc_to_symbols(panther_db, hgnc_symbol_db)
    panther_db <- update_gene_symbols(panther_db, genes)
    panther_db <- update_gene_to_ensembl(panther_db, symbol_ensembl_db)

    if (save_tsv) {
      write_tsv(panther_db, paste0(term_dir, "Panther_", as.character(Sys.Date()), ".tsv"))
    }

    return(panther_db)
  } else {
    cat("Error: Invalid source URL. Please check the URL and try again.\n")
    return(invisible(NULL))
  }
}

# process_nci() ----------------------------------------------------------------
# Fetches all NCI-PID networks from NDEx (accountName = "nci-pid"), extracts
# gene node names, and maps to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_nci <- function(term_dir, genes, symbol_ensembl_db, save_tsv = FALSE) {
  ndexcon <- ndex_connect()

  start <- 0
  size <- 600
  all_networks <- list()
  repeat {
    networks_chunk <- ndex_find_networks(ndexcon, start = start, size = size, accountName = "nci-pid")
    if (length(networks_chunk) == 0) break
    all_networks <- append(all_networks, list(as.data.frame(networks_chunk)))
    start <- start + size
  }

  common_columns <- Reduce(intersect, lapply(all_networks, colnames))
  all_networks_df <- do.call(rbind, lapply(all_networks, function(df) df[, common_columns]))

  extract_genes <- function(network_id) {
    ndex_get_network(ndexcon, network_id)$nodes$n
  }

  term_gene_list <- list()
  for (i in seq_len(nrow(all_networks_df))) {
    gene_names <- extract_genes(all_networks_df[i, "externalId"])
    term_gene_list <- append(term_gene_list, list(
      data.frame(Term = all_networks_df[i, "name"], Gene = gene_names, stringsAsFactors = FALSE)
    ))
  }

  term_gene_df <- do.call(rbind, term_gene_list) %>%
    filter(!is.na(Gene) & Gene != "" & !grepl("^[0-9]", Gene))

  term_gene_df <- update_gene_symbols(term_gene_df, genes)
  term_gene_df <- update_gene_to_ensembl(term_gene_df, symbol_ensembl_db)

  if (save_tsv) {
    write_tsv(term_gene_df, paste0(term_dir, "NCI-PID_v2_", as.character(Sys.Date()), ".tsv"))
  }

  return(term_gene_df)
}

# process_humancyc() -----------------------------------------------------------
# Downloads the combined Pathway Commons GMT (pc-hgnc.gmt.gz), filters to
# HumanCyc entries, extracts pathway names from the metadata column, and maps
# HGNC symbols to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_humancyc <- function(term_dir, genes, symbol_ensembl_db, save_tsv = FALSE) {
  index_url <- "https://download.baderlab.org/PathwayCommons/PC2/"

  if (!httr::http_error(index_url)) {

    # Version discovery: dirs are "v{N}/" — sort by integer after stripping "v"
    page         <- rvest::read_html(index_url)
    links        <- rvest::html_nodes(page, "a") %>% rvest::html_attr("href")
    version_dirs <- grep("^v[0-9]+/$", links, value = TRUE)

    if (length(version_dirs) == 0) {
      cat("Error: Could not identify Pathway Commons version from index.\n")
      return(invisible(NULL))
    }

    latest  <- version_dirs[which.max(as.integer(sub("^v", "", sub("/$", "", version_dirs))))]
    gmt_url <- paste0(index_url, latest, "pc-hgnc.gmt.gz")

    # Download combined GMT to tempfile
    temp_gz <- tempfile(fileext = ".gmt.gz")
    download.file(gmt_url, temp_gz, mode = "wb")

    if (!file.exists(temp_gz) || file.size(temp_gz) == 0) {
      cat("Error: Pathway Commons GMT download is empty.\n")
      return(invisible(NULL))
    }

    # Read gzip directly — no extraction to disk
    con   <- gzcon(file(temp_gz, open = "rb"))
    lines <- readLines(con)
    close(con)
    unlink(temp_gz)

    # Filter HumanCyc entries (col1 starts with "humancyc:")
    lines <- lines[grepl("^humancyc:", lines)]

    # Parse GMT: col1 = pathway ID (skip), col2 = metadata with name, col3+ = HGNC symbols
    # col2 format: "name: PATHWAY_NAME; datasource: humancyc; organism: 9606; idtype: hgnc.symbol"
    humancyc_df <- do.call(rbind, lapply(lines, function(x) {
      parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
      term  <- sub(".*name: ([^;]+);.*", "\\1", parts[2])
      data.frame(Term = term, Gene = parts[-(1:2)], stringsAsFactors = FALSE)
    }))

    humancyc_df <- humancyc_df %>%
      filter(!is.na(Gene) & Gene != "" & !grepl("^[0-9]", Gene))

    humancyc_df <- update_gene_symbols(humancyc_df, genes)
    humancyc_df <- update_gene_to_ensembl(humancyc_df, symbol_ensembl_db)

    if (save_tsv) {
      write_tsv(humancyc_df, paste0(term_dir, "HumanCyc_", Sys.Date(), ".tsv"))
    }

    return(humancyc_df)

  } else {
    cat("Error: Pathway Commons index unavailable. Please check the URL and try again.\n")
    return(invisible(NULL))
  }
}

# process_acsn() ---------------------------------------------------------------
# Scrapes the ACSN downloads page for *_master.gmt files (cell signalling
# networks), parses GMT format, and maps symbols to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_acsn <- function(url, term_dir, genes, symbol_ensembl_db, save_tsv = FALSE) {

  page_content <- GET(url, user_agent("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"))

  if (status_code(page_content) == 200) {
    webpage <- rvest::read_html(page_content, encoding = "UTF-8")

    gmt_links <- webpage %>%
      html_nodes("a") %>%
      html_attr("href") %>%
      grep("\\_master.gmt$", ., value = TRUE) %>%
      .[-1] %>%
      paste0("https://acsn.curie.fr/ACSN2/", .)

    acsn_data <- lapply(gmt_links, function(link) {
      gmt <- readLines(link) %>% lapply(function(x) strsplit(x, "\t", fixed = TRUE)[[1]])
      do.call(rbind, lapply(gmt, function(x) data.frame(Term = x[1], Gene = x[-1], stringsAsFactors = FALSE)))
    })

    acsn_processed <- do.call(rbind, acsn_data) %>%
      filter(!is.na(Gene) & Gene != "" & !grepl("^[0-9]", Gene))

    acsn_processed <- update_gene_symbols(acsn_processed, genes)
    acsn_processed <- update_gene_to_ensembl(acsn_processed, symbol_ensembl_db)

    if (save_tsv) {
      write_tsv(acsn_processed, paste0(term_dir, "ACSN_v2_", as.character(Sys.Date()), ".tsv"))
    }

    return(acsn_processed)
  } else {
    cat("ACSN server unavailable. Status code:", status_code(page_content), "\n")
    return(invisible(NULL))
  }
}

# process_biocarta() -----------------------------------------------------------
# Downloads the BioCarta GMT from the MSigDB c2.cp.biocarta collection
# (latest release, auto-discovered), and maps symbols to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_biocarta <- function(term_dir, genes, symbol_ensembl_db, save_tsv = FALSE) {
  index_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"

  if (!httr::http_error(index_url)) {
    page         <- rvest::read_html(index_url)
    links        <- html_nodes(page, "a") %>% html_attr("href")
    version_dirs <- grep("^[0-9]{4}\\.[0-9]+\\.Hs/$", links, value = TRUE)

    if (length(version_dirs) == 0) {
      cat("Error: Could not identify MSigDB release version from index.\n")
      return(invisible(NULL))
    }

    latest      <- sort(version_dirs, decreasing = TRUE)[1]
    version_str <- sub("/$", "", latest)
    gmt_url     <- paste0(index_url, latest, "c2.cp.biocarta.v", version_str, ".symbols.gmt")

    temp_gmt <- tempfile(fileext = ".gmt")
    download.file(gmt_url, temp_gmt, mode = "wb")

    if (!file.exists(temp_gmt) || file.size(temp_gmt) == 0) {
      cat("Error: BioCarta GMT download is empty.\n")
      return(invisible(NULL))
    }

    lines <- readLines(temp_gmt)
    unlink(temp_gmt)

    biocarta_df <- do.call(rbind, lapply(lines, function(x) {
      parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
      data.frame(Term = parts[1], Gene = parts[-(1:2)], stringsAsFactors = FALSE)
    }))

    biocarta_df <- biocarta_df %>%
      filter(!is.na(Gene) & Gene != "" & !grepl("^[0-9]", Gene))

    biocarta_df <- update_gene_symbols(biocarta_df, genes)
    biocarta_df <- update_gene_to_ensembl(biocarta_df, symbol_ensembl_db)

    if (save_tsv) {
      write_tsv(biocarta_df, paste0(term_dir, "BioCarta_", Sys.Date(), ".tsv"))
    }

    return(biocarta_df)

  } else {
    cat("Error: MSigDB release index unavailable. Please check the URL and try again.\n")
    return(invisible(NULL))
  }
}

# process_msigdb() -------------------------------------------------------------
# Downloads four MSigDB collections (Hallmark, RegulatoryTargets, Immunologic,
# CellTypeSignatures) from the latest release and maps symbols to Ensembl IDs.
#
# Returns: named list with elements "MSigDB-Hallmark", "MSigDB-RegulatoryTargets",
#          "MSigDB-Immunologic", "MSigDB-CellTypeSignatures",
#          each a data.frame(Term, ensembl_gene_id); or invisible(NULL) on error
process_msigdb <- function(term_dir, genes, symbol_ensembl_db, save_tsv = FALSE) {
  index_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"

  if (!httr::http_error(index_url)) {
    page <- rvest::read_html(index_url)
    links <- html_nodes(page, "a") %>% html_attr("href")
    version_dirs <- grep("^[0-9]{4}\\.[0-9]+\\.Hs/$", links, value = TRUE)

    if (length(version_dirs) == 0) {
      cat("Error: Could not identify MSigDB release version from index.\n")
      return(invisible(NULL))
    }

    latest <- sort(version_dirs, decreasing = TRUE)[1]
    version_str <- sub("/$", "", latest)
    base_url <- paste0(index_url, latest)

    collections <- list(
      "MSigDB-Hallmark"           = "h.all",
      "MSigDB-RegulatoryTargets"  = "c3.all",
      "MSigDB-Immunologic"        = "c7.all",
      "MSigDB-CellTypeSignatures" = "c8.all"
    )

    result <- list()

    for (coll_name in names(collections)) {
      gmt_prefix <- collections[[coll_name]]
      gmt_url <- paste0(base_url, gmt_prefix, ".v", version_str, ".symbols.gmt")

      temp_gmt <- tempfile(fileext = ".gmt")
      download.file(gmt_url, temp_gmt, mode = "wb")

      if (!file.exists(temp_gmt) || file.size(temp_gmt) == 0) {
        cat("Error: MSigDB", coll_name, "download is empty.\n")
        next
      }

      lines <- readLines(temp_gmt)
      unlink(temp_gmt)

      msigdb_df <- do.call(rbind, lapply(lines, function(x) {
        parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
        data.frame(Term = parts[1], Gene = parts[-(1:2)], stringsAsFactors = FALSE)
      }))

      msigdb_df <- msigdb_df %>%
        filter(!is.na(Gene) & Gene != "" & !grepl("^[0-9]", Gene))

      msigdb_df <- update_gene_symbols(msigdb_df, genes)
      msigdb_df <- update_gene_to_ensembl(msigdb_df, symbol_ensembl_db)

      if (save_tsv) {
        write_tsv(msigdb_df, paste0(term_dir, coll_name, "_", Sys.Date(), ".tsv"))
      }

      result[[coll_name]] <- msigdb_df
    }

    if (length(result) == 0) return(invisible(NULL))
    return(result)

  } else {
    cat("Error: MSigDB release index unavailable. Please check the URL and try again.\n")
    return(invisible(NULL))
  }
}

# process_netpath() ------------------------------------------------------------
# Downloads the NetPath GMT from the Bader Lab current release and maps
# symbols to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_netpath <- function(term_dir, genes, symbol_ensembl_db, save_tsv = FALSE) {
  base_url <- "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Pathways/"

  if (!httr::http_error(base_url)) {
    page  <- rvest::read_html(base_url)
    links <- rvest::html_nodes(page, "a") %>% rvest::html_attr("href")
    netpath_file <- grep("^Human_NetPath_.*_symbol\\.gmt$", links, value = TRUE)

    if (length(netpath_file) == 0) {
      cat("Error: Could not find NetPath GMT file in Bader Lab directory.\n")
      return(invisible(NULL))
    }

    gmt_url  <- paste0(base_url, netpath_file[1])
    temp_gmt <- tempfile(fileext = ".gmt")
    download.file(gmt_url, temp_gmt, mode = "wb")

    if (!file.exists(temp_gmt) || file.size(temp_gmt) == 0) {
      cat("Error: NetPath GMT download is empty.\n")
      return(invisible(NULL))
    }

    lines <- readLines(temp_gmt)
    unlink(temp_gmt)

    netpath_df <- do.call(rbind, lapply(lines, function(x) {
      parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
      data.frame(Term = parts[1], Gene = parts[-(1:2)], stringsAsFactors = FALSE)
    }))

    netpath_df <- netpath_df %>%
      filter(!is.na(Gene) & Gene != "" & !grepl("^[0-9]", Gene))

    netpath_df <- update_gene_symbols(netpath_df, genes)
    netpath_df <- update_gene_to_ensembl(netpath_df, symbol_ensembl_db)

    if (save_tsv) {
      write_tsv(netpath_df, paste0(term_dir, "NetPath_", Sys.Date(), ".tsv"))
    }

    return(netpath_df)
  } else {
    cat("Error: Bader Lab GeneSets directory unavailable. Please check the URL and try again.\n")
    return(invisible(NULL))
  }
}

# process_reactome() -----------------------------------------------------------
# Downloads the Reactome GMT zip, parses standard GMT format (col1=Term,
# col2=URL skip, col3+=symbols), and maps to Ensembl IDs.
#
# Returns: data.frame(Term, ensembl_gene_id), or invisible(NULL) on error
process_reactome <- function(url, term_dir, genes, symbol_ensembl_db, save_tsv = FALSE) {
  if (!httr::http_error(url)) {
    temp_zip <- tempfile(fileext = ".zip")
    download.file(url, temp_zip, mode = "wb")

    if (!file.exists(temp_zip) || file.size(temp_zip) == 0) {
      cat("Error: Reactome GMT download is empty.\n")
      return(invisible(NULL))
    }

    gmt_path <- unzip(temp_zip, exdir = tempdir())
    unlink(temp_zip)

    lines <- readLines(gmt_path[1])
    reactome_df <- do.call(rbind, lapply(lines, function(x) {
      parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
      data.frame(Term = parts[1], Gene = parts[-(1:2)], stringsAsFactors = FALSE)
    }))
    unlink(gmt_path)

    reactome_df <- reactome_df %>%
      filter(!is.na(Gene) & Gene != "" & !grepl("^[0-9]", Gene))

    reactome_df <- update_gene_symbols(reactome_df, genes)
    reactome_df <- update_gene_to_ensembl(reactome_df, symbol_ensembl_db)

    if (save_tsv) {
      write_tsv(reactome_df, paste0(term_dir, "Reactome_", Sys.Date(), ".tsv"))
    }

    return(reactome_df)
  } else {
    cat("Error: Invalid Reactome source URL. Please check the URL and try again.\n")
    return(invisible(NULL))
  }
}

# ==============================================================================
# Build geneset database
# ==============================================================================
# Run all process_* functions, collect results into a named list, and save.
# Failed sources (returned NULL) are silently excluded from the compiled RDS.
# ==============================================================================

dir.create("rds/common", recursive = TRUE, showWarnings = FALSE)
if (SAVE_TSVS) dir.create(term_dir, showWarnings = FALSE, recursive = TRUE)

today      <- as.character(Sys.Date())
geneset_db <- list()

# ── PathBank ──────────────────────────────────────────────────────────────────
pathbank_url <- "http://www.pathbank.org/downloads/pathbank_all_proteins.csv.zip"
res <- process_pathbank(pathbank_url, term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("PathBank_", today)]] <- res

# ── Jensen Compartments ───────────────────────────────────────────────────────
jencom_url <- "https://download.jensenlab.org/human_compartment_integrated_full.tsv"
res <- process_jensen_compartment(jencom_url, term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("Jensen-Compartments_", today)]] <- res

# ── Gene Ontology ─────────────────────────────────────────────────────────────
go_url <- "https://current.geneontology.org/annotations/goa_human.gaf.gz"
res <- process_go(go_url, term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) {
  names(res) <- paste0(names(res), "_", today)
  geneset_db <- c(geneset_db, res)
}

# ── KEGG ──────────────────────────────────────────────────────────────────────
res <- process_kegg(term_dir, genes, symbol_entrez_db, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("KEGG_", today)]] <- res

# ── WikiPathways ──────────────────────────────────────────────────────────────
wiki_url <- "https://data.wikipathways.org/current/gmt/"
res <- process_wiki(wiki_url, term_dir, genes, symbol_entrez_db, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("WikiPathways_", today)]] <- res

# ── Panther ───────────────────────────────────────────────────────────────────
panther_url <- "http://data.pantherdb.org/ftp/pathway/current_release/"
res <- process_panther(panther_url, term_dir, hgnc_symbol_db, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("Panther_", today)]] <- res

# ── NCI-PID ───────────────────────────────────────────────────────────────────
res <- process_nci(term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("NCI-PID_", today)]] <- res

# ── ACSN ──────────────────────────────────────────────────────────────────────
acsn_url <- "https://acsn.curie.fr/ACSN2/downloads.html"
res <- process_acsn(acsn_url, term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("ACSN_", today)]] <- res

# ── Reactome ──────────────────────────────────────────────────────────────────
reactome_url <- "https://reactome.org/download/current/ReactomePathways.gmt.zip"
res <- process_reactome(reactome_url, term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("Reactome_", today)]] <- res

# ── BioCarta ──────────────────────────────────────────────────────────────────
res <- process_biocarta(term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("BioCarta_", today)]] <- res

# ── MSigDB (Hallmark, RegulatoryTargets, Immunologic, CellTypeSignatures) ─────
res <- process_msigdb(term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) {
  names(res) <- paste0(names(res), "_", today)
  geneset_db <- c(geneset_db, res)
}

# ── HumanCyc ──────────────────────────────────────────────────────────────────
res <- process_humancyc(term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("HumanCyc_", today)]] <- res

# ── NetPath ───────────────────────────────────────────────────────────────────
res <- process_netpath(term_dir, genes, symbol_ensembl_db, save_tsv = SAVE_TSVS)
if (!is.null(res)) geneset_db[[paste0("NetPath_", today)]] <- res

message("Compiled geneset_db with ", length(geneset_db), " sources: ",
        paste(names(geneset_db), collapse = ", "))

if (SAVE_OUTPUT_RDS) {
  saveRDS(geneset_db, file = output_rds[["geneset_database"]])
  message("Saved geneset_database.rds to ", output_rds[["geneset_database"]])
}