# =============================================================================
# 00_functions.R
# =============================================================================
# Purpose  : Shared utility functions sourced across all chapter scripts.
# Project  : PhD Thesis — Microglial Identity
# Author   : Szymon
# Updated  : 2026
#
# Usage    : source("microglial-identity/common/R/00_functions.R")
#
# Notes    : - No function bodies are modified here; this file is the
#              canonical, annotated version for the reproducibility package.
#            - Functions are grouped by analytical stage (see Table of
#              Contents below).
#            - map_terms_to_genes() is a single consolidated function that
#              handles both ORA (enricher) and fgsea output. Pass
#              test = "fgsea" when the input has a 'pathway' column instead
#              of 'ID'; omit it (or pass NULL) for standard ORA results.
#
# Dependencies (not exhaustive — individual scripts may require additional):
#   edgeR, DESeq2, SummarizedExperiment, org.Hs.eg.db, AnnotationDbi,
#   sva, DaMiR2, WGCNA, ggplot2, ggrepel, patchwork, ggdendro,
#   dplyr, tidyr, purrr, tibble, broom, data.table, Hmisc,
#   Matrix, matrixStats, progress, pheatmap, qvalue, ineq
# =============================================================================


# =============================================================================
# TABLE OF CONTENTS
# =============================================================================
#  1. Quality Control ..................................... line ~60
#     plot_qc_boxplot()
#     qc_bar()
#     plot_qc_density()
#     SexDetect()
#     custom_DaMiR.sampleFilt()
#
#  2. Normalisation ....................................... line ~400
#     DaMiR_normalisation_custom()
#
#  3. Surrogate Variable Analysis ......................... line ~500
#     DaMiR.SV_modified()
#
#  4. Dimensionality Reduction / PCA ..................... line ~610
#     plot_pca()
#     twplot_pca()
#     plot_scree()
#     plot_pca_loadings()
#     compute_correlation()
#     compute_kruskal_wallis()
#
#  5. WGCNA ............................................... line ~800
#     Kmodule_sparse_progress()
#     restructureModules()
#     create_quantile_ranked_heatmap()
#     .wgcna_theme()          [internal style helper]
#     plot_soft_threshold()
#     plot_wgcna_dendrogram()
#
#  6. Enrichment Analysis ................................. line ~1100
#     custom_enricher()
#     add_odds_ratio()
#     map_terms_to_genes()    [ORA / enricher variant]
#     map_terms_to_genes()    [fgsea variant — overwrites above]
#     summarize_term_gene_list()
#     jaccard_index()
#     calculate_jaccard_index_with_expression()
#     run_enrichment_analysis()
#     count_gene_set_terms()
#
#  7. Visualisation Utilities ............................. line ~1450
#     modify_labels()
#     gt_table()
# =============================================================================


# =============================================================================
# 1. QUALITY CONTROL
# =============================================================================

# plot_qc_boxplot() ------------------------------------------------------------
# Draws a per-sample/per-group boxplot with overlaid jittered points for a
# single QC metric. Points below `threshold` are coloured black to flag them
# as potential outliers.
#
# Two colouring modes are supported:
#
#   Mode A — cell-type palette (colours provided):
#     Pass a data frame with columns 'cell_type' and 'color' via `colours`.
#     Above-threshold points inherit the cell-type colour from that lookup;
#     below-threshold points are overridden to black. data_frame must have a
#     'cell_type' column to join on.
#
#   Mode B — simple two-colour (colours = NULL, default):
#     Above-threshold points are coloured "darkgreen"; below-threshold points
#     are coloured "black". No cell-type metadata is required.
#
# An optional dashed horizontal threshold line is added when `threshold` is
# not NULL.
#
# Note on group label formatting:
#   This function does NOT modify group labels internally. Any label
#   formatting (e.g. stripping database accession suffixes) should be applied
#   to the data frame before calling this function. This ensures that all
#   layers — boxplots, jittered points, and any external overlays such as
#   triangle markers — use consistent labels and align correctly on the x-axis.
#
# Args:
#   data_frame   : data frame containing QC metrics; must have a 'cell_type'
#                  column when `colours` is provided
#   metric       : character — name of the numeric QC column to plot
#   group_by     : character — name of the column to use as the x-axis grouping
#   threshold    : numeric — value below which points are flagged (NULL = no line)
#   colours      : optional data frame with columns 'cell_type' and 'color'
#                  (Mode A); NULL activates Mode B (default)
#   sort_by_mean : logical — if TRUE (default), groups are ordered on the
#                  x-axis by descending mean metric value
#
# Returns: a ggplot object
plot_qc_boxplot <- function(data_frame, metric, group_by, threshold = NULL,
                            colours = NULL, sort_by_mean = TRUE) {
  
  if (is.null(data_frame) || is.null(metric) || is.null(group_by)) {
    stop("Please provide a data frame, metric, and group_by arguments")
  }
  if (!is.null(colours) && !("cell_type" %in% colnames(data_frame))) {
    stop("'colours' was provided but data_frame has no 'cell_type' column to join on")
  }
  
  metric_sym   <- rlang::sym(metric)
  group_by_sym <- rlang::sym(group_by)
  
  # Optionally order x-axis groups by descending mean metric
  if (sort_by_mean) {
    mean_metric <- data_frame %>%
      group_by(!!group_by_sym) %>%
      summarise(mean_metric = mean(!!metric_sym, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_metric))
    data_frame[[group_by]] <- factor(data_frame[[group_by]],
                                     levels = mean_metric[[group_by]])
  }
  
  # Assign point colours -------------------------------------------------------
  if (!is.null(colours)) {
    # Mode A: join cell-type palette, then override below-threshold with black
    data_frame <- left_join(data_frame, colours, by = "cell_type")
    if (!is.null(threshold)) {
      data_frame <- mutate(data_frame,
                           color = ifelse(!!metric_sym < threshold, "black", color))
    }
  } else {
    # Mode B: simple black / darkgreen scheme
    if (!is.null(threshold)) {
      data_frame$color <- ifelse(data_frame[[metric]] < threshold, "black", "darkgreen")
    } else {
      data_frame$color <- "black"
    }
  }
  
  p <- ggplot(data_frame, aes_string(x = group_by, y = metric)) +
    geom_boxplot(fill = NA, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.75, size = 1.5, aes(color = color)) +
    scale_color_identity() +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x  = element_line(color = "black", size = 0.2),
          axis.line.y  = element_line(color = "black", size = 0.2))
  
  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept = threshold, linetype = "dashed",
                        color = "black", size = 0.4)
  }
  
  p <- p + theme(axis.text.x  = element_text(angle = 45, hjust = 1),
                 plot.margin  = margin(5.5, 5.5, 5.5, 50, "pt"))
  
  return(p)
}

# qc_bar() ---------------------------------------------------------------------
# Draws a stacked bar chart showing the count of pathological vs.
# non-pathological samples within each group. Bars are ordered by total
# non-pathological count (descending). Intended as a companion plot to
# plot_qc_boxplot() for discrete tissue-type QC flags.
#
# Args:
#   data_frame : data frame with a logical/boolean `metric` column and a
#                grouping column
#   metric     : character — name of the logical column (TRUE = non-pathological)
#   group_by   : character — name of the grouping column
#
# Returns: a ggplot object
qc_bar <- function(data_frame, metric, group_by) {
  if(is.null(data_frame) || is.null(metric) || is.null(group_by)) {
    stop("Please provide a data frame, metric, and group_by arguments")
  }
  
  data_frame[[group_by]] <- sub("_[^_]*$", "", data_frame[[group_by]])
  metric_sym <- rlang::sym(metric)
  group_by_sym <- rlang::sym(group_by)
  
  total_controls <- data_frame %>%
    group_by(!!group_by_sym) %>%
    summarise(total_controls = sum(!!metric_sym, na.rm = TRUE)) %>%
    arrange(desc(total_controls))
  
  data_frame[[group_by]] <- factor(data_frame[[group_by]], levels = total_controls[[group_by]])
  
  p <- ggplot(data_frame, aes_string(x = group_by, fill = metric_sym)) +
    geom_bar(position = "stack", alpha = 0.5) +
    scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "black"),
                      labels = c("Pathological", "Non-pathological"),
                      name = "Tissue Type") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2))
  
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 plot.margin = margin(5.5, 5.5, 5.5, 50, "pt"))
  
  return(p)
}


# plot_qc_density() ------------------------------------------------------------
# Draws a density curve for a single QC metric across all samples, with an
# optional vertical dashed threshold line. Annotates the number of samples
# falling below (Fail) and at/above (Pass) the threshold.
#
# Args:
#   data_frame : data frame with at least one numeric column matching `metric`
#   metric     : character — name of the numeric column to plot
#   threshold  : numeric — threshold for pass/fail annotation (NULL = no line)
#
# Returns: a ggplot object
plot_qc_density <- function(data_frame, metric, threshold) {
  if(is.null(data_frame) || is.null(metric)) {
    stop("Please provide a data frame and metric arguments")
  }
  
  metric_sym <- rlang::sym(metric)
  
  p <- ggplot(data_frame, aes(x = !!metric_sym)) +
    geom_density(color = "black") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2))
  
  density_data <- density(data_frame[[metric_sym]], na.rm = TRUE)
  max_density_value <- max(density_data$y, na.rm = T)
  
  if (!is.null(threshold)) {
    p <- p + geom_vline(xintercept = threshold, linetype = "dashed", color = "black", size = 0.4)
    count_below <- sum(data_frame[[metric]] < threshold, na.rm = T)
    count_above <- sum(data_frame[[metric]] >= threshold, na.rm = T)
    p <- p + annotate("text", x = (threshold + min(data_frame[[metric]], na.rm = T)) / 2, y = max_density_value * 0.9, label = paste("Fail:", count_below))
    p <- p + annotate("text", x = (max(data_frame[[metric]], na.rm = T) + threshold) / 2, y = max_density_value * 0.9, label = paste("Pass:", count_above))
  }
  
  p <- p + ylab("Density")
  return(p)
}


# SexDetect() ------------------------------------------------------------------
# Imputes biological sex for each sample using CPM expression of Y-chromosome
# genes. PCA is performed on log2-transformed Y-gene CPM values; samples are
# classified as Male or Female based on the sign of PC1. The imputed sex is
# appended to the metadata table as a new column 'sex_imputed' (inserted after
# the existing 'sex' column).
#
# Args:
#   data     : a count matrix (rows = genes with Ensembl IDs, cols = samples)
#   metadata : data frame with a 'sex' column and rownames matching colnames(data)
#
# Returns: metadata data frame with an additional 'sex_imputed' column
#
# Note: Requires org.Hs.eg.db for chromosome mapping and calls plot_pca()
#       (defined below), which assigns `pca` to the global environment.
SexDetect <- function(data, metadata) {
  cpm_counts <- cpm(data)
  
  rows_chrom <- mapIds(org.Hs.eg.db, column = c('CHR'), keytype = 'ENSEMBL', keys = gsub('\\..*','',rownames(cpm_counts)), unique = T)
  rows_y <- which(rows_chrom == "Y")
  y_cpm_counts <- cpm_counts[rows_y,]
  
  meow <- plot_pca(log2(y_cpm_counts + 0.5), metadata, intgroup = "cell_type", ntop = nrow(y_cpm_counts), scale = FALSE)
  
  sex_data <- as.data.frame(pca$x)
  sex_imputed <- ifelse(sex_data$PC1 > 0, "Male", "Female")
  metadata <- add_column(metadata, sex_imputed = sex_imputed, .after = "sex")
  
  return(metadata)
}


# custom_DaMiR.sampleFilt() ----------------------------------------------------
# Multi-study sample QC filter extending DaMiR2::DaMiR.sampleFilt().
# Applies four QC criteria evaluated SIMULTANEOUSLY on the full sample set.
# All metrics are computed before any sample is removed, so each sample is
# assessed against the same reference distribution. Samples failing any
# criterion are removed in a single pass at the end.
#
#   1. Study-size-weighted inter-study connectivity Z-score
#      Each study contributes one equal vote regardless of sample count,
#      preventing large studies from dominating the connectivity estimate.
#      Computed on the full correlation matrix.
#
#   2. Within-study Z-score
#      Outlier detection is relative to each study's own distribution.
#      Studies with <= small_study_threshold samples are exempt (too few
#      samples to compute a stable Z-score). Computed on the full
#      correlation matrix.
#
#   3. Near-duplicate detection
#      Sample pairs with pairwise correlation >= th.duplicate_cor are
#      flagged. Within each duplicate pair, the sample with the lower mean
#      connectivity to ALL other samples is removed. Evaluated on the full
#      correlation matrix.
#
#   4. WGCNA standardised connectivity (Z.k)
#      Global connectivity QC computed on the full sample set (not on
#      survivors of the preceding filters, as in the original sequential
#      implementation). Catches samples with globally low connectivity
#      that may not be flagged by study-level criteria.
#
# Note on simultaneous vs sequential evaluation:
#   In the original sequential implementation, metrics in steps 2–4 were
#   recomputed on survivors of the preceding step. This meant that samples
#   removed early never received later metrics (leaving NAs in the QC table)
#   and that the reference distribution shifted between steps. The simultaneous
#   approach uses a single correlation matrix computed once on all samples,
#   ensuring every sample has a complete set of QC metrics and is evaluated
#   against the same reference. Z.k values will therefore differ slightly from
#   the sequential implementation because they are no longer computed on a
#   post-filter subset.
#
# Produces five diagnostic ggplots (printed automatically) and returns a
# named list with the filtered SummarizedExperiment, per-sample QC table,
# WGCNA connectivity values, and duplicate-pair metadata.
#
# Args:
#   data                  : SummarizedExperiment (normalised counts)
#   study_df              : data frame with a 'study' column; rownames = sample IDs
#   th.inter_study_z      : Z-score cutoff for inter-study connectivity (default -2)
#   th.within_study_z     : Z-score cutoff for within-study connectivity (default -2)
#   th.wgcna_z            : Z.k cutoff for WGCNA connectivity (default -2.5)
#   th.duplicate_cor      : pairwise correlation threshold for duplicate detection (default 0.99)
#   small_study_threshold : studies with <= this many samples skip within-study Z (default 3)
#   type                  : correlation method — "spearman" (default) or "pearson"
#
# Returns: list with elements:
#   $filtered_data     — SummarizedExperiment, samples passing all filters
#   $sample_qc_summary — tibble, full per-sample QC metrics and filter outcomes
#                        (all four metrics populated for every sample)
#   $wgcna_connectivity— tibble, mean IAC and Z.k for all samples
#   $keep_samples      — character vector of retained sample IDs
#   $duplicate_pairs   — data frame of all flagged near-duplicate pairs
#   $removed_duplicates— character vector of samples removed as duplicates
custom_DaMiR.sampleFilt <- function(data,
                                    study_df,
                                    th.inter_study_z      = -2,
                                    th.within_study_z     = -2,
                                    th.wgcna_z            = -2.5,
                                    th.duplicate_cor      = 0.99,
                                    small_study_threshold = 3,
                                    type                  = c("spearman", "pearson")) {
  
  # ── Input checks ────────────────────────────────────────────────────────────
  if (missing(data))     stop("'data' argument must be provided")
  if (missing(study_df)) stop("'study_df' argument must be provided")
  if (!(is(data, "SummarizedExperiment")))
    stop("'data' must be a 'SummarizedExperiment' object")
  if (!("class" %in% colnames(colData(data))))
    stop("'class' info is lacking! Include the variable 'class' in colData(data).")
  if (!("study" %in% colnames(study_df)))
    stop("'study_df' must contain a 'study' column.")
  if (any(is.na(assay(data))))       stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(assay(data)))) stop("Inf values are not allowed in the 'data' matrix")
  if (all(assay(data) == 0))         stop("All genes have 0 values")
  if (all((assay(data) %% 1) == 0))
    warning("It seems you are using raw counts! This function works with normalised data.")
  if (missing(type)) type <- type[1]
  if (length(type) > 1) stop("length(type) must be equal to 1")
  if (!(type %in% c("pearson", "spearman")))
    stop("'type' must be 'pearson' or 'spearman'")
  
  count_data     <- assay(data)
  n_samples_init <- ncol(count_data)
  all_samples    <- colnames(count_data)
  
  cat("=== custom_DaMiR.sampleFilt: multi-study QC (simultaneous mode) ===\n")
  cat("Starting with", n_samples_init, "samples.\n\n")
  
  # ── Step 1: Full correlation matrix on ALL samples ───────────────────────────
  # Computed once on the complete sample set. All downstream metrics derive
  # from this single matrix, ensuring a consistent reference for every sample.
  cat("Step 1: Computing inter-sample correlation matrix (", type, ") on full sample set...\n")
  if (type == "spearman") {
    cormatrix <- abs(Hmisc::rcorr(as.matrix(count_data), type = "spearman")$r)
  } else {
    cormatrix <- abs(Hmisc::rcorr(as.matrix(count_data), type = "pearson")$r)
  }
  diag(cormatrix) <- NA   # exclude self-correlation from all means
  
  # ── Step 2: Study membership lookup ─────────────────────────────────────────
  study_vec   <- setNames(study_df[all_samples, "study"], all_samples)
  study_sizes <- table(study_vec)
  cat("Study sizes:\n"); print(study_sizes); cat("\n")
  
  # ── Step 3: Study-size-weighted inter-study connectivity ─────────────────────
  # For each sample, compute mean correlation to each OTHER study separately,
  # then average those per-study means — giving every study an equal vote
  # regardless of sample count. Result is Z-scored globally across all samples.
  cat("Step 2: Computing study-size-weighted inter-study connectivity on full sample set...\n")
  unique_studies <- unique(study_vec)
  
  inter_study_conn <- sapply(all_samples, function(focal) {
    focal_study   <- study_vec[focal]
    other_studies <- unique_studies[unique_studies != focal_study]
    if (length(other_studies) == 0) return(NA)
    per_study_means <- sapply(other_studies, function(s) {
      s_samples <- names(study_vec[study_vec == s])
      mean(cormatrix[focal, s_samples], na.rm = TRUE)
    })
    mean(per_study_means, na.rm = TRUE)
  })
  inter_study_z        <- as.numeric(scale(inter_study_conn))
  names(inter_study_z) <- all_samples
  
  # ── Step 4: Within-study Z-score on ALL samples ──────────────────────────────
  # Mean correlation to own study-mates, Z-scored within each study.
  # Studies with <= small_study_threshold samples are exempt: too few
  # observations for a stable Z-score, so within_study_z is set to NA
  # and these samples automatically pass the within-study filter.
  cat("Step 3: Computing within-study Z-scores on full sample set...\n")
  within_study_mean_corr <- sapply(all_samples, function(focal) {
    focal_study <- study_vec[focal]
    study_mates <- names(study_vec[study_vec == focal_study & names(study_vec) != focal])
    if (length(study_mates) == 0) return(NA)
    mean(cormatrix[focal, study_mates], na.rm = TRUE)
  })
  within_study_z <- ave(
    within_study_mean_corr,
    study_vec,
    FUN = function(x) {
      if (length(x) <= 2 || all(is.na(x))) return(rep(NA, length(x)))
      as.numeric(scale(x))
    }
  )
  names(within_study_z) <- all_samples
  
  # ── Step 5: WGCNA standardised connectivity (Z.k) on ALL samples ────────────
  # Mean inter-array correlation (IAC) computed from the full correlation matrix,
  # then Z-scored globally. Because this is computed before any removal, Z.k
  # reflects each sample's connectivity to the complete dataset rather than
  # a post-filter subset. Values will therefore differ from a sequential
  # implementation where Z.k is recomputed on survivors only.
  cat("Step 4: Computing WGCNA standardised connectivity (Z.k) on full sample set...\n")
  mean_IAC_full        <- rowMeans(cormatrix, na.rm = TRUE)
  Z.k_full             <- as.numeric(scale(mean_IAC_full))
  names(Z.k_full)      <- all_samples
  
  wgcna_df <- tibble(
    sample_name = all_samples,
    mean_IAC    = mean_IAC_full,
    Z.k         = Z.k_full
  )
  
  # ── Step 6: Near-duplicate detection on ALL samples ──────────────────────────
  # Sample pairs with pairwise correlation >= th.duplicate_cor are flagged.
  # Within each pair, the sample with the lower mean connectivity to all other
  # samples is removed. Pairs are processed in descending correlation order;
  # once a sample is removed it cannot be involved in further comparisons.
  cat("Step 5: Detecting near-duplicate samples (pairwise correlation >=",
      th.duplicate_cor, ") on full sample set...\n")
  
  dup_pairs <- which(cormatrix >= th.duplicate_cor, arr.ind = TRUE)
  dup_pairs <- dup_pairs[dup_pairs[, 1] < dup_pairs[, 2], , drop = FALSE]
  
  dup_pair_df <- if (nrow(dup_pairs) > 0) {
    data.frame(
      sample_a     = rownames(cormatrix)[dup_pairs[, 1]],
      sample_b     = colnames(cormatrix)[dup_pairs[, 2]],
      pairwise_cor = cormatrix[dup_pairs],
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(sample_a = character(), sample_b = character(),
               pairwise_cor = numeric(), stringsAsFactors = FALSE)
  }
  
  mean_conn_full     <- rowMeans(cormatrix, na.rm = TRUE)
  removed_duplicates <- character(0)
  
  if (nrow(dup_pair_df) > 0) {
    cat("  Found", nrow(dup_pair_df), "near-duplicate pair(s):\n")
    dup_pair_df_sorted <- dup_pair_df[order(-dup_pair_df$pairwise_cor), ]
    active_samples     <- all_samples
    
    for (i in seq_len(nrow(dup_pair_df_sorted))) {
      sa <- dup_pair_df_sorted$sample_a[i]
      sb <- dup_pair_df_sorted$sample_b[i]
      if (!(sa %in% active_samples) || !(sb %in% active_samples)) next
      
      if (mean_conn_full[sa] >= mean_conn_full[sb]) {
        drop_sample <- sb; keep_sample_dup <- sa
      } else {
        drop_sample <- sa; keep_sample_dup <- sb
      }
      cat(sprintf("    Pair (r = %.4f): keeping %s, removing %s\n",
                  dup_pair_df_sorted$pairwise_cor[i], keep_sample_dup, drop_sample))
      removed_duplicates <- c(removed_duplicates, drop_sample)
      active_samples     <- active_samples[active_samples != drop_sample]
    }
    cat("  Total removed as near-duplicates:", length(removed_duplicates), "\n")
  } else {
    cat("  No near-duplicate pairs found above threshold", th.duplicate_cor, "\n")
  }
  
  # ── Step 7: Build per-sample QC summary — all metrics populated ──────────────
  # Every sample receives values for all four QC metrics. Samples exempt from
  # the within-study filter (study_n <= small_study_threshold) have NA for
  # within_study_z and automatically pass that criterion.
  qc_df <- tibble(
    sample_name          = all_samples,
    study                = study_vec[all_samples],
    study_n              = as.integer(study_sizes[study_vec[all_samples]]),
    inter_study_conn     = inter_study_conn[all_samples],
    inter_study_z        = inter_study_z[all_samples],
    within_study_corr    = within_study_mean_corr[all_samples],
    within_study_z       = within_study_z[all_samples],
    mean_IAC             = mean_IAC_full[all_samples],
    Z.k                  = Z.k_full[all_samples],
    removed_as_duplicate = all_samples %in% removed_duplicates
  )
  
  # ── Step 8: Apply all filters simultaneously ─────────────────────────────────
  # All four pass/fail decisions are made in a single step using the metrics
  # computed above. No metric is recomputed between filter passes.
  # Within-study exempt samples (NA within_study_z) pass the within-study
  # criterion unconditionally.
  cat("\nStep 6: Applying all filters simultaneously...\n")
  
  pass_inter  <- is.na(qc_df$inter_study_z) | (qc_df$inter_study_z  >= th.inter_study_z)
  pass_within <- qc_df$study_n <= small_study_threshold |
    is.na(qc_df$within_study_z) |
    (qc_df$within_study_z >= th.within_study_z)
  pass_wgcna  <- qc_df$Z.k >= th.wgcna_z
  pass_dedup  <- !qc_df$removed_as_duplicate
  
  qc_df <- qc_df %>%
    mutate(
      pass_correlation_filters = pass_inter & pass_within,
      pass_dedup_filter        = pass_dedup,
      pass_wgcna_filter        = pass_wgcna,
      pass_all_filters         = pass_inter & pass_within & pass_dedup & pass_wgcna
    )
  
  keep_samples_final <- qc_df$sample_name[qc_df$pass_all_filters]
  
  cat("  Removed by inter-study Z filter   :", sum(!pass_inter), "\n")
  cat("  Removed by within-study Z filter  :", sum(!pass_within & pass_inter), "\n")
  cat("  Removed as near-duplicates        :", length(removed_duplicates), "\n")
  cat("  Removed by WGCNA Z.k filter       :", sum(!pass_wgcna & pass_inter & pass_within & pass_dedup), "\n")
  cat("  Samples retained                  :", length(keep_samples_final), "\n\n")
  
  # ── Step 9: Filter the SummarizedExperiment ─────────────────────────────────
  data_filt <- data[, keep_samples_final]
  
  cat("=== Summary ===\n")
  cat("  Initial samples :", n_samples_init, "\n")
  cat("  Final samples   :", length(keep_samples_final), "\n")
  cat("  Total removed   :", n_samples_init - length(keep_samples_final), "\n\n")
  
  # ── Diagnostic plots ────────────────────────────────────────────────────────
  p1 <- ggplot(qc_df, aes(x = inter_study_z, fill = pass_correlation_filters)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "white") +
    geom_vline(xintercept = th.inter_study_z, linetype = "dashed",
               color = "#AA4643", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#44BB99", "FALSE" = "#AA4643"),
                      labels = c("TRUE" = "Pass", "FALSE" = "Fail"), name = "Filter") +
    labs(title = "Inter-study connectivity (study-weighted Z-score)",
         x = "Z-score", y = "Count") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2))
  
  p2 <- ggplot(qc_df, aes(x = within_study_z, fill = pass_correlation_filters)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "white") +
    geom_vline(xintercept = th.within_study_z, linetype = "dashed",
               color = "#AA4643", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#44BB99", "FALSE" = "#AA4643"),
                      labels = c("TRUE" = "Pass", "FALSE" = "Fail"), name = "Filter") +
    labs(title = paste0("Within-study Z-score (studies with n \u2264 ",
                        small_study_threshold, " are exempt)"),
         x = "Within-study Z-score", y = "Count") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2))
  
  p3 <- ggplot(qc_df, aes(x = Z.k, fill = pass_wgcna_filter)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "white") +
    geom_vline(xintercept = th.wgcna_z, linetype = "dashed",
               color = "#AA4643", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#44BB99", "FALSE" = "#AA4643"),
                      labels = c("TRUE" = "Pass", "FALSE" = "Fail"), name = "Filter") +
    labs(title = "WGCNA standardised connectivity (Z.k) — computed on full sample set",
         x = "Z.k", y = "Count") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2))
  
  p4 <- ggplot(qc_df,
               aes(x = reorder(study, within_study_corr, FUN = median, na.rm = TRUE),
                   y = within_study_corr,
                   fill = factor(study_n <= small_study_threshold,
                                 levels = c(FALSE, TRUE),
                                 labels = c("Normal", "Small (exempt)")))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.8,
                aes(color = pass_correlation_filters)) +
    scale_color_manual(values = c("TRUE" = "#44BB99", "FALSE" = "#AA4643"),
                       labels = c("TRUE" = "Pass", "FALSE" = "Fail"),
                       name = "Corr. filter") +
    scale_fill_manual(values = c("Normal" = "#77AADD", "Small (exempt)" = "#EEDD88"),
                      name = "Study size") +
    labs(title = "Within-study correlation by study",
         x = "Study", y = "Mean within-study correlation") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  max_pairwise_cor <- apply(cormatrix, 1, max, na.rm = TRUE)
  dup_plot_df <- data.frame(
    sample_name      = names(max_pairwise_cor),
    max_pairwise_cor = max_pairwise_cor,
    is_duplicate     = names(max_pairwise_cor) %in% removed_duplicates,
    stringsAsFactors = FALSE
  )
  p_dup <- ggplot(dup_plot_df, aes(x = max_pairwise_cor, fill = is_duplicate)) +
    geom_histogram(bins = 40, alpha = 0.7, color = "white") +
    geom_vline(xintercept = th.duplicate_cor, linetype = "dashed",
               color = "#AA4643", linewidth = 0.6) +
    scale_fill_manual(values = c("TRUE" = "#AA4643", "FALSE" = "#44BB99"),
                      labels = c("TRUE" = "Removed (duplicate)", "FALSE" = "Retained"),
                      name = "Status") +
    labs(title = paste0("Near-duplicate detection (threshold: r \u2265 ", th.duplicate_cor, ")"),
         x = "Max pairwise correlation to any other sample", y = "Count") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2))
  
  print(p1); print(p2); print(p3); print(p4); print(p_dup)
  
  # ── Return ───────────────────────────────────────────────────────────────────
  return(list(
    filtered_data      = data_filt,
    sample_qc_summary  = qc_df,
    wgcna_connectivity = wgcna_df,
    keep_samples       = keep_samples_final,
    duplicate_pairs    = dup_pair_df,
    removed_duplicates = removed_duplicates
  ))
}

# =============================================================================
# 2. NORMALISATION
# =============================================================================

# DaMiR_normalisation_custom() -------------------------------------------------
# Extended wrapper around DaMiR2 normalisation routines. Performs low-count
# gene filtering (rowSums), optional hypervariant gene removal (CV-based), and
# normalisation by one of three methods:
#   - "vst"    : DESeq2::varianceStabilizingTransformation
#   - "rlog"   : DESeq2::rlog (slow for large datasets)
#   - "logcpm" : edgeR::cpm(log = TRUE, prior.count = 1)
#
# Args:
#   data       : SummarizedExperiment with raw integer counts in assay slot
#   minCounts  : minimum count threshold for filtering (default 10)
#   fSample    : minimum fraction of samples that must meet minCounts (default 0.5)
#   hyper      : "yes" or "no" — whether to remove hypervariant genes (default "yes")
#   th.cv      : CV threshold for hypervariant gene removal (default 3)
#   type       : normalisation method — "vst", "rlog", or "logcpm" (default "vst")
#   nFitType   : DESeq2 dispersion fitting — "parametric", "local", or "mean" (default "parametric")
#
# Returns: SummarizedExperiment with normalised values in the assay slot
DaMiR_normalisation_custom <- function(data, minCounts = 10, fSample = 0.5, hyper = c("yes", 
                                                                                      "no"), th.cv = 3, type = c("vst", "rlog", "logcpm"), nFitType = c("parametric", 
                                                                                                                                                        "local", "mean")) {
  
  if (missing(data)) 
    stop("'data' argument must be provided")
  if (missing(hyper)) {
    hyper <- hyper[1]
  }
  if (missing(type)) {
    type <- type[1]
  }
  if (missing(nFitType)) {
    nFitType <- nFitType[1]
  }
  if (!(is(data, "SummarizedExperiment"))) 
    stop("'data' must be a 'SummarizedExperiment' object")
  if (!(is.numeric(minCounts))) 
    stop("'minCounts' must be numeric")
  if (!(is.numeric(fSample))) 
    stop("'fSample' must be numeric")
  if (!(is.numeric(th.cv))) 
    stop("'th.cv' must be numeric")
  if (any(is.na(assay(data)))) 
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(assay(data)))) 
    stop("Inf values are not allowed in the 'data' matrix")
  if (any(assay(data) < 0)) 
    stop("'data' contains negative values")
  if (all(assay(data) == 0)) 
    stop("All genes have 0 counts.")
  if (any((assay(data)%%1) != 0)) 
    stop("Some values are not raw counts. Check the dataset")
  if (all(!(colnames(assay(data)) %in% rownames(colData(data))))) 
    stop("colnames(assay(data)) must be equal to rownames(colData(data))")
  if (!("class" %in% colnames(colData(data)))) 
    stop("'class' info is lacking!\n         Include the variable 'class' in colData(data) and label it 'class'!")
  if (!(nFitType %in% c("parametric", "local", "mean"))) 
    stop("'nFitType' must be 'parametric', 'local' or 'mean'")
  if (length(type) > 1) 
    stop("length(type) must be equal to 1")
  if (!(all(type %in% c("vst", "rlog", "logcpm")))) 
    stop("'type' must be 'vst', 'rlog' or 'logcpm' ")
  init_lenght <- dim(data)[1]
  minSamplesExpr <- round((dim(data)[2]) * fSample)
  data <- data[rowSums(assay(data) >= minCounts) > minSamplesExpr, ]
  exprs_counts <- assay(data)
  cat(init_lenght - dim(data)[1], "Features have been filtered out by espression.", 
      dim(data)[1], "Features remained.", "\n")
  if (hyper == "yes") {
    init_lenght_cv <- dim(data)[1]
    classes <- levels(data@colData$class)
    cv_value <- matrix(nrow = dim(data)[1], ncol = length(classes))
    for (i in seq_len(length(classes))) {
      matr2cv <- exprs_counts[, which(levels(data@colData$class)[i] == 
                                        data@colData$class)]
      cv_value[, i] <- as.matrix(apply(matr2cv, 1, ineq, type = "var"))
    }
    index2keep <- 0
    for (j in seq_len(dim(cv_value)[1])) {
      index2keep[j] <- all(cv_value[j, ] < th.cv)
    }
    exprs_counts <- exprs_counts[which(index2keep == 1), ]
    cat(init_lenght_cv - dim(exprs_counts)[1], "'Hypervariant' Features have been filtered out.", 
        dim(exprs_counts)[1], "Features remained.", "\n")
  } else if (hyper == "no") {
  } else {
    stop("'hyper' must be set with 'yes' or 'no' ")
  }
  if (type == "vst") {
    cat("Performing Normalization by 'vst'", "with dispersion parameter: ", nFitType, "\n")
    data_norm <- varianceStabilizingTransformation(exprs_counts, fitType = nFitType)
  } else if (type == "rlog") {
    cat("Performing Normalization by 'rlog'", "with dispersion parameter: ", nFitType, "\n",
        "For large dataset it could be very time-consuming.", "\n")
    data_norm <- rlog(exprs_counts, fitType = nFitType)
    rownames(data_norm) <- rownames(exprs_counts)
  } else if (type == "logcpm") {
    cat("Performing Normalization by 'log2cpm'", "\n")
    data_norm <- cpm(exprs_counts, log = TRUE, prior.count = 1)
    rownames(data_norm) <- rownames(exprs_counts)
  } else {
    stop("Please set 'vst or 'rlog' as normalization type.")
  }
  data_norm <- SummarizedExperiment(assays = as.matrix(data_norm), 
                                    colData = as.data.frame(colData(data)))
  return(data_norm)
}


# =============================================================================
# 3. SURROGATE VARIABLE ANALYSIS
# =============================================================================

# DaMiR.SV_modified() ----------------------------------------------------------
# Modified version of DaMiR2::DaMiR.SV(). Identifies surrogate variables (SVs)
# using the sva package to capture unmeasured sources of variation (e.g. batch,
# technical noise). Three methods for determining the number of SVs are
# supported: "fve" (fraction of variance explained, default), "leek", and "be".
#
# Improvements over the original:
#   - 'second.var' can be excluded from the model without erroring
#   - Both 'class' and 'second.var' are validated for sufficient levels before
#     being added to the model
#   - A ggplot scree plot is returned alongside the SV matrix (for "fve" only)
#
# Args:
#   data        : SummarizedExperiment (normalised counts)
#   method      : SV detection method — "fve" (default), "leek", or "be"
#   th.fve      : cumulative variance threshold for "fve" method (default 0.95)
#   second.var  : optional additional covariate to protect in the full model
#   include.class : logical — whether to include 'class' in the model (default TRUE)
#
# Returns: list with elements:
#   $sv_matrix — matrix of SVs (rows = samples, cols = SVs), or NULL if none found
#   $plot      — ggplot scree plot (fve method only), or NULL
DaMiR.SV_modified <- function (data, method = c("fve", "leek", "be"), th.fve = 0.95, 
                               second.var = NULL, include.class = TRUE) 
{
  if (missing(data)) 
    stop("'data' argument must be provided")
  if (missing(method)) {
    method <- method[1]
  }
  if (!(is(data, "SummarizedExperiment"))) 
    stop("'data' must be a 'SummarizedExperiment' object")
  if (!(is.numeric(th.fve))) 
    stop("'th.fve' must be numeric")
  if (any(is.na(assay(data)))) 
    stop("NA values are not allowed in the 'data' matrix")
  if (any(is.infinite(assay(data)))) 
    stop("Inf values are not allowed in the 'data' matrix")
  if (all(assay(data) == 0)) 
    stop("All genes have 0 counts.")
  if (all((assay(data)%%1) == 0)) 
    stop("It seems that you are using raw counts!\nThis function works with normalized data")
  if (!("class" %in% colnames(colData(data)))) 
    stop("'class' info is lacking!\nInclude the variable 'class' in colData(data) and label it 'class'!")
  if (th.fve > 1 | th.fve < 0) 
    stop("'th.fve' must be between 0 and 1")
  data_filt <- assay(data)
  
  variables <- list()
  if (include.class) {
    class_var <- colData(data)$class
    if (is.factor(class_var) && nlevels(class_var) < 2) {
      warning("'class' variable has less than 2 levels and will be excluded from the model.")
    } else if (is.character(class_var) && length(unique(class_var)) < 2) {
      warning("'class' variable has less than 2 unique values and will be excluded from the model.")
    } else {
      variables$class <- class_var
    }
  }
  if (!is.null(second.var)) {
    if (is.factor(second.var) && nlevels(second.var) < 2) {
      warning("'second.var' has less than 2 levels and will be excluded from the model.")
    } else if (is.numeric(second.var) && length(unique(second.var)) < 2) {
      warning("'second.var' is constant and will be excluded from the model.")
    } else {
      variables$second.var <- second.var
    }
  }
  
  if (length(variables) == 0) {
    mod <- model.matrix(~1, data = colData(data))
  } else {
    mod_data <- as.data.frame(variables)
    mod <- model.matrix(~., data = mod_data)
  }
  
  mod0 <- cbind(mod[, 1])
  
  if (method == "fve") {
    if (round(ncol(data)) > 100) {
      n.sv.max <- 50
    } else {
      n.sv.max <- round(ncol(data)/2)
    }
    invisible(capture.output(svaobj <- sva(data_filt, mod, mod0, n.sv.max)))
    pprob <- svaobj$pprob.gam * (1 - svaobj$pprob.b)
    dats <- data_filt * pprob
    dats <- dats - rowMeans(dats)
    uu <- eigen(t(dats) %*% dats)
    uu_val2 <- uu$values^2/sum(uu$values^2)
    x_val <- seq_len(ncol(dats))
    uu_val2_dec <- round(uu_val2, 3)
    pve_cumsum <- cumsum(uu_val2_dec)
    expl_var_plot <- as.data.frame(cbind(x_val, uu_val2, pve_cumsum))
    n.sv <- order(which(pve_cumsum <= th.fve), decreasing = TRUE)[1]
    if (is.na(n.sv) | (n.sv == 0)) {
      sv_matrix <- NULL
      p <- NULL
      cat("No SV identified.\n")
    } else {
      if (n.sv > n.sv.max) n.sv <- n.sv.max
      expl_var_plot <- expl_var_plot[seq_len(n.sv + 5),]
      p <- ggplot(expl_var_plot, aes(x_val, pve_cumsum)) + 
        geom_point(color = "black") + 
        geom_text(aes(label = x_val), hjust = 0.5, vjust = 2, size = 3) + 
        geom_point(data = expl_var_plot[n.sv, ], aes(x_val, pve_cumsum), color = "#cf3434", size = 3) + 
        xlab("SV") + ylab("Fraction of Variance Explained") + 
        ggtitle("Fraction of Variance Explained") + 
        theme_bw() + 
        theme(panel.border = element_blank(),
              axis.line.x = element_line(color = "black", size = 0.2),
              axis.line.y = element_line(color = "black", size = 0.2))
      sv_matrix <- as.matrix(svaobj$sv[, seq_len(n.sv)])
      cat("The number of SVs identified, which explain", th.fve * 100, "% of Variance, is:", n.sv, "\n")
    }
  } else if (method == "leek") {
    n.sv <- num.sv(data_filt, mod = mod, method = "leek")
    if (is.na(n.sv) | (n.sv == 0)) {
      sv_matrix <- NULL
      p <- NULL
      cat("No SV identified.\n")
    } else {
      invisible(capture.output(svaobj <- sva(data_filt, mod, mod0, n.sv)))
      sv_matrix <- as.matrix(svaobj$sv)
      p <- NULL
      cat("The number of SVs identified by SVA (with method = 'leek') is:", n.sv, "\n")
    }
  } else if (method == "be") {
    n.sv <- num.sv(data_filt, mod = mod, method = "be")
    if (is.na(n.sv) | (n.sv == 0)) {
      sv_matrix <- NULL
      p <- NULL
      cat("No SV identified.\n")
    } else {
      invisible(capture.output(svaobj <- sva(data_filt, mod, mod0, n.sv)))
      sv_matrix <- as.matrix(svaobj$sv)
      p <- NULL
      cat("The number of SVs identified by SVA (with method = 'be') is:", n.sv, "\n")
    }
  } else {
    stop("Please set 'fve', 'leek', or 'be' as SV identification method.")
  }
  return(list(sv_matrix = sv_matrix, plot = p))
}


# =============================================================================
# 4. DIMENSIONALITY REDUCTION / PCA
# =============================================================================

# plot_pca() -------------------------------------------------------------------
# General-purpose PCA plot. Selects the top `ntop` most variable genes
# (by row variance), runs prcomp(), and plots the chosen PC pair coloured by
# a metadata variable. Categorical variables use discrete colours; continuous
# variables use a red-white-blue gradient.
#
# Side effects: assigns `pca` and `percentVar` to the global environment so
#               downstream functions (e.g. SexDetect, plot_scree) can use them.
#
# Args:
#   object           : numeric matrix (rows = genes, cols = samples)
#   metadata         : data frame with sample-level metadata; rownames = colnames(object)
#   intgroup         : character — metadata column to use for colouring (default "condition")
#   ntop             : number of most-variable genes to use (default 500)
#   pcs              : integer vector of length 2 — which PCs to plot (default c(1,2))
#   returnData       : if TRUE, return the PC data frame instead of a plot
#   annotate_points  : character vector of sample names to label with ggrepel
#   scale            : logical — whether to scale variables before PCA (default TRUE)
#
# Returns: ggplot object, or a data frame (if returnData = TRUE)
plot_pca <- function (object, metadata, intgroup = "condition", ntop = 500, pcs = c(1, 2), returnData = FALSE, annotate_points = NULL, scale = TRUE) {
  rv <- matrixStats::rowVars(object, useNames = TRUE)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(object[select, ]), scale. = scale)
  assign("pca", pca, envir = globalenv())
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  assign("percentVar", percentVar, envir = globalenv())
  
  if (!all(intgroup %in% colnames(metadata))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- metadata[, intgroup, drop = FALSE]
  is_categorical <- is.factor(metadata[[intgroup]]) || is.character(metadata[[intgroup]])
  
  if (is_categorical) {
    group <- if (ncol(intgroup.df) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
      metadata[[intgroup]]
    }
    groupColor <- "group"
  } else {
    group <- metadata[[intgroup]]
    groupColor <- intgroup
  }
  
  d <- data.frame(pc_1 = pca$x[, pcs[1]], pc_2 = pca$x[, pcs[2]], group = group, name = colnames(object)) %>% cbind(metadata)
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcs[1]:pcs[2]]
    return(d)
  }
  
  colours <- c("#77AADD", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00", 
               "#EEDD88", "#EE8866", "#FFAABB", "#DDDDDD", "#CC99CC")  
  vivid_colours <- c("#AA4643", "#4572A7")
  
  if (is_categorical) {
    p <- ggplot(data = d, aes_string(x = "pc_1", y = "pc_2", color = groupColor)) +
      geom_point(size = 3) +
      xlab(paste0("PC", pcs[1], ": ", round(percentVar[pcs[1]] * 100), "% variance")) + theme_bw() + scale_color_manual(values = colours)
  } else {
    p <- ggplot(data = d, aes_string(x = "pc_1", y = "pc_2", fill =  groupColor)) +
      geom_point(size = 3, shape = 21, color = "grey") +
      xlab(paste0("PC", pcs[1], ": ", round(percentVar[pcs[1]] * 100), "% variance")) + theme_bw() + scale_fill_gradient2(low = vivid_colours[2], high = vivid_colours[1], mid = "white", midpoint = mean(d$group))
  }
  
  p <- p + theme(panel.border = element_blank(),
                 axis.line.x = element_line(color = "black", size = 0.2),
                 axis.line.y = element_line(color = "black", size = 0.2)) +
    ylab(paste0("PC", pcs[2], ": ", round(percentVar[pcs[2]] * 100), "% variance")) +
    coord_fixed()
  
  if (!is.null(annotate_points)) {
    p <- p + ggrepel::geom_text_repel(aes_string(label = "name"), 
                                      data = subset(d, name %in% annotate_points), 
                                      size = 3.5, 
                                      box.padding = 1.5,
                                      point.padding = 0.7,
                                      color = "black",
                                      max.overlaps = Inf,
                                      segment.color = "black")
  }
  
  return(p)
}


# twplot_pca() -----------------------------------------------------------------
# Lightweight PCA plot variant. Like plot_pca() but without scaling, without
# the categorical/continuous colour branching, and without ggrepel annotation
# support. Intended for quick exploratory two-way PCA plots where the grouping
# variable is always categorical.
#
# Side effects: assigns `pca` and `percentVar` to the global environment.
#
# Args:
#   object     : numeric matrix (rows = genes, cols = samples)
#   metadata   : data frame with sample-level metadata
#   intgroup   : character — metadata column for colouring (default "condition")
#   ntop       : number of most-variable genes to use (default 500)
#   title      : plot title string (default NA)
#   pcs        : integer vector of length 2 (default c(1,2))
#   returnData : if TRUE, return the data frame instead of plotting
#
# Returns: ggplot object, or a data frame (if returnData = TRUE)
twplot_pca <- function (object, metadata, intgroup = "condition", ntop = 500, title = NA, pcs = c(1,2),
                        returnData = FALSE) {
  rv <- rowVars(object, useNames = TRUE)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(object[select, ]), scale. = F)
  assign("pca", pca, envir = globalenv())
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  assign("percentVar", percentVar, envir = globalenv())
  if (!all(intgroup %in% colnames(metadata))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- metadata[, intgroup]
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    metadata[[intgroup]]
  }
  d <- data.frame(pc_1 = pca$x[, pcs[1]], pc_2 = pca$x[, pcs[2]],
                  group = group, intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcs[1]:pcs[2]]
    return(d)
  }
  ggplot(data = d, aes_string(x = "pc_1", y = "pc_2", color = "group")) +
    geom_point(size = 3) + xlab(paste0("PC", pcs[1], ": ", round(percentVar[pcs[1]] * 100), "% variance")) +
    ylab(paste0("PC", pcs[2], ": ", round(percentVar[pcs[2]] * 100), "% variance")) + coord_fixed() +
    ggtitle(title)
}


# plot_scree() -----------------------------------------------------------------
# Draws a scree plot (bar chart of per-PC variance + cumulative variance line)
# from a prcomp object. Designed to be called after plot_pca() using the `pca`
# object assigned to the global environment.
#
# Args:
#   pca        : prcomp object
#   components : integer vector of PC indices to include (default: all)
#   ...        : additional arguments (currently unused)
#
# Returns: a ggplot object
plot_scree <- function(pca, components = 1:length(pca$sdev), ...) {
  
  colour_vector <- c("#77AADD", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00", "#EEDD88", "#EE8866", "#FFAABB", "#DDDDDD") 
  percentVar <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
  
  plot_data <- data.frame(PC = components, Variance = percentVar[components])
  plot_data$CumulativeVariance <- cumsum(plot_data$Variance)
  
  ggplot(plot_data, aes(x = factor(PC))) +
    geom_bar(aes(y = Variance), stat = "identity", fill = colour_vector[2], alpha = 0.8) +
    geom_line(aes(y = CumulativeVariance, group = 1), color = colour_vector[1]) +
    geom_point(aes(y = CumulativeVariance), color = colour_vector[1], size = 2) +
    xlab("Principal Component") +
    ylab("Explained Variance (%)") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2))
}


# plot_pca_loadings() ----------------------------------------------------------
# Visualises PCA gene loadings for one or more components as a scatter plot
# with ggrepel labels. The top N positive and top N negative loading genes are
# shown for each requested component. Points are filled by loading magnitude
# using a blue-white-red gradient.
#
# Args:
#   pca_obj        : prcomp object
#   components     : integer vector of components to plot (default c(1,2))
#   n_loadings     : number of top positive and negative genes to show per PC (default 10)
#   col            : length-3 colour vector for gradient [low, mid, high] (default blue-white-red)
#   shape          : ggplot point shape (default 21 = filled circle with border)
#   shapeSizeRange : size range for points (default c(10,10) = constant size)
#   xlim / ylim    : optional axis limits
#   labSize        : ggrepel label size (default 3)
#   drawConnectors : logical — whether to draw ggrepel connector lines (default TRUE)
#   title / xlab / ylab : plot labels
#   legendPosition : legend placement (default "right")
#   borderWidth / borderColour : point border aesthetics
#   returnPlot     : if TRUE, return the plot object; if FALSE, print it (default TRUE)
#
# Returns: ggplot object (or NULL if returnPlot = FALSE)
plot_pca_loadings <- function(pca_obj, components = c(1, 2), n_loadings = 10,
                              col = c( "#4572A7", "white", "#AA4643"), shape = 21, shapeSizeRange = c(10, 10),
                              xlim = NULL, ylim = NULL, labSize = 3, drawConnectors = TRUE, 
                              title = "", xlab = "Principal Component", ylab = "Component Loading",
                              legendPosition = "right", borderWidth = 0.8, borderColour = "black",
                              returnPlot = TRUE) {
  
  loadings <- pca_obj$rotation[, components, drop = FALSE]
  plot_data <- data.frame(Variable = character(), PC = character(), Loading = numeric())
  
  for (i in seq_along(components)) {
    top_positive_indices <- order(loadings[, i], decreasing = TRUE)[1:n_loadings]
    top_negative_indices <- order(loadings[, i], decreasing = FALSE)[1:n_loadings]
    combined_indices <- unique(c(top_positive_indices, top_negative_indices))
    retained_loadings <- loadings[combined_indices, i, drop = FALSE]
    retained_data <- data.frame(Variable = rownames(retained_loadings), 
                                PC = paste0("PC", components[i]), 
                                Loading = retained_loadings[, 1])
    plot_data <- rbind(plot_data, retained_data)
  }
  
  plot <- ggplot(plot_data, aes(x = PC, y = Loading, fill = Loading)) +
    geom_point(shape = shape) +
    labs(title = title, x = xlab, y = ylab) +
    scale_size(range = shapeSizeRange) +
    scale_fill_gradient2(low = col[1], mid = col[2], high = col[3], midpoint = 0) +
    theme_bw() +
    theme(legend.position = legendPosition,
          panel.border = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.2),
          axis.line.y = element_line(color = "black", size = 0.2))
  
  if (!is.null(xlim)) plot <- plot + xlim(xlim)
  if (!is.null(ylim)) plot <- plot + ylim(ylim)
  
  plot <- plot + ggrepel::geom_text_repel(aes(label = Variable), 
                                          size = labSize,
                                          box.padding = unit(0.35, "lines"),
                                          point.padding = unit(0.3, "lines"),
                                          segment.color = "grey50",
                                          segment.size = 0.5,
                                          max.overlaps = Inf,
                                          force = 2)
  
  if (returnPlot) return(plot) else print(plot)
}


# compute_correlation() --------------------------------------------------------
# Computes Pearson or Spearman correlation between each of the first `num_pcs`
# PCs and each variable in `pc_metadata`. Returns a tidy long data frame ready
# for heatmap plotting, with BH-adjusted p-values and a significance flag.
#
# Args:
#   pc_eigengenes      : data frame of PC scores (rows = samples, cols = PCs)
#   pc_metadata        : data frame of numeric metadata variables
#   num_pcs            : integer — how many PCs to test
#   correlation_method : "spearman" or "pearson"
#
# Returns: long tibble with columns: pc, var, estimate, p.value, p.adjusted, significant
compute_correlation <- function(pc_eigengenes, pc_metadata, num_pcs, correlation_method) {
  if (correlation_method != "spearman" && correlation_method != "pearson") {
    stop("Invalid correlation method. Please use 'spearman' or 'pearson'.")
  }
  
  compute_cor <- function(pc, metadata_var) {
    cor_test <- cor.test(pc, metadata_var, method = correlation_method)
    tidy(cor_test)
  }
  
  first_pcs <- pc_eigengenes[, 1:num_pcs]
  
  correlation_results <- map_df(1:num_pcs, function(pc_index) {
    map_df(names(pc_metadata), function(var_name) {
      cor_data <- compute_cor(first_pcs[[pc_index]], pc_metadata[[var_name]])
      tibble(pc = paste("PC", pc_index, sep=""),
             var = var_name,
             estimate = cor_data$estimate,
             p.value = cor_data$p.value)
    })
  })
  
  correlation_results <- correlation_results %>%
    mutate(p.adjusted = p.adjust(p.value, method = "BH"))
  
  heatmap_data <- correlation_results %>%
    pivot_wider(names_from = var, values_from = estimate, id_cols = pc) %>%
    pivot_longer(cols = -pc, names_to = "var", values_to = "estimate") %>%
    left_join(correlation_results %>%
                pivot_wider(names_from = var, values_from = p.adjusted, id_cols = pc) %>%
                pivot_longer(cols = -pc, names_to = "var", values_to = "p.adjusted"),
              by = c("pc", "var")) %>%
    mutate(significant = ifelse(p.adjusted < 0.05, "bold", "normal"),
           pc = factor(pc, levels = rev(unique(pc))))
  
  return(heatmap_data)
}


# compute_kruskal_wallis() -----------------------------------------------------
# Computes Kruskal-Wallis tests between each of the first `num_pcs` PCs and
# each categorical variable in `pc_metadata`. Returns a tidy long data frame
# analogous to the output of compute_correlation(), with chi-square statistics
# and BH-adjusted p-values.
#
# Args:
#   pc_eigengenes : data frame of PC scores (rows = samples, cols = PCs)
#   pc_metadata   : data frame of categorical metadata variables
#   num_pcs       : integer — how many PCs to test
#
# Returns: long tibble with columns: pc, var, chi.square, p.value, p.adjusted, significant
compute_kruskal_wallis <- function(pc_eigengenes, pc_metadata, num_pcs) {
  compute_kw <- function(pc, metadata_var) {
    my_data <- data.frame(metric = pc, group = metadata_var)
    kw_test <- kruskal.test(metric ~ group, data = my_data)
    tidy(kw_test)
  }
  
  first_pcs <- pc_eigengenes[, 1:num_pcs]
  
  test_results <- map_df(1:num_pcs, function(pc_index) {
    map_df(names(pc_metadata), function(var_name) {
      kw_data <- compute_kw(first_pcs[[pc_index]], pc_metadata[[var_name]])
      tibble(pc = paste("PC", pc_index, sep=""),
             var = var_name,
             chi.square = kw_data$statistic,
             p.value = kw_data$p.value)
    })
  })
  
  test_results <- test_results %>%
    mutate(p.adjusted = p.adjust(p.value, method = "BH"))
  
  heatmap_data <- test_results %>%
    pivot_wider(names_from = var, values_from = chi.square, id_cols = pc) %>%
    pivot_longer(cols = -pc, names_to = "var", values_to = "chi.square") %>%
    left_join(test_results %>%
                pivot_wider(names_from = var, values_from = p.adjusted, id_cols = pc) %>%
                pivot_longer(cols = -pc, names_to = "var", values_to = "p.adjusted"),
              by = c("pc", "var")) %>%
    mutate(significant = ifelse(p.adjusted < 0.05, "bold", "normal"),
           pc = factor(pc, levels = rev(unique(pc))))
  
  return(heatmap_data)
}


# =============================================================================
# 5. WGCNA
# =============================================================================

# Kmodule_sparse_progress() ----------------------------------------------------
# Iterative module refinement using sparse matrix operations. Starting from an
# initial module assignment (e.g. from dynamicTreeCut), each gene is reassigned
# to the module with the highest mean adjacency connectivity, using sparse
# matrix multiplication for efficiency. Iterations continue until module
# assignments stabilise or MaxIterations is reached. A progress bar is displayed.
#
# Args:
#   dynamicMods   : integer vector of initial module assignments (from cutreeDynamic)
#   adjacency     : numeric adjacency matrix (nGenes × nGenes)
#   MaxIterations : maximum number of refinement iterations
#
# Returns: character vector of final module colour assignments (via labels2colors)
Kmodule_sparse_progress <- function(dynamicMods, adjacency, MaxIterations){
  adjacency <- as(adjacency, "dgCMatrix")
  nGenes <- nrow(adjacency)
  rownames(adjacency) <- 1:nGenes
  colnames(adjacency) <- 1:nGenes
  
  nIteration <- 0
  continue.change <- TRUE 
  cluster.vector <- dynamicMods
  initialClusterColors <- clusterColors <- labels2colors(cluster.vector)
  
  pb <- progress_bar$new(
    format = "Iteration :current/:total [:bar] :percent Elapsed: :elapsedfull ETA: :eta",
    total = MaxIterations, clear = FALSE, width = 60
  )
  
  while(continue.change){
    nIteration <- nIteration + 1
    clusterColors_prev <- clusterColors
    colorlevels <- unique(clusterColors)
    nCluster <- length(colorlevels)
    
    membership <- sparseMatrix(
      i = 1:nGenes,
      j = match(clusterColors, colorlevels),
      x = 1,
      dims = c(nGenes, nCluster),
      dimnames = list(NULL, colorlevels)
    )
    module_sizes <- colSums(membership)
    
    error.matrix <- adjacency %*% membership
    error.matrix <- sweep(error.matrix, 2, module_sizes, FUN = "/")
    
    max_indices <- max.col(as.matrix(error.matrix), ties.method = "first")
    clusterColors <- colorlevels[max_indices]
    
    nChange <- sum(clusterColors != initialClusterColors)
    nLastChange <- sum(clusterColors != clusterColors_prev)
    continue.change <- nIteration < MaxIterations && nLastChange > 0 
    
    pb$tick()
  }
  
  dynamicColors <- clusterColors
  return(dynamicColors)
}


# restructureModules() ---------------------------------------------------------
# Filters module gene sets to retain only genes with statistically significant
# module membership (after Bonferroni correction), then recomputes module
# eigengenes on the filtered gene sets.
#
# Args:
#   all_data          : data frame with columns Module, ENSEMBL, Module_membership,
#                       and Memership_Pvalue (note: spelling as in original pipeline)
#   pvalueThreshold   : adjusted p-value cutoff for module membership (default 0.05)
#   datExpr           : expression matrix (rows = samples, cols = genes; Ensembl IDs)
#
# Returns: list with elements:
#   $MEs          — data frame of recomputed module eigengenes (rows = samples)
#   $module_genes — named list of per-module data frames (significant genes only)
restructureModules <- function(all_data, pvalueThreshold = 0.05, datExpr = datExpr){
  
  newMEs <- data.frame(row.names = rownames(datExpr))
  moduleList <- unique(all_data$Module)
  module_genes <- list()
  
  for(module in moduleList){
    sigGenes <- all_data %>% 
      filter(Module == module) %>%
      filter(Module_membership > 0) %>%
      add_column(MM_adjusted_p = p.adjust(.$Memership_Pvalue, method = "bonferroni")) %>%
      filter(MM_adjusted_p < pvalueThreshold)
    
    module_genes[[module]] <- sigGenes
    
    geneExpr <- datExpr[, sigGenes$ENSEMBL]
    moduleEigengene <- moduleEigengenes(geneExpr, colors = rep(module, nrow(sigGenes)))$eigengenes 
    names(moduleEigengene) <- module
    newMEs <- cbind(newMEs, moduleEigengene)
  }
  
  return(list(MEs = newMEs, module_genes = module_genes))
}


# create_quantile_ranked_heatmap() ---------------------------------------------
# Draws a pheatmap in which gene expression values are replaced by their
# within-gene quantile rank (ecdf), so all genes are scaled to [0,1] and
# are directly comparable across samples regardless of absolute expression
# magnitude. Also returns per-sample CV values.
#
# Args:
#   data_matrix   : numeric matrix (rows = genes, cols = samples)
#   sample_vector : character vector of sample IDs to subset from data_matrix
#   colors        : colour palette for pheatmap (default: blue-white-red, 100 levels)
#   main_title    : string title for the heatmap
#   ...           : additional arguments passed to pheatmap()
#
# Returns: list with elements:
#   [[1]] — data frame of per-sample CV values
#   [[2]] — pheatmap object
create_quantile_ranked_heatmap <- function(data_matrix, sample_vector, colors = colorRampPalette(c("#4572A7", "white", "#AA4643"))(100), main_title = "", ...){
  
  subset_matrix <- data_matrix[, colnames(data_matrix) %in% sample_vector]
  ranked_matrix <- t(apply(subset_matrix, 1, function(x) ecdf(x)(x)))
  cv_values <- apply(subset_matrix, 2, function(x) {
    if (median(x) == 0) return(NA)
    return(sd(x) / median(x))
  })
  cv_df <- data.frame(name_id = colnames(subset_matrix), cv = cv_values)
  
  if (is.null(colors)) {
    colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
  }
  
  pheatmap_object <- pheatmap(ranked_matrix, color = colors, main = main_title, ...)
  return(list(cv_df, pheatmap_object))
}


# ── WGCNA plotting helpers ────────────────────────────────────────────────────

# .wgcna_theme() ---------------------------------------------------------------
# Internal ggplot2 theme helper shared by the WGCNA plotting functions.
# Not intended to be called directly by analysis scripts.
.wgcna_theme <- function() {
  theme_bw(base_family = "sans", base_size = 14) +
    theme(
      panel.border     = element_blank(),
      axis.line.x      = element_line(color = "black", linewidth = 0.2),
      axis.line.y      = element_line(color = "black", linewidth = 0.2),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0, family = "sans"),
      axis.title       = element_text(size = 12, family = "sans"),
      axis.text        = element_text(size = 11, family = "sans"),
      legend.position  = "none"
    )
}

.muted_red  <- "#AA4643"
.muted_blue <- "#4572A7"


# plot_soft_threshold() --------------------------------------------------------
# Replaces the base-R two-panel WGCNA soft-thresholding diagnostic. Draws
# scale-free topology R² (panel A) and mean connectivity (panel B) as a
# side-by-side patchwork ggplot.
#
# Args:
#   sft            : output of WGCNA::pickSoftThreshold()
#   powers         : the power vector tested (default 1:22)
#   r2_cutoff      : reference line on panel A (default 0.90)
#   selected_power : if not NULL, adds a vertical dotted line at this power
#   title_A / title_B : panel labels (default "A" / "B")
#
# Returns: patchwork object (two side-by-side ggplots)
plot_soft_threshold <- function(sft,
                                powers         = 1:22,
                                r2_cutoff      = 0.90,
                                selected_power = NULL,
                                title_A        = "A",
                                title_B        = "B") {
  
  df <- data.frame(
    power        = sft$fitIndices[, 1],
    signed_r2    = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    mean_connect = sft$fitIndices[, 5]
  )
  
  pA <- ggplot(df, aes(x = power, y = signed_r2)) +
    geom_hline(yintercept = r2_cutoff,
               color = "black", linewidth = 0.2, linetype = "solid") +
    annotate("text",
             x     = min(powers) + 0.2,
             y     = r2_cutoff + 0.03,
             label = paste0("R² = ", r2_cutoff),
             hjust = 0, vjust = 0, size = 2.5, color = "black") +
    geom_point(size = 1.5, color = .muted_blue, alpha = 0.8) +
    geom_text(aes(label = power), vjust = -0.8, size = 3.5, color = "grey30") +
    labs(
      title = title_A,
      x     = "Soft Threshold (power)",
      y     = "Scale Free Topology Model Fit, signed R²"
    ) +
    scale_x_continuous(breaks = powers) +
    .wgcna_theme()
  
  if (!is.null(selected_power)) {
    pA <- pA +
      geom_vline(xintercept = selected_power,
                 color = .muted_red, linewidth = 0.6, linetype = "dotted")
  }
  
  pB <- ggplot(df, aes(x = power, y = mean_connect)) +
    geom_point(size = 1.5, color = .muted_blue, alpha = 0.8) +
    geom_text(aes(label = power),
              vjust = -0.8, size = 2.8, color = "grey30") +
    labs(
      title = title_B,
      x     = "Soft Threshold (power)",
      y     = "Mean Connectivity"
    ) +
    scale_x_continuous(breaks = powers) +
    .wgcna_theme()
  
  if (!is.null(selected_power)) {
    pB <- pB +
      geom_vline(xintercept = selected_power,
                 color = .muted_red, linewidth = 0.6, linetype = "dotted")
  }
  
  pA + pB
}


# plot_wgcna_dendrogram() ------------------------------------------------------
# Replaces plotDendroAndColors(). Draws the gene dendrogram (using ggdendro)
# with one or more colour-bar tracks below it, stacked via patchwork.
#
# Args:
#   geneTree     : hclust object
#   color_mat    : a colour vector OR matrix/data.frame (one column per track)
#   track_labels : character vector of track names (one per column)
#   dendro_title : panel letter / title shown above the dendrogram (default "C")
#   bar_height   : relative height of each colour bar track (default 0.06)
#
# Returns: patchwork object (dendrogram + colour bars stacked vertically)
plot_wgcna_dendrogram <- function(geneTree,
                                  color_mat,
                                  track_labels = NULL,
                                  dendro_title = "C",
                                  bar_height   = 0.06) {
  
  if (!requireNamespace("ggdendro", quietly = TRUE))
    stop("Please install ggdendro: install.packages('ggdendro')")
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Please install patchwork: install.packages('patchwork')")
  
  ddata   <- ggdendro::dendro_data(geneTree, type = "rectangle")
  seg_df  <- ggdendro::segment(ddata)
  n_genes <- nrow(ddata$labels)
  
  branch_df <- seg_df[seg_df$yend > 0, ]
  
  y_min <- min(branch_df$yend) * 0.98
  y_max <- max(branch_df$y)
  
  p_dendro <- ggplot(branch_df) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.2, color = "grey20") +
    scale_x_continuous(expand = expansion(mult = 0.005)) +
    coord_cartesian(ylim = c(y_min, y_max * 1.05)) +
    labs(title = dendro_title, x = NULL, y = "Height") +
    .wgcna_theme() +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank(),
      plot.margin  = margin(t = 4, r = 4, b = 0, l = 4)
    )
  
  if (is.vector(color_mat) && !is.list(color_mat)) {
    color_mat <- matrix(color_mat, ncol = 1)
    if (is.null(track_labels)) track_labels <- "Module"
  }
  
  color_mat <- as.matrix(color_mat)
  n_tracks  <- ncol(color_mat)
  
  if (is.null(track_labels))
    track_labels <- if (!is.null(colnames(color_mat))) colnames(color_mat) else paste0("Track ", seq_len(n_tracks))
  
  gene_order <- order.dendrogram(as.dendrogram(geneTree))
  
  bar_plots <- lapply(seq_len(n_tracks), function(i) {
    
    colors_ordered <- color_mat[gene_order, i]
    
    tile_df <- data.frame(
      x     = seq_along(colors_ordered),
      y     = 0.5,
      color = colors_ordered
    )
    
    ggplot(tile_df, aes(x = x, y = y, fill = color)) +
      geom_tile(height = 1, width = 1) +
      scale_fill_identity() +
      scale_x_continuous(expand = expansion(mult = 0.005)) +
      scale_y_continuous(expand = c(0, 0)) +
      annotate("text",
               x = -n_genes * 0.01, y = 0.5,
               label = track_labels[i],
               hjust = 1, size = 2.8, color = "grey20") +
      labs(x = NULL, y = NULL) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(
        plot.margin = margin(t = 0, r = 4, b = 1, l = 80),
        axis.text   = element_blank(),
        axis.ticks  = element_blank()
      )
  })
  
  bar_rel <- rep(bar_height, n_tracks)
  heights <- c(1 - sum(bar_rel), bar_rel)
  
  patchwork::wrap_plots(
    c(list(p_dendro), bar_plots),
    ncol    = 1,
    heights = heights
  )
}


# =============================================================================
# 6. ENRICHMENT ANALYSIS
# =============================================================================

# custom_enricher() ------------------------------------------------------------
# Flexible over-representation analysis (ORA) using the hypergeometric test.
# Supports both standard (all-DEG) and directional (up/down split) enrichment.
# When gene_up and gene_down are provided, an UpDownRatio column is appended
# showing the count of up- and down-regulated genes per term.
#
# Args:
#   gene              : character vector of DEG Ensembl IDs (combined up+down)
#   universe          : character vector of all background gene IDs
#   geneSets          : data frame (col1 = term ID, col2 = gene ID) or named list
#   TERM2NAME         : optional named vector mapping term IDs to descriptions
#   pvalueCutoff      : raw p-value filter (default 0.05)
#   pAdjustMethod     : p.adjust method (default "BH")
#   minGSSize         : minimum overlap with DEGs (default 10)
#   maxGSSize         : maximum gene set size in background (default 500)
#   qvalueCutoff      : q-value filter applied to p.adjust column (default 0.2)
#   minGSSize_background : minimum gene set size in background (default 10)
#   gene_up / gene_down  : optional directional DEG vectors for UpDownRatio
#
# Returns: data frame of enriched terms, or NULL if none pass filters
custom_enricher <- function(gene,
                            universe,
                            geneSets,
                            TERM2NAME = NULL,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize = 10,
                            maxGSSize = 500,
                            qvalueCutoff = 0.2,
                            minGSSize_background = 10,
                            gene_up = NULL,
                            gene_down = NULL) {
  gene <- as.character(unique(gene))
  universe <- as.character(unique(universe))
  gene <- intersect(gene, universe)
  
  if (length(gene) == 0) {
    stop("None of the input genes are in the universe.")
  }
  
  analysis_type <- if (is.null(gene_up) || is.null(gene_down)) {
    "Standard Hypergeometric Enrichment"
  } else {
    "Directional Enrichment (Up/Down)"
  }
  
  if (is.data.frame(geneSets)) {
    if (ncol(geneSets) < 2) {
      stop("geneSets data frame must have at least two columns: term IDs and gene IDs.")
    }
    termIDs <- as.character(geneSets[[1]])
    geneIDs <- as.character(geneSets[[2]])
    geneSets <- split(geneIDs, termIDs)
  }
  
  geneSets <- lapply(geneSets, function(gs) intersect(gs, universe))
  term_ids <- names(geneSets)
  M <- sapply(geneSets, function(gs) length(intersect(gs, universe)))
  
  valid_terms_background <- M >= minGSSize_background & M <= maxGSSize
  geneSets <- geneSets[valid_terms_background]
  term_ids <- term_ids[valid_terms_background]
  M <- M[valid_terms_background]
  
  if (length(geneSets) == 0) {
    message("No gene sets have size between ", minGSSize_background, " and ", maxGSSize, " after filtering based on background.")
    return(NULL)
  }
  
  overlaps <- lapply(geneSets, function(gs) intersect(gs, gene))
  k <- sapply(overlaps, length)
  
  if (is.null(gene_up) || is.null(gene_down)) {
    valid_terms <- k >= minGSSize
    geneSets <- geneSets[valid_terms]
    overlaps <- overlaps[valid_terms]
    k <- k[valid_terms]
    term_ids <- term_ids[valid_terms]
    M <- M[valid_terms]
    
    if (length(term_ids) == 0) {
      message("No gene sets have at least ", minGSSize, " overlapping genes with DEGs after filtering.")
      return(NULL)
    }
    
    valid_terms_nonzero <- k > 0
    geneSets <- geneSets[valid_terms_nonzero]
    overlaps <- overlaps[valid_terms_nonzero]
    k <- k[valid_terms_nonzero]
    term_ids <- term_ids[valid_terms_nonzero]
    M <- M[valid_terms_nonzero]
    
    if (length(term_ids) == 0) {
      message("No enriched terms found after removing terms with zero overlap...")
      return(NULL)
    }
    
    N <- length(universe)
    n <- length(gene)
    pvalues <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    GeneRatio <- paste(k, n, sep = "/")
    BgRatio <- paste(M, N, sep = "/")
    RichFactor <- k / M
    FoldEnrichment <- (k * N) / (n * M)
    mu <- (M * n) / N
    sigma <- sqrt(mu * (N - n) * (N - M) / (N * (N - 1)))
    zScore <- (k - mu) / sigma
    geneID <- sapply(overlaps, function(genes) paste(genes, collapse = "/"))
    Count <- k
    
    result_df <- data.frame(
      ID = term_ids,
      Description = if (!is.null(TERM2NAME)) TERM2NAME[term_ids] else term_ids,
      GeneRatio = GeneRatio,
      BgRatio = BgRatio,
      RichFactor = RichFactor,
      FoldEnrichment = FoldEnrichment,
      zScore = zScore,
      pvalue = pvalues,
      geneID = geneID,
      Count = Count,
      stringsAsFactors = FALSE
    )
  } else {
    gene_up <- as.character(unique(gene_up))
    gene_down <- as.character(unique(gene_down))
    gene_up <- intersect(gene_up, universe)
    gene_down <- intersect(gene_down, universe)
    
    overlaps_up <- lapply(geneSets, function(gs) intersect(gs, gene_up))
    k_up <- sapply(overlaps_up, length)
    overlaps_down <- lapply(geneSets, function(gs) intersect(gs, gene_down))
    k_down <- sapply(overlaps_down, length)
    
    valid_terms <- k >= minGSSize
    geneSets <- geneSets[valid_terms]
    overlaps <- overlaps[valid_terms]
    overlaps_up <- overlaps_up[valid_terms]
    overlaps_down <- overlaps_down[valid_terms]
    k <- k[valid_terms]
    k_up <- k_up[valid_terms]
    k_down <- k_down[valid_terms]
    term_ids <- term_ids[valid_terms]
    M <- M[valid_terms]
    
    if (length(term_ids) == 0) {
      message("No gene sets have at least ", minGSSize, " overlapping genes with DEGs after filtering.")
      return(NULL)
    }
    
    valid_terms_nonzero <- k > 0
    geneSets <- geneSets[valid_terms_nonzero]
    overlaps <- overlaps[valid_terms_nonzero]
    overlaps_up <- overlaps_up[valid_terms_nonzero]
    overlaps_down <- overlaps_down[valid_terms_nonzero]
    k <- k[valid_terms_nonzero]
    k_up <- k_up[valid_terms_nonzero]
    k_down <- k_down[valid_terms_nonzero]
    term_ids <- term_ids[valid_terms_nonzero]
    M <- M[valid_terms_nonzero]
    
    if (length(term_ids) == 0) {
      message("No enriched terms found after removing terms with zero overlap...")
      return(NULL)
    }
    
    N <- length(universe)
    n <- length(gene)
    pvalues <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    GeneRatio <- paste(k, n, sep = "/")
    BgRatio <- paste(M, N, sep = "/")
    RichFactor <- k / M
    FoldEnrichment <- (k * N) / (n * M)
    mu <- (M * n) / N
    sigma <- sqrt(mu * (N - n) * (N - M) / (N * (N - 1)))
    zScore <- (k - mu) / sigma
    geneID <- sapply(overlaps, function(genes) paste(genes, collapse = "/"))
    Count <- k
    UpDownRatio <- paste(k_up, k_down, sep = ":")
    
    result_df <- data.frame(
      ID = term_ids,
      Description = if (!is.null(TERM2NAME)) TERM2NAME[term_ids] else term_ids,
      GeneRatio = GeneRatio,
      BgRatio = BgRatio,
      RichFactor = RichFactor,
      FoldEnrichment = FoldEnrichment,
      zScore = zScore,
      pvalue = pvalues,
      geneID = geneID,
      Count = Count,
      UpDownRatio = UpDownRatio,
      stringsAsFactors = FALSE
    )
  }
  
  result_df$p.adjust <- p.adjust(result_df$pvalue, method = pAdjustMethod)
  
  if (requireNamespace("qvalue", quietly = TRUE)) {
    qobj <- tryCatch(
      qvalue::qvalue(p = result_df$pvalue, lambda = 0.05, pi0.method = "bootstrap"),
      error = function(e) NULL
    )
    if (inherits(qobj, "qvalue")) {
      result_df$qvalue <- qobj$qvalues
    } else {
      result_df$qvalue <- NA
    }
  } else {
    result_df$qvalue <- NA
  }
  
  if (!"qvalue" %in% names(result_df)) {
    result_df$qvalue <- NA
  }
  
  result_df <- result_df[result_df$pvalue <= pvalueCutoff, ]
  result_df <- result_df[result_df$p.adjust <= qvalueCutoff, ]
  
  if (nrow(result_df) == 0) {
    message("No enriched terms passed the p-value and q-value cutoffs.")
    return(NULL)
  }
  
  result_df <- result_df[order(result_df$pvalue), ]
  
  result_columns <- c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
                      "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue",
                      "geneID", "Count")
  
  if ("UpDownRatio" %in% names(result_df)) {
    result_columns <- c(result_columns, "UpDownRatio")
  }
  
  result_df <- result_df[, result_columns]
  rownames(result_df) <- NULL
  attr(result_df, "analysis_type") <- analysis_type
  
  return(result_df)
}


# add_odds_ratio() -------------------------------------------------------------
# Computes odds ratios and a combined significance score from enrichment results.
#
# Args:
#   enrichment_results : data frame with GeneRatio, BgRatio, and qvalue columns
#
# Returns: the input data frame with three new columns:
#   OddsRatio, log_qvalue (intermediate), combined_score
add_odds_ratio <- function(enrichment_results) {
  enrichment_results <- as.data.frame(enrichment_results)
  
  enrichment_results$GeneCount <- as.numeric(sapply(strsplit(as.character(enrichment_results$GeneRatio), "/"), "[[", 1))
  enrichment_results$TotalGeneCount <- as.numeric(sapply(strsplit(as.character(enrichment_results$GeneRatio), "/"), "[[", 2))
  enrichment_results$BgGeneCount <- as.numeric(sapply(strsplit(as.character(enrichment_results$BgRatio), "/"), "[[", 1))
  enrichment_results$TotalBgCount <- as.numeric(sapply(strsplit(as.character(enrichment_results$BgRatio), "/"), "[[", 2))
  
  enrichment_results$OddsRatio <- with(enrichment_results, {
    deg_in_term <- GeneCount                
    deg_not_in_term <- TotalGeneCount - GeneCount
    non_deg_in_term <- BgGeneCount - GeneCount
    non_deg_not_in_term <- TotalBgCount - TotalGeneCount - non_deg_in_term
    
    for (i in seq_along(GeneCount)) {
      if (GeneCount[i] == BgGeneCount[i]) {
        deg_in_term[i] <- deg_in_term[i] + 0.5
        deg_not_in_term[i] <- deg_not_in_term[i] + 0.5
        non_deg_in_term[i] <- non_deg_in_term[i] + 0.5
        non_deg_not_in_term[i] <- non_deg_not_in_term[i] + 0.5
      }
    }
    
    (deg_in_term / deg_not_in_term) / (non_deg_in_term / non_deg_not_in_term)
  })
  
  log_qvalue <- -log10(enrichment_results$qvalue)
  enrichment_results$combined_score <- 2 / ((1 / enrichment_results$OddsRatio) + (1 / log_qvalue))
  
  return(enrichment_results)
}


# map_terms_to_genes() ---------------------------------------------------------
# Maps enriched term IDs back to the individual expressed Ensembl gene IDs
# that drove the enrichment signal. Handles both standard ORA (enricher) and
# fgsea output via the `test` parameter.
#
# Args:
#   enrichment_results : data frame / data.table with columns 'library' and
#                        either 'ID' (ORA) or 'pathway' (fgsea)
#   data_adjust        : expression matrix; rownames are Ensembl IDs
#   library_directory  : path to directory containing per-library TSV files
#   test               : pass "fgsea" if using 'pathway' column; NULL for ORA
#
# Returns: data.table with columns Term, ensembl_gene_id, Library
map_terms_to_genes <- function(enrichment_results, data_adjust, library_directory, test = NULL) {
  
  ensembl_gene_ids <- gsub('\\..*', '', rownames(data_adjust))
  enrichment_results_dt <- as.data.table(enrichment_results)
  unique_libraries <- unique(enrichment_results_dt$library)
  
  if (!is.null(test) && test == "fgsea") {
    enrichment_results_dt <- rename(enrichment_results_dt, "ID" = "pathway")
  }
  
  term_gene_mappings_list <- vector("list", length(unique_libraries))
  names(term_gene_mappings_list) <- unique_libraries
  
  for (lib in unique_libraries) {
    library_file_path <- file.path(library_directory, paste0(lib, ".tsv"))
    library_data <- fread(library_file_path, header = TRUE, sep = "\t")
    
    if (!("ensembl_gene_id" %in% names(library_data))) {
      stop("Column 'ensembl_gene_id' not found in the library data.")
    }
    
    terms_in_library <- enrichment_results_dt[library == lib, unique(ID)]
    library_data_filtered <- library_data[
      ensembl_gene_id %in% ensembl_gene_ids & Term %in% terms_in_library,
      .(Term, ensembl_gene_id)
    ]
    library_data_filtered[, Library := lib]
    term_gene_mappings_list[[lib]] <- library_data_filtered
  }
  
  term_gene_mapping <- rbindlist(term_gene_mappings_list, use.names = TRUE)
  return(term_gene_mapping)
}


# summarize_term_gene_list() ---------------------------------------------------
# Collapses a term-gene mapping data frame into one row per (Term, Library)
# with a list-column of unique gene IDs.
#
# Args:
#   term_gene_mapping : data frame with columns Term, ensembl_gene_id, Library
#
# Returns: tibble with columns Term, Library, genes (list-column)
summarize_term_gene_list <- function(term_gene_mapping) {
  term_gene_summary <- term_gene_mapping %>%
    group_by(Term, Library) %>%
    summarise(genes = list(unique(ensembl_gene_id))) %>%
    ungroup()
  
  return(term_gene_summary)
}


# jaccard_index() --------------------------------------------------------------
# Computes the Jaccard similarity index between two sets.
#
# Args:
#   set1, set2 : vectors (gene IDs or any comparable elements)
#
# Returns: numeric scalar in [0, 1]
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}


# calculate_jaccard_index_with_expression() ------------------------------------
# Computes the Jaccard similarity index between two gene sets after filtering
# both to only genes present in `expressed_genes`.
#
# Args:
#   set1, set2       : vectors of gene IDs
#   expressed_genes  : character vector of expressed gene IDs in the dataset
#
# Returns: numeric scalar in [0, 1]
calculate_jaccard_index_with_expression <- function(set1, set2, expressed_genes) {
  filtered_set1 <- intersect(set1, expressed_genes)
  filtered_set2 <- intersect(set2, expressed_genes)
  intersection_size <- length(intersect(filtered_set1, filtered_set2))
  union_size <- length(union(filtered_set1, filtered_set2))
  return(intersection_size / union_size)
}


# run_enrichment_analysis() ----------------------------------------------------
# Batch ORA across all TSV library files in `library_directory` using
# clusterProfiler::enricher(). Combines up- and down-regulated DEGs into a
# single gene vector per library.
#
# Args:
#   degs_up           : character vector of up-regulated Ensembl gene IDs
#   degs_down         : character vector of down-regulated Ensembl gene IDs
#   background_genes  : character vector of all background gene IDs
#   library_directory : path to directory containing per-library TSV files
#
# Returns: data frame of all significant enriched terms across all libraries
run_enrichment_analysis <- function(degs_up, degs_down, background_genes, library_directory) {
  enrichment_results_all <- data.frame()
  library_files <- list.files(library_directory, full.names = TRUE)
  total_libraries <- length(library_files)
  progress_bar <- txtProgressBar(min = 0, max = total_libraries, style = 3)
  
  for (library_file in library_files) {
    library_name <- gsub(".tsv", "", basename(library_file))
    library_data <- fread(library_file, header = TRUE, sep = "\t")
    library_data <- rename(library_data, "Gene" = "ensembl_gene_id")
    
    overlapping_genes_up <- intersect(degs_up, library_data$Gene)
    overlapping_genes_down <- intersect(degs_down, library_data$Gene)
    
    if (length(overlapping_genes_up) == 0 & length(overlapping_genes_down) == 0) next
    
    print(paste("Running enrichment analysis for", library_name))
    combined_degs <- union(degs_up, degs_down)
    
    enrichment_analysis <- enricher(gene = combined_degs,
                                    universe = background_genes,
                                    TERM2GENE = library_data,
                                    qvalueCutoff = 0.05)
    
    if (nrow(enrichment_analysis) > 0) {
      enrichment_analysis <- add_odds_ratio(enrichment_analysis)
      enrichment_analysis$library <- library_name
      
      enrichment_analysis$num_genes_up <- sapply(enrichment_analysis$ID, function(term) {
        term_genes <- unique(library_data$Gene[library_data$Term == term])
        sum(term_genes %in% overlapping_genes_up)
      })
      
      enrichment_analysis$num_genes_down <- sapply(enrichment_analysis$ID, function(term) {
        term_genes <- unique(library_data$Gene[library_data$Term == term])
        sum(term_genes %in% overlapping_genes_down)
      })
      
      enrichment_results_all <- rbindlist(list(enrichment_results_all, enrichment_analysis))
    }
    
    setTxtProgressBar(progress_bar, match(library_file, library_files))
  }
  
  return(enrichment_results_all)
}


# count_gene_set_terms() -------------------------------------------------------
# Counts the total number of unique terms across all library TSV files.
#
# Args:
#   library_directory : path to directory containing per-library TSV files
#
# Returns: list with elements $total_libraries and $total_terms
count_gene_set_terms <- function(library_directory) {
  total_libraries <- 0
  total_terms <- 0
  library_files <- list.files(library_directory, full.names = TRUE)
  progress_bar <- txtProgressBar(min = 0, max = length(library_files), style = 3)
  
  for (library_file in library_files) {
    library_name <- gsub(".tsv", "", basename(library_file))
    library_data <- fread(library_file, header = TRUE, sep = "\t")
    library_data <- rename(library_data, "Gene" = "ensembl_gene_id")
    unique_terms <- length(unique(library_data$Term))
    total_libraries <- total_libraries + 1
    total_terms <- total_terms + unique_terms
    setTxtProgressBar(progress_bar, match(library_file, library_files))
  }
  
  close(progress_bar)
  return(list(total_libraries = total_libraries, total_terms = total_terms))
}


# =============================================================================
# 7. VISUALISATION UTILITIES
# =============================================================================

# modify_labels() --------------------------------------------------------------
# Strips the last underscore-delimited suffix from a vector of labels.
#
# Args:
#   labels : character vector (or factor coercible to character)
#
# Returns: character vector of cleaned labels
modify_labels <- function(labels) {
  sapply(labels, function(label) {
    parts <- strsplit(as.character(label), "_")[[1]]
    if (length(parts) > 1) {
      paste(parts[-length(parts)], collapse = "_")
    } else {
      label
    }
  })
}


# gt_table() -------------------------------------------------------------------
# Formats the top `nrow` rows of an enrichment results data frame as a styled
# gt table for HTML/RMarkdown output.
#
# Args:
#   df   : data frame with standard enrichment result columns
#   nrow : number of top rows to display (default 20)
#
# Returns: a gt table object
gt_table <- function(df, nrow = 20) {
  df %>% .[seq(nrow),] %>%
    gt() %>%
    tab_header(title = "Enrichment Analysis") %>%
    cols_label(
      Pathway = "Pathway",
      Library = "Library",
      `Gene ratio` = "Gene Ratio",
      `Bg ratio` = "Bg Ratio",
      `OddsRatio` = "Odds Ratio",
      `q value` = "q adjust",
      `Up:Down ratio` = "Up:Down ratio",
      `Significance score` = "Signif. Score",
      `Louvain cluster` = "Louvain Cluster"
    ) %>%
    fmt_number(columns = c(`Significance score`), decimals = 2) %>%
    fmt_scientific(columns = `q value`, decimals = 2) %>%
    tab_style(style = list(cell_text(weight = "bold")),
              locations = cells_column_labels(everything())) %>%
    tab_options(
      table.font.size = px(10 * 1.0),
      column_labels.font.size = px(10 * 1.0),
      data_row.padding = px(2)
    ) %>%
    cols_width(
      vars(No) ~ px(80 * 0.3),
      vars(Pathway, Library) ~ px(200 * 0.6),
      vars(`q value`) ~ px(100 * 0.6),
      everything() ~ px(80 * 0.6)
    )
}

# =============================================================================
# END OF 00_functions.R
# =============================================================================