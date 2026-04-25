#!/usr/bin/env Rscript

# ============================================================
# Clonal GSEA Hallmark Analysis
# scPhylogenomics Workflow Helper Script
# Eric Wafula | 2026
# ============================================================
# Description:
#   Performs GSEA using MSigDB Hallmark gene sets on all DE
#   result files produced by run-clonal-deg-analysis.R.
#   For each DE file, runs GSEA per group, writes enrichment
#   tables beside the DE files, and generates comparative
#   GSEA dotplots for significant pathways across groups.
#
# Input DE result structure expected:
#   <deg_dir>/
#     all_samples.by_clone.deg.tsv        -> group: clone_id
#     all_clones.by_sample.deg.tsv        -> group: clone_id (sample_id)
#     per_sample/<SAMPLE>.deg.tsv         -> group: sample_clone_id
#     per_clone/<CLONE>.deg.tsv           -> group: sample_clone_id
#
# Usage:
#   Rscript run-clonal-gsea-analysis.R \
#     --deg_dir <path/to/deg/results> \
#     --plots_dir <path/to/output/plots>  \
#     [--min_gsea_size 10] \
#     [--max_gsea_size 500] \
#     [--nperm 10000] \
#     [--padj_threshold 0.05]
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
  library(optparse)
})

set.seed(123)

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option(c("-d", "--deg_dir"), type = "character", default = NULL,
              help = "Path to DE results directory produced by run-clonal-deg-analysis.R",
              metavar = "DIR"),
  make_option(c("-o", "--plots_dir"), type = "character", default = NULL,
              help = "Output directory for comparative GSEA dotplots",
              metavar = "DIR"),
  make_option("--min_gsea_size", type = "integer", default = 10,
              help = "Minimum gene set size for fgsea [default: %default]",
              metavar = "INT"),
  make_option("--max_gsea_size", type = "integer", default = 500,
              help = "Maximum gene set size for fgsea [default: %default]",
              metavar = "INT"),
  make_option("--nperm", type = "integer", default = 10000,
              help = "Number of permutations for fgsea [default: %default]",
              metavar = "INT"),
  make_option("--padj_threshold", type = "double", default = 0.05,
              help = "Adjusted p-value threshold for significant pathways in dotplot
                      [default: %default]",
              metavar = "DOUBLE")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Validate required arguments ---
if (is.null(opt$deg_dir))   stop("--deg_dir is required.")
if (is.null(opt$plots_dir)) stop("--plots_dir is required.")
if (!dir.exists(opt$deg_dir)) stop("DEG results directory not found: ", opt$deg_dir)

# --- Create output directories ---
dir.create(opt$plots_dir,                              recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$plots_dir, "per_sample"),     recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$plots_dir, "per_clone"),      recursive = TRUE, showWarnings = FALSE)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

# --- Build ranked gene list from DEG table for one group ---
# Uses signed -log10(p_val) * sign(avg_log2FC) as ranking metric.
# This preserves both direction and significance, which is more
# appropriate for Seurat FindAllMarkers output than raw log2FC alone.
# Given Seurat’s single‑cell DE idiosyncrasies, 
# weighting by p_val (or better, p_val_adj) is a reasonable choice
build_ranked_list <- function(deg_df) {
  deg_df %>%
    filter(!is.na(gene), !is.na(avg_log2FC), !is.na(p_val)) %>%
    mutate(
      # p_val_safe = pmax(p_val, 1e-300),   # avoid -Inf from log10(0)
      p_adj_safe = pmax(p_val_adj, 1e-300),      
      # rank_score = sign(avg_log2FC) * (-log10(p_val_safe))
      rank_score = sign(avg_log2FC) * (-log10(p_adj_safe)) + (avg_log2FC * 0.01)
    ) %>%
    arrange(desc(rank_score)) %>%
    distinct(gene, .keep_all = TRUE) %>%  # keep highest-ranked if duplicated
    pull(rank_score, name = gene)
}

# --- Run fgsea on a ranked list ---
run_gsea <- function(ranked_list, hallmark_sets, min_size, max_size, nperm) {
  if (length(ranked_list) < min_size) {
    message("    [SKIP] Fewer than ", min_size, " genes in ranked list.")
    return(NULL)
  }
  
  result <- tryCatch({
    fgsea::fgsea(
      pathways   = hallmark_sets,
      stats      = ranked_list,
      minSize    = min_size,
      maxSize    = max_size,
      nPermSimple = nperm
    ) %>%
      as_tibble() %>%
      mutate(
        pathway = gsub("HALLMARK_", "", pathway),
        pathway = gsub("_", " ", pathway)
      ) %>%
      arrange(padj)
  }, error = function(e) {
    message("    [ERROR] fgsea failed: ", e$message)
    NULL
  })
  
  return(result)
}

# --- Flatten leadingEdge list column to comma-separated string for writing ---
flatten_gsea <- function(gsea_df) {
  gsea_df %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(unlist(x), collapse = ",")))
}

# --- Process one DE file: run GSEA per group, save results, return combined ---
process_deg_file <- function(deg_file, group_col, hallmark_sets,
                              min_size, max_size, nperm, padj_threshold) {
  
  label    <- tools::file_path_sans_ext(basename(deg_file))
  out_dir  <- dirname(deg_file)
  
  message("  [INFO] Processing: ", basename(deg_file))
  message("  [INFO] Grouping by: ", group_col)
  
  # Load DEG table
  deg_df <- tryCatch({
    readr::read_tsv(deg_file, show_col_types = FALSE)
  }, error = function(e) {
    message("  [ERROR] Could not read ", basename(deg_file), ": ", e$message)
    return(NULL)
  })
  
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    message("  [SKIP] Empty or unreadable file: ", basename(deg_file))
    return(NULL)
  }
  
  # Validate required columns
  required_cols <- c("gene", "avg_log2FC", "p_val", group_col)
  missing_cols  <- setdiff(required_cols, colnames(deg_df))
  if (length(missing_cols) > 0) {
    message("  [SKIP] Missing columns in ", basename(deg_file), ": ",
            paste(missing_cols, collapse = ", "))
    return(NULL)
  }
  
  groups         <- sort(unique(deg_df[[group_col]]))
  all_gsea_list  <- list()
  
  for (grp in groups) {
    message("    >>> Group: ", grp)
    
    grp_df <- deg_df %>% filter(.data[[group_col]] == grp)
    
    if (nrow(grp_df) < min_size) {
      message("    [SKIP] Too few genes (", nrow(grp_df), ") for group: ", grp)
      next
    }
    
    # Build ranked gene list
    ranked <- build_ranked_list(grp_df)
    
    if (length(ranked) == 0) {
      message("    [SKIP] Ranked list empty for group: ", grp)
      next
    }
    
    message("    [INFO] Ranked genes: ", length(ranked))
    
    # Run GSEA
    gsea_result <- run_gsea(ranked, hallmark_sets, min_size, max_size, nperm)
    
    if (is.null(gsea_result) || nrow(gsea_result) == 0) {
      message("    [SKIP] No GSEA results for group: ", grp)
      next
    }
    
    # Add group label
    gsea_result[[group_col]] <- grp
    
    # Save per-group GSEA results beside the DE file
    gsea_out_file <- file.path(out_dir,
                               paste0(label, ".", grp, ".gsea.hallmark.tsv.gz"))
    flatten_gsea(gsea_result) %>%
      readr::write_tsv(gsea_out_file)
    
    n_sig <- sum(gsea_result$padj < padj_threshold, na.rm = TRUE)
    message("    [SUCCESS] ", nrow(gsea_result), " pathways (",
            n_sig, " significant) written to ", basename(gsea_out_file))
    
    all_gsea_list[[as.character(grp)]] <- gsea_result
  }
  
  if (length(all_gsea_list) == 0) return(NULL)
  
  # Combine all groups for dotplot
  combined <- bind_rows(all_gsea_list)
  return(combined)
}

# --- Generate and save comparative GSEA barplot ---
# Shows union of pathways significant in at least one group.
# Every group x pathway combination is shown:
#   - Colored bar with *  : significant in that group (padj < threshold)
#   - Colored bar no *    : present but not significant in that group
#   - Flat line at NES=0  : pathway not detected in that group at all
plot_gsea_dotplot <- function(combined_gsea, group_col, label, plots_out_dir,
                               padj_threshold) {
  
  if (is.null(combined_gsea) || nrow(combined_gsea) == 0) {
    message("  [SKIP] No GSEA results to plot for: ", label)
    return(invisible(NULL))
  }
  
  # --- Identify union of pathways significant in at least one group ---
  sig_pathways <- combined_gsea %>%
    filter(padj < padj_threshold & !is.na(padj)) %>%
    pull(pathway) %>%
    unique()
  
  if (length(sig_pathways) == 0) {
    message("  [SKIP] No significant pathways (padj < ", padj_threshold,
            ") to plot for: ", label)
    return(invisible(NULL))
  }
  
  # --- Build complete group x pathway grid ---
  # Missing combinations (pathway not detected in a group) get NES=0, padj=NA
  all_groups <- sort(unique(as.character(combined_gsea[[group_col]])))
  
  full_grid <- expand.grid(
    pathway     = sig_pathways,
    group_label = all_groups,
    stringsAsFactors = FALSE
  )
  
  # Join observed results onto full grid
  observed <- combined_gsea %>%
    filter(pathway %in% sig_pathways) %>%
    mutate(group_label = as.character(.data[[group_col]]))
  
  plot_df <- full_grid %>%
    left_join(observed %>% select(pathway, group_label, NES, padj),
              by = c("pathway", "group_label")) %>%
    mutate(
      NES          = replace_na(NES, 0),
      padj         = replace_na(padj, 1),
      significance = ifelse(padj < padj_threshold, "*", ""),
      group_label  = factor(group_label, levels = all_groups),
      pathway      = factor(pathway, levels = rev(sort(unique(sig_pathways))))
    )
  
  # --- Dynamic plot dimensions ---
  plot_height <- max(4, length(sig_pathways) * 0.35 + 2)
  plot_width  <- max(6, length(all_groups) * 3 + 3)
  
  # Subset significant rows for geom_text to avoid hjust length mismatch
  sig_df <- plot_df %>% filter(significance == "*")
  
  # --- Build plot ---
  p <- ggplot2::ggplot(plot_df,
                       aes(x = NES, y = pathway, fill = NES)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey30") +
    geom_text(
      data  = sig_df,
      aes(label = significance,
          x     = ifelse(NES >= 0, NES + 0.05, NES - 0.05)),
      hjust = ifelse(sig_df$NES >= 0, 0, 1),
      size  = 4, color = "black"
    ) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, name = "NES") +
    facet_wrap(~ group_label, nrow = 1) +
    theme_bw() +
    theme(
      axis.title.y       = element_blank(),
      axis.title.x       = element_text(size = 9),
      axis.text.y        = element_text(size = 8),
      axis.text.x        = element_text(size = 8),
      strip.text         = element_text(size = 9, face = "bold"),
      strip.background   = element_rect(fill = "grey92"),
      plot.title         = element_text(size = 10, face = "bold"),
      panel.grid.major.y = element_line(color = "grey90"),
      legend.position    = "right"
    ) +
    labs(
      x     = "Normalized Enrichment Score (NES)",
      title = paste("GSEA Hallmark —", label,
                    "\n(padj <", padj_threshold,
                    "| * = significant | no bar = not detected)")
    )
  
  out_plot <- file.path(plots_out_dir, paste0(label, ".gsea.hallmark.barplot.pdf"))
  ggplot2::ggsave(out_plot, p, width = plot_width, height = plot_height)
  message("  [PLOT] Barplot saved to ", basename(out_plot))
  
  return(invisible(p))
}


# ============================================================
# MAIN
# ============================================================
message("\n========== CLONAL GSEA HALLMARK ANALYSIS ==========")

# --- Step 1: Load MSigDB Hallmark gene sets ---
message("\n[Step 1] Loading MSigDB Hallmark gene sets...")
hs_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
message("  [INFO] Loaded ", length(hs_hallmark), " Hallmark gene sets.")

# --- Step 2: Discover all DE result files ---
message("\n[Step 2] Discovering DE result files in: ", opt$deg_dir)

# Map each file pattern to its group column
deg_files_map <- list(
  list(
    files     = list.files(opt$deg_dir, pattern = "all_samples\\..*\\.deg\\.tsv\\.gz$",
                           full.names = TRUE),
    group_col = "clone_id",
    plots_out = opt$plots_dir,
    tier      = "all_samples"
  ),
  list(
    files     = list.files(opt$deg_dir, pattern = "all_clones\\..*\\.deg\\.tsv\\.gz$",
                           full.names = TRUE),
    group_col = "sample_id",
    plots_out = opt$plots_dir,
    tier      = "all_clones"
  ),
  list(
    files     = list.files(file.path(opt$deg_dir, "per_sample"),
                           pattern = "\\.deg\\.tsv\\.gz$", full.names = TRUE),
    group_col = "sample_clone_id",
    plots_out = file.path(opt$plots_dir, "per_sample"),
    tier      = "per_sample"
  ),
  list(
    files     = list.files(file.path(opt$deg_dir, "per_clone"),
                           pattern = "\\.deg\\.tsv\\.gz$", full.names = TRUE),
    group_col = "sample_clone_id",
    plots_out = file.path(opt$plots_dir, "per_clone"),
    tier      = "per_clone"
  )
)

# Filter out summary files
deg_files_map <- lapply(deg_files_map, function(x) {
  x$files <- x$files[!grepl("\\.summary\\.tsv\\.gz$", x$files)]
  x
})

total_files <- sum(sapply(deg_files_map, function(x) length(x$files)))
message("  [INFO] Found ", total_files, " DE result files to process.")

# --- Step 3: Process each tier ---
message("\n[Step 3] Running GSEA and generating dotplots...")

for (tier_info in deg_files_map) {
  
  if (length(tier_info$files) == 0) {
    message("\n  [SKIP] No files found for tier: ", tier_info$tier)
    next
  }
  
  message("\n  >>> Tier: ", tier_info$tier, " (", length(tier_info$files), " files)")
  dir.create(tier_info$plots_out, recursive = TRUE, showWarnings = FALSE)
  
  for (deg_file in tier_info$files) {
    
    # Skip summary files defensively
    if (grepl("summary", basename(deg_file))) next
    
    label <- tools::file_path_sans_ext(basename(deg_file))
    # Strip .deg suffix from label for cleaner file names
    label <- sub("\\.deg$", "", label)
    
    message("\n  --- File: ", basename(deg_file), " ---")
    
    combined_gsea <- process_deg_file(
      deg_file      = deg_file,
      group_col     = tier_info$group_col,
      hallmark_sets = hs_hallmark,
      min_size      = opt$min_gsea_size,
      max_size      = opt$max_gsea_size,
      nperm         = opt$nperm,
      padj_threshold = opt$padj_threshold
    )
    
    plot_gsea_dotplot(
      combined_gsea  = combined_gsea,
      group_col      = tier_info$group_col,
      label          = label,
      plots_out_dir  = tier_info$plots_out,
      padj_threshold = opt$padj_threshold
    )
  }
}

message("\n========== GSEA ANALYSIS COMPLETE ==========")
message("GSEA tables written beside DE result files in: ", opt$deg_dir)
message("Dotplots written to: ", opt$plots_dir)
