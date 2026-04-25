#!/usr/bin/env Rscript

# ============================================================
# Clonal Differential Gene Expression Analysis
# scPhylogenomics Workflow Helper Script
# Eric Wafula | 2026
# ============================================================
# Description:
#   Subsets a Seurat object to concordance-cleaned clone cells
#   from the phylogeny-inference module and runs differential
#   expression analysis (FindAllMarkers) by clone identity,
#   both across all samples combined and within each sample
#   individually.
#
# Usage:
#   Rscript run-clonal-deg-analysis.R \
#     --seurat_object <path/to/seurat.rds> \
#     --concordant_clones <path/to/concordant.clones.tsv> \
#     --output_dir <path/to/output/dir> \
#     [--min_cells 10] \
#     [--logfc_threshold 0.25] \
#     [--min_pct 0.1] \
#     [--test_use wilcox] \
#     [--assay RNA]
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(optparse)
})

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option(c("-s", "--seurat_object"), type = "character", default = NULL,
              help = "Path to Seurat v5 object (.rds) from scPhylogenomics workflow
                      (ploidy-inference or cell-typing module)",
              metavar = "FILE"),
  make_option(c("-c", "--concordant_clones"), type = "character", default = NULL,
              help = "Path to concordant clones table (.tsv or .tsv.gz) from
                      phylogeny-inference module (columns: cell_id, clone_id)",
              metavar = "FILE"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory for DEG results",
              metavar = "DIR"),
  make_option("--min_cells", type = "integer", default = 10,
              help = "Minimum number of cells required per clone to include in
                      DEG analysis [default: %default]",
              metavar = "INT"),
  make_option("--logfc_threshold", type = "double", default = 0.25,
              help = "Minimum log2 fold-change threshold for FindAllMarkers
                      [default: %default]",
              metavar = "DOUBLE"),
  make_option("--min_pct", type = "double", default = 0.1,
              help = "Minimum fraction of cells expressing a gene in either
                      group for FindAllMarkers [default: %default]",
              metavar = "DOUBLE"),
  make_option("--test_use", type = "character", default = "wilcox",
              help = "Statistical test for FindAllMarkers. Options: wilcox,
                      bimod, t, negbinom, poisson, LR, MAST, DESeq2
                      [default: %default]",
              metavar = "character"),
  make_option("--assay", type = "character", default = "RNA",
              help = "Seurat assay to use for DEG analysis [default: %default]",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Validate required arguments ---
if (is.null(opt$seurat_object))    stop("--seurat_object is required.")
if (is.null(opt$concordant_clones)) stop("--concordant_clones is required.")
if (is.null(opt$output_dir))       stop("--output_dir is required.")
if (!file.exists(opt$seurat_object)) stop("Seurat object file not found: ", opt$seurat_object)
if (!file.exists(opt$concordant_clones)) stop("Concordant clones file not found: ", opt$concordant_clones)

# --- Create output directory ---
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper: Run FindAllMarkers and write results
# ============================================================
run_deg <- function(obj, ident_col, label, output_dir, 
                    logfc_threshold, min_pct, test_use, min_cells) {
  
  # Set identity class
  Idents(obj) <- ident_col
  
  # Check which identities have enough cells
  cell_counts <- table(Idents(obj))
  valid_idents <- names(cell_counts[cell_counts >= min_cells])
  
  if (length(valid_idents) < 2) {
    message("  [SKIP] ", label, ": fewer than 2 clones with >= ", min_cells, 
            " cells. Skipping DEG analysis.")
    return(invisible(NULL))
  }
  
  dropped <- setdiff(names(cell_counts), valid_idents)
  if (length(dropped) > 0) {
    message("  [INFO] ", label, ": dropping ", paste(dropped, collapse = ", "),
            " (< ", min_cells, " cells).")
    obj <- subset(obj, idents = valid_idents)
    Idents(obj) <- ident_col
  }
  
  message("  [INFO] ", label, ": running FindAllMarkers on ",
          length(valid_idents), " clones across ",
          ncol(obj), " cells...")
  
  markers <- tryCatch({
    FindAllMarkers(
      obj,
      only.pos        = FALSE,
      logfc.threshold = logfc_threshold,
      min.pct         = min_pct,
      test.use        = test_use,
      verbose         = FALSE
    )
  }, error = function(e) {
    message("  [ERROR] ", label, ": FindAllMarkers failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(markers) || nrow(markers) == 0) {
    message("  [INFO] ", label, ": no significant markers found.")
    return(invisible(NULL))
  }
  
  # Add useful annotation columns
  markers <- markers %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    rename(!!ident_col := cluster)
  
  # Write results
  out_file <- file.path(output_dir, paste0(label, ".deg.tsv.gz"))
  readr::write_tsv(markers, out_file)
  
  message("  [SUCCESS] ", label, ": ", nrow(markers), " DEGs written to ",
          basename(out_file))
  
  # Write per-clone summary
  summary_df <- markers %>%
    group_by(.data[[ident_col]]) %>%
    summarise(
      n_deg_total  = n(),
      n_deg_up     = sum(avg_log2FC > 0),
      n_deg_down   = sum(avg_log2FC < 0),
      n_sig_padj   = sum(p_val_adj < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_file <- file.path(output_dir, paste0(label, ".deg.summary.tsv.gz"))
  readr::write_tsv(summary_df, summary_file)
  message("  [SUCCESS] ", label, ": clone DEG summary written to ",
          basename(summary_file))
  
  return(invisible(markers))
}


# ============================================================
# MAIN
# ============================================================

# --- Step 1: Load inputs ---
message("\n========== CLONAL DEG ANALYSIS ==========")
message("\n[Step 1] Loading inputs...")

message("  [INFO] Loading Seurat object: ", basename(opt$seurat_object))
obj <- readRDS(opt$seurat_object)
message("  [INFO] Seurat object loaded: ", ncol(obj), " cells x ", nrow(obj), " genes.")

message("  [INFO] Loading concordant clones: ", basename(opt$concordant_clones))
concordant_df <- readr::read_tsv(opt$concordant_clones, show_col_types = FALSE)
concordant_df$cell_id  <- as.character(concordant_df$cell_id)
concordant_df$clone_id <- as.character(concordant_df$clone_id)
message("  [INFO] Concordant clones loaded: ", nrow(concordant_df), " cells across ",
        length(unique(concordant_df$clone_id)), " clones (",
        paste(sort(unique(concordant_df$clone_id)), collapse = ", "), ").")

# --- Step 2: Subset Seurat object to concordant cells ---
message("\n[Step 2] Subsetting Seurat object to concordant cells...")

# Validate that cell_id column exists in Seurat metadata
if (!"cell_id" %in% colnames(obj@meta.data)) {
  stop("Seurat object metadata does not contain a 'cell_id' column.")
}

common_cells <- intersect(obj@meta.data$cell_id, concordant_df$cell_id)
message("  [INFO] Cells in Seurat object:    ", ncol(obj))
message("  [INFO] Cells in concordant table: ", nrow(concordant_df))
message("  [INFO] Cells in common:           ", length(common_cells))

if (length(common_cells) < 10) {
  stop("Fewer than 10 cells overlap between Seurat object and concordant clones table. ",
       "Check that cell_id values match between inputs.")
}

# Subset to concordant cells using cell_id
obj <- subset(obj, cell_id %in% common_cells)
message("  [INFO] Seurat object subsetted to ", ncol(obj), " concordant cells.")

# --- Step 3: Add clone_id to Seurat metadata ---
message("\n[Step 3] Adding clone_id to Seurat metadata...")

# Match clone_id to cells by cell_id
clone_lookup           <- concordant_df$clone_id
names(clone_lookup)    <- concordant_df$cell_id
obj$clone_id           <- clone_lookup[obj@meta.data$cell_id]

# Validate
n_na <- sum(is.na(obj$clone_id))
if (n_na > 0) {
  message("  [WARNING] ", n_na, " cells could not be matched to a clone_id. ",
          "These will be excluded from DEG analysis.")
  obj <- subset(obj, !is.na(clone_id))
}
message("  [INFO] clone_id added to metadata. Clone breakdown:")
for (cln in sort(unique(obj$clone_id))) {
  message("      Clone ", cln, ": ", sum(obj$clone_id == cln), " cells")
}

# --- Step 4: Set assay and normalize if needed ---
message("\n[Step 4] Setting assay and preparing for DEG analysis...")
DefaultAssay(obj) <- opt$assay

# Normalize if data layer is not populated
if (!("data" %in% Layers(obj))) {
  message("  [INFO] Normalizing data layer...")
  obj <- NormalizeData(obj, verbose = FALSE)
}
message("  [INFO] Using assay: ", opt$assay)

# --- Step 5: Create sample_clone_id metadata column ---
message("\n[Step 5] Creating sample_clone_id metadata column...")

if (!"sample_id" %in% colnames(obj@meta.data)) {
  message("  [WARNING] 'sample_id' column not found in metadata. ",
          "Skipping per-sample DEG analysis.")
  run_per_sample <- FALSE
} else {
  obj$sample_clone_id <- paste(obj$sample_id, obj$clone_id, sep = "_")
  message("  [INFO] sample_clone_id created. Examples: ",
          paste(head(unique(obj$sample_clone_id), 5), collapse = ", "))
  run_per_sample <- TRUE
}

# --- Step 6: DEG analysis across ALL samples combined by clone_id ---
message("\n[Step 6] Running DEG analysis across all samples by clone_id...")

run_deg(
  obj             = obj,
  ident_col       = "clone_id",
  label           = "all_samples.by_clone",
  output_dir      = opt$output_dir,
  logfc_threshold = opt$logfc_threshold,
  min_pct         = opt$min_pct,
  test_use        = opt$test_use,
  min_cells       = opt$min_cells
)

# --- Step 7: Per-sample DEG analysis by sample_clone_id ---
if (run_per_sample) {
  sample_ids <- sort(unique(obj$sample_id))
  message("\n[Step 7] Running per-sample DEG analysis for ", 
          length(sample_ids), " samples...")
  
  # Create per-sample output subdirectory
  per_sample_dir <- file.path(opt$output_dir, "per_sample")
  dir.create(per_sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (sid in sample_ids) {
    message("\n  >>> Sample: ", sid)
    
    # Subset to this sample
    obj_sample <- subset(obj, sample_id == sid)
    message("  [INFO] ", sid, ": ", ncol(obj_sample), " cells.")
    
    # Run DEG within this sample using sample_clone_id as identity
    run_deg(
      obj             = obj_sample,
      ident_col       = "sample_clone_id",
      label           = paste0(sid, ".by_sample_clone"),
      output_dir      = per_sample_dir,
      logfc_threshold = opt$logfc_threshold,
      min_pct         = opt$min_pct,
      test_use        = opt$test_use,
      min_cells       = opt$min_cells
    )
  }
}

# --- Step 8: DEG analysis across ALL clones combined by sample_id ---
message("\n[Step 8] Running DEG analysis across all clones by sample_id...")

if (run_per_sample) {
  run_deg(
    obj             = obj,
    ident_col       = "sample_id",
    label           = "all_clones.by_sample",
    output_dir      = opt$output_dir,
    logfc_threshold = opt$logfc_threshold,
    min_pct         = opt$min_pct,
    test_use        = opt$test_use,
    min_cells       = opt$min_cells
  )
} else {
  message("  [SKIP] sample_id not available. Skipping by-sample DEG analysis.")
}

# --- Step 9: Per-clone DEG analysis by sample_clone_id ---
if (run_per_sample) {
  clone_ids <- sort(unique(obj$clone_id))
  message("\n[Step 9] Running per-clone DEG analysis for ",
          length(clone_ids), " clones...")

  # Create per-clone output subdirectory
  per_clone_dir <- file.path(opt$output_dir, "per_clone")
  dir.create(per_clone_dir, recursive = TRUE, showWarnings = FALSE)

  for (cln in clone_ids) {
    message("\n  >>> Clone: ", cln)

    # Subset to this clone
    obj_clone <- subset(obj, clone_id == cln)
    message("  [INFO] Clone ", cln, ": ", ncol(obj_clone), " cells.")

    # Run DEG within this clone using sample_clone_id as identity
    run_deg(
      obj             = obj_clone,
      ident_col       = "sample_clone_id",
      label           = paste0("clone_", cln, ".by_clone_sample"),
      output_dir      = per_clone_dir,
      logfc_threshold = opt$logfc_threshold,
      min_pct         = opt$min_pct,
      test_use        = opt$test_use,
      min_cells       = opt$min_cells
    )
  }
}

message("\n========== DEG ANALYSIS COMPLETE ==========")
message("Results written to: ", opt$output_dir)
