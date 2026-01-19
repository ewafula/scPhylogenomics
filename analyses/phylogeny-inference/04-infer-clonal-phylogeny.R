#!/usr/bin/env Rscript

# ============================================================
# Plot clonal phylogenetic trees from SNP MSA
# Eric Wafula | 2025
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(phytools)
  library(tidyverse)
  library(ggplot2)
  library(adephylo)
  library(optparse)
  library(rprojroot)
})

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option(opt_str = "--project", type = "character", default = NULL,
              help = "A valid scPhylogenomics project name",
              metavar = "character"),
  make_option(opt_str = "--refseq", action = "store_true", default = FALSE,
              help = "Whether to include reference sequence [default: %default]"),
  make_option(opt_str = "--rescale", action = "store_true", default = FALSE,
              help = "Whether to rescale branch lengths to compress extremes [default: %default]"),
  make_option(opt_str = "--percentile_outlier", type = "double", default = 0.99,
              help = "Determine likely outlier super long branches to exclude [default: %default]
                      The default means the value above which the top 1% of tip distances fall",
              metavar = "percentile"),
  make_option(opt_str = "--annot_type", type = "character", default = "Clone",
              help = "Clone metadata annotation - [default: %default]
                      - Available choices - Clone, Lineage, and Sample 
                      - Sample choice only works with merged samples alignment",
              metavar = "character"),
  make_option(opt_str = "--cell_category", type = "character", default = NULL,
              help = "Optional cell annotation categories file (if not provided, skipped)",
              metavar = "character")
)

# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project <- opt$project 
refseq <- opt$refseq
rescale <- opt$rescale
percentile_outlier <- opt$percentile_outlier
annot_type <- opt$annot_type
cell_category <- opt$cell_category

stopifnot(opt$annot_type %in% c("Clone", "Lineage", "Sample"))
if (opt$annot_type == "Lineage" && is.null(opt$cell_category)) {
  stop("--cell_category required when annot_type = Lineage")
}

# Set seed for reproducibility
set.seed(123)

# ---------------------------
# Paths
# ---------------------------
# Establish base directory
root_dir   <- rprojroot::find_root(has_dir(".git"))

# Set data directories
# Set data directories
module_dir <- file.path(root_dir, "analyses", "phylogeny-inference")
results_dir <- file.path(module_dir, "results", project)
inputs_dir <- file.path(module_dir, "inputs")
plots_dir <- file.path(module_dir, "plots", project)

# Create plots directory if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# ---------------------------
# Load trees
# ---------------------------
# List all tree files
tree_files <- list.files(results_dir, pattern = "\\.tree$", full.names = TRUE)
stopifnot(length(tree_files) > 0)

# Reference sequence inferred in the snv-module
reference_sequence <- "REFERENCE"

# ============================================================
# MAIN LOOP
# ============================================================
# Process each tree file
for (tree_file in tree_files) {
  
  message("\n========== Processing tree: ", basename(tree_file), " ==========")
  
  # Load tree
  tree <- ggtree::read.tree(tree_file)
  
  ### Rooting ###
  # Check if the reference sequence exists
  if (reference_sequence %in% tree$tip.label) {
    tree <- ape::root(tree, reference_sequence, resolve.root = TRUE)
    
    # Status flag to keep or exclude root sequence
    if (!refseq) tree <- ape::drop.tip(tree, reference_sequence)
  } else {
    # Fallback: Root by Midpoint using the phytools package
    message(paste("WARNING: Reference sequence", reference_sequence, "not found in", 
                  tree_file, ". Rooted by midpoint instead."))
    tree <- phytools::midpoint.root(tree)
  }
  
  ### Remove long branches ###
  # Compute root-to-tip distances
  tip_dists <- adephylo::distRoot(tree)
  
  # Define outlier tips
  threshold <- quantile(tip_dists, percentile_outlier)
  outliers  <- names(tip_dists[tip_dists > threshold])
  
  # Drop outlier tips
  message("Removing ", length(outliers), " long-branch tips")
  tree_cleaned <- drop.tip(tree, outliers)
  
  # Get treefile tree method and  basename
  tree_method <- sub(".*\\.(fasttree|iqtree)\\.tree$", "\\1", basename(tree_file))
  base_name   <- sub("\\.(fasttree|iqtree)\\.tree$", "", basename(tree_file))
  
  # Check that at least one clone file exists
  clone_files <- c(
    hierarchical = file.path(results_dir, paste0(base_name, ".clones.hierarchical.tsv.gz")),
    manifold     = file.path(results_dir, paste0(base_name, ".clones.manifold.tsv.gz"))
  )
  clone_files <- clone_files[file.exists(clone_files)]
  stopifnot(length(clone_files) > 0)
  
  # ------------------------------------------------------------
  # Per clone clustering method
  # ------------------------------------------------------------
  for (clone_type in names(clone_files)) {
    
    message("\n  --- Processing ", clone_type, " clustering ---")
    
    tree_current <- tree_cleaned
    clone_df <- readr::read_tsv(clone_files[[clone_type]], show_col_types = FALSE)
    
    common_cells <- intersect(tree_current$tip.label, clone_df$cell_id)
    tree_current <- ape::keep.tip(tree_current, common_cells)
    
    message("  [DEBUG] tips after filtering: ", length(tree_current$tip.label))
    
    # Annotation mode
    if (annot_type == "Clone") {
      annot_id <- "clone_id"
      tree_type <- "clonal"
    } else if (annot_type == "Sample") {
      annot_id <- "sample_id"
      tree_type <- "sample"
    } else {
      annot_id <- "cell_category"
      tree_type <- "lineage"
      clone_df <- readr::read_tsv(file.path(inputs_dir, cell_category), show_col_types = FALSE)
    }
    
    # Build metadata
    tree_data <- tibble(label = tree_current$tip.label) %>%
      inner_join(clone_df, by = c("label" = "cell_id"))
    
    # force discrete annotation
    tree_data[[annot_id]] <- as.factor(tree_data[[annot_id]])
    stopifnot(is.factor(tree_data[[annot_id]]))
    
    clone_levels <- levels(tree_data[[annot_id]])
    clone_colors <- setNames(scales::hue_pal()(length(clone_levels)), clone_levels)
    
    # Optional branch rescaling
    if (rescale) {
      tree_current$edge.length <- log1p(tree_current$edge.length)
    }
    
    # Plot
    p <- ggtree(tree_current, layout = "rectangular") %<+% tree_data +
      geom_tippoint(aes(color = .data[[annot_id]]), size = 1.4) +
      scale_color_manual(values = clone_colors) +
      theme_tree2() +
      theme(
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank()
      ) +
      guides(color = guide_legend(title = opt$annot_type))
    
    p$data$label <- ifelse(p$data$isTip, "", p$data$label)
    
    # Clade labels
    for (cl in clone_levels) {
      tips <- tree_data %>% filter(.data[[annot_id]] == cl) %>% pull(label)
      if (length(tips) > 1) {
        node <- MRCA(tree_current, tips)
        p <- p + geom_cladelabel(
          node = node, label = cl,
          color = clone_colors[cl],
          fontsize = 3, align = TRUE
        )
      }
    }
    
    print(p)
    
    # Save outputs
    plot_file <- file.path(
      plots_dir,
      paste0(base_name, ".", tree_method, ".", clone_type, ".", tree_type, "_tree.png")
    )
    ggsave(plot_file, p, width = 10, height = 8, dpi = 300)
    
    tree_out <- file.path(
      results_dir,
      paste0(base_name, ".", tree_method, ".", clone_type, ".", tree_type, ".treefile")
    )
    write.tree(tree_current, tree_out)
    
    message("  Saved plot: ", basename(plot_file))
  }
}

message("\n========== ALL TREES PROCESSED SUCCESSFULLY ==========")
