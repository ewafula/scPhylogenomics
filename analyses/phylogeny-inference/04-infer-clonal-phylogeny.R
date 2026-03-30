#!/usr/bin/env Rscript

# ============================================================
# Plot Clonal Phylogenetic Trees with Concordance Cleaning
# Eric Wafula & Sayaka | 2026
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
  library(phangorn) 
})

# ---------------------------
# Command-line options
# ---------------------------
option_list <- list(
  make_option(c("-p", "--project"), type = "character", default = NULL,
              help = "A valid scPhylogenomics project name", metavar = "character"),
  make_option("--tree_method", type = "character", default = NULL,
              help = "Phylogeny method to process (iqtree or fasttree). If NULL, processes all.", 
              metavar = "character"),
  make_option("--cluster_method", type = "character", default = NULL,
              help = "Clustering method to process (manifold or hierarchical). If NULL, processes all.", 
              metavar = "character"),
  make_option("--refseq", action = "store_true", default = FALSE,
              help = "Include reference sequence [default: %default]"),
  make_option("--rescale", action = "store_true", default = FALSE,
              help = "Rescale branch lengths (log1p) [default: %default]"),
  make_option("--percentile_outlier", type = "double", default = 0.99,
              help = "Percentile for long-branch removal [default: %default]"),
  make_option("--annot_type", type = "character", default = "Clone",
              help = "Annotation type: Clone, Lineage, or Sample [default: %default]"),
  make_option("--cell_category", type = "character", default = NULL,
              help = "Optional cell annotation categories file", metavar = "character"),
  # --- CONCORDANCE CLEANING ---
  make_option("--concordance_clean", action = "store_true", default = FALSE,
              help = "Perform tree-clone concordance cleaning [default: %default]"),
  make_option("--min_clone_frac", type = "double", default = 0.5,
              help = "Minimum clone fraction for a node to be considered dominant [default: %default]"),
  make_option("--merge_clusters", type = "character", default = NULL,
              help = "Specify cluster IDs to be merged into groups. 
                      Format: '1:4,5:3:7' (merges 1 & 4; and 5, 3, & 7). 
                      New IDs are named by concatenating the original IDs.
                      [default: %default]", 
              metavar = "ID_LIST")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project)) stop("Project name is required.")

# ---------------------------
# Helper Function: Concordance Cleaning
# Following Sayaka's approach: find largest clade dominated by each clone,
# keep only concordant tips — no ace() required.
#
# Returns a list with two elements:
#   $tree          - pruned phylogenetic tree (concordant cells only)
#   $clean_clone_df - data frame with cell_id and clone_id for all cells
#                     retained in the cleaned tree (for DEG analysis)
# ---------------------------
clean_tree_by_concordance <- function(tree, clone_df, min_clone_frac = 0.5) {
  message("  [Cleaning] Running concordance-based tree cleaning...")
  
  clone_df          <- as.data.frame(clone_df)
  clone_df$cell_id  <- as.character(clone_df$cell_id)
  clone_df$clone_id <- as.character(clone_df$clone_id)
  
  # --- Drop tips with no clone assignment (e.g. REFERENCE) ---
  common_cells  <- intersect(tree$tip.label, clone_df$cell_id)
  initial_count <- length(common_cells)
  tips_to_drop  <- setdiff(tree$tip.label, common_cells)
  
  if (length(tips_to_drop) > 0) {
    message("  [INFO] Dropping ", length(tips_to_drop), 
            " tips with no clone assignment (e.g. REFERENCE).")
    tree <- drop.tip(tree, tips_to_drop)
  }
  
  if (length(tree$tip.label) < 3) {
    message("  [ERROR] Insufficient tips for cleaning.")
    return(list(tree = tree, 
                clean_clone_df = clone_df[clone_df$cell_id %in% tree$tip.label, 
                                          c("cell_id", "clone_id")]))
  }
  
  message("  [INFO] Working with ", length(tree$tip.label), " tips.")
  
  # --- Build clone lookup ---
  working_clones           <- clone_df[clone_df$cell_id %in% tree$tip.label, ]
  rownames(working_clones) <- working_clones$cell_id
  clone_ids                <- unique(working_clones$clone_id)
  
  message("  [INFO] Clones present: ", paste(sort(clone_ids), collapse = ", "))
  
  # --- Compute node depth (number of descendant tips per node) ---
  Ntip      <- length(tree$tip.label)
  node_ids  <- (Ntip + 1):(Ntip + tree$Nnode)
  depth_all <- node.depth(tree)
  tip_count <- depth_all[node_ids]
  
  # --- For each internal node compute clone fraction among descendants ---
  message("  [INFO] Computing clone fractions per node (this may take a moment)...")
  
  node_clone_frac <- lapply(node_ids, function(node) {
    desc_tip_indices <- Descendants(tree, node, type = "tips")[[1]]
    desc_cells       <- tree$tip.label[desc_tip_indices]
    desc_clones      <- working_clones[desc_cells, "clone_id"]
    tbl              <- table(desc_clones)
    frac             <- tbl / sum(tbl)
    as.data.frame(frac, stringsAsFactors = FALSE)
  })
  
  # --- Find best node per clone: largest clade where clone fraction > min_clone_frac ---
  clean_cells <- c()
  
  for (cln in clone_ids) {
    best_node <- NULL
    best_size <- 0
    
    for (i in seq_along(node_ids)) {
      frac_df <- node_clone_frac[[i]]
      cln_row <- frac_df[frac_df$desc_clones == cln, ]
      
      if (nrow(cln_row) > 0 && cln_row$Freq > min_clone_frac) {
        size <- tip_count[i]
        if (size > best_size) {
          best_size <- size
          best_node <- node_ids[i]
        }
      }
    }
    
    if (!is.null(best_node)) {
      message("  [INFO] Clone ", cln, ": best node ", best_node,
              " with ", best_size, " descendant tips.")
      desc_tip_indices <- Descendants(tree, best_node, type = "tips")[[1]]
      desc_cells       <- tree$tip.label[desc_tip_indices]
      # Keep only tips actually assigned to this clone
      concordant  <- desc_cells[working_clones[desc_cells, "clone_id"] == cln]
      clean_cells <- c(clean_cells, concordant)
    } else {
      message("  [INFO] Clone ", cln, ": no dominant node found, skipping.")
    }
  }
  
  clean_cells <- unique(clean_cells)
  
  if (length(clean_cells) < 3) {
    message("  [SKIP] Too few concordant cells found. Returning unfiltered tree.")
    return(list(tree  = tree,
                clean_clone_df = working_clones[, c("cell_id", "clone_id")]))
  }
  
  tree <- keep.tip(tree, clean_cells)
  
  # --- Build cleaned clone assignment table reflecting cells in final tree ---
  # This is the key output for Sayaka's DEG analysis:
  # only cell_id and clone_id for cells retained after concordance filtering
  clean_clone_df           <- working_clones[clean_cells, c("cell_id", "clone_id")]
  rownames(clean_clone_df) <- NULL
  
  final_count  <- length(tree$tip.label)
  percent_kept <- round((final_count / initial_count) * 100, 2)
  
  message("  [SUMMARY] Concordance Cleaning Results:")
  message("    - Initial cells: ", initial_count)
  message("    - Final cells:   ", final_count)
  message("    - Retained:      ", percent_kept, "%")
  message("    - Clone breakdown:")
  for (cln in sort(unique(clean_clone_df$clone_id))) {
    n <- sum(clean_clone_df$clone_id == cln)
    message("        Clone ", cln, ": ", n, " cells")
  }
  
  return(list(tree = tree, clean_clone_df = clean_clone_df))
}

# ---------------------------
# Setup Paths
# ---------------------------
root_dir    <- rprojroot::find_root(has_dir(".git"))
results_dir <- file.path(root_dir, "analyses", "phylogeny-inference", "results", opt$project)
plots_dir   <- file.path(root_dir, "analyses", "phylogeny-inference", "plots", opt$project)
inputs_dir  <- file.path(root_dir, "analyses", "phylogeny-inference", "inputs")

# Process tree files based on --tree_method
pattern <- "\\.(fasttree|iqtree)\\.tree$"
if (!is.null(opt$tree_method)) {
  pattern <- paste0(".*\\.", opt$tree_method, "\\.tree$")
}
tree_files         <- list.files(results_dir, pattern = pattern, full.names = TRUE)
reference_sequence <- "REFERENCE"

# ============================================================
# MAIN LOOP
# ============================================================
for (tree_path in tree_files) {
  base_name   <- sub("\\.(fasttree|iqtree)\\.tree$", "", basename(tree_path))
  tree_method <- sub(".*\\.(fasttree|iqtree)\\.tree$", "\\1", basename(tree_path))
  
  # Clustering methods available
  clone_methods <- c(
    hierarchical = file.path(results_dir, paste0(base_name, ".clones.hierarchical.tsv.gz")),
    manifold     = file.path(results_dir, paste0(base_name, ".clones.manifold.tsv.gz"))
  )
  
  # Filter by --cluster_method if provided
  if (!is.null(opt$cluster_method)) {
    clone_methods <- clone_methods[names(clone_methods) == opt$cluster_method]
  }
  clone_methods <- clone_methods[file.exists(clone_methods)]
  
  for (clone_type in names(clone_methods)) {
    message("\n>>> Processing ", tree_method, " + ", clone_type, 
            " [Mode: ", opt$annot_type, "]")
    
    cleaned_tree_file   <- file.path(results_dir, 
                                     paste0(base_name, ".", tree_method, ".", 
                                            clone_type, ".clone.treefile"))
    cleaned_clones_file <- file.path(results_dir,
                                     paste0(base_name, ".", tree_method, ".",
                                            clone_type, ".concordant.clones.tsv.gz"))
    
    # --- Load tree baseline ---
    if (opt$annot_type %in% c("Sample", "Lineage") && file.exists(cleaned_tree_file)) {
      message("  [STATUS] Loading existing cleaned tree baseline.")
      tree_current <- read.tree(cleaned_tree_file)
    } else {
      message("  [STATUS] Loading tree baseline.")
      tree_current <- read.tree(tree_path)
      
      if (reference_sequence %in% tree_current$tip.label) {
        tree_current <- root(tree_current, reference_sequence, resolve.root = TRUE)
        if (!opt$refseq) tree_current <- drop.tip(tree_current, reference_sequence)
      } else {
        tree_current <- midpoint.root(tree_current)
      }
      tip_dists    <- distRoot(tree_current)
      tree_current <- drop.tip(tree_current, 
                               names(tip_dists[tip_dists > quantile(tip_dists, 
                                                                    opt$percentile_outlier)]))
    }
    
    # --- Load clone assignments ---
    clone_df <- readr::read_tsv(clone_methods[[clone_type]], show_col_types = FALSE)
    
    # --- Merge clusters if specified ---
    if (!is.null(opt$merge_clusters)) {
      groups <- strsplit(opt$merge_clusters, ",")[[1]]
      for (g in groups) {
        ids    <- strsplit(g, ":")[[1]]
        new_id <- paste(ids, collapse = "")
        clone_df$clone_id[clone_df$clone_id %in% ids] <- new_id
      }
    }
    
    # --- Concordance cleaning ---
    if (opt$annot_type == "Clone" && opt$concordance_clean) {
      result         <- clean_tree_by_concordance(tree_current, clone_df, opt$min_clone_frac)
      tree_current   <- result$tree
      clean_clone_df <- result$clean_clone_df
      
      # Save cleaned tree (Newick format)
      write.tree(tree_current, cleaned_tree_file)
      message("  [SUCCESS] Cleaned tree saved to ", basename(cleaned_tree_file))
      
      # Save concordant clone assignment table (for DEG analysis)
      # Columns: cell_id, clone_id
      # Rows: only cells retained in the concordance-cleaned tree
      readr::write_tsv(clean_clone_df, cleaned_clones_file)
      message("  [SUCCESS] Concordant clone table saved to ", basename(cleaned_clones_file))
      message("            Rows: ", nrow(clean_clone_df), " cells | ",
              "Clones: ", paste(sort(unique(clean_clone_df$clone_id)), collapse = ", "))
    }
    
    # --- Annotation selection ---
    if (opt$annot_type == "Lineage" && !is.null(opt$cell_category)) {
      annot_df  <- readr::read_tsv(file.path(inputs_dir, opt$cell_category), 
                                   show_col_types = FALSE)
      plot_df   <- annot_df %>% filter(cell_id %in% tree_current$tip.label)
      color_var <- names(plot_df)[2]
    } else if (opt$annot_type == "Sample") {
      plot_df   <- clone_df %>% 
        mutate(sample_id = sub("_[^_]+$", "", cell_id)) %>% 
        filter(cell_id %in% tree_current$tip.label)
      color_var <- "sample_id"
    } else {
      plot_df   <- clone_df %>% filter(cell_id %in% tree_current$tip.label)
      color_var <- "clone_id"
    }
    
    # --- Rescale branch lengths ---
    if (opt$rescale) tree_current$edge.length <- log1p(tree_current$edge.length)
    
    # --- Plot ---
    p <- ggtree(tree_current) %<+% plot_df +
      geom_tippoint(aes(color = as.factor(.data[[color_var]])), size = 1.2) +
      theme_tree2() +
      labs(title = paste(opt$project, clone_type, opt$annot_type), 
           color = opt$annot_type)
    
    out_suffix <- tolower(opt$annot_type)
    out_plot   <- file.path(plots_dir, 
                            paste0(base_name, ".", tree_method, ".", 
                                   clone_type, ".", out_suffix, "_tree.png"))
    ggsave(out_plot, p, width = 10, height = 8)
    message("  [PLOT] Saved to ", basename(out_plot))
  }
}

message("\n========== PROCESSING COMPLETE ==========")
