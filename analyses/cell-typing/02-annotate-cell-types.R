# Reference Cell Typing
 
# Eric Wafula
# 2025

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(scCustomize))

# Functions
annotate_cells <- function(sample_obj, matrix_normized, reference, cell_types, 
                           project, annot_method) {
  # reference cell typing using singleR
  cell_predictions <- SingleR::SingleR(test = matrix_normized, ref = reference,
                                       labels = cell_types, de.method = 'wilcox')
  
  # delta heatmap for cell type assignment and dynamically adjust plot width
  max_rowname_length <- max(nchar(rownames(cell_predictions)))
  extra_width <- max_rowname_length * 5 # Scale the width dynamically (tweak factor as needed)
  extra_width <- extra_width / 25 # Convert to inches (R uses inches for file dimensions)
  
  pdf(file =
        file.path(plots_dir, 
                  glue::glue("{project}-{annot_method}-delta-heatmap.pdf")), 
      width = 8 + extra_width, height = 8)
  SingleR::plotScoreHeatmap(cell_predictions, show.pruned = TRUE)
  dev.off()
  
  # add cell type annotations to Seurat object
  cell_predictions <- cell_predictions %>% as.data.frame() %>% 
    dplyr::rename(cell_type = `pruned.labels`)
  sample_obj <-
    Seurat::AddMetaData(sample_obj, cell_predictions, col.name = "cell_type")

  return(sample_obj)
}

# Set up optparse options
option_list <- list(
  make_option(opt_str = "--project", type = "character", default = NULL,
              help = "Project for samples in the Seurat  object",
              metavar = "character"),
  make_option(opt_str = "--assay", type = "character", default = "RNA",
              help = "Project Seurat object assay with normalized counts - [default \"%default\"]
                      - Available choices - RNA and SCT",
              metavar = "character"),
  make_option(opt_str = "--annot_method", type = "character", default = "celldex",
              help = "Automatic cell type annotation method - [default \"%default\"]
                      - Available choices - celldex, custom, and mapping",
              metavar = "character"),
  make_option(opt_str = "--ref_data", type = "character", default = "hpca",
              help = "Reference dataset - [default \"%default\"]
                      - Available celldex reference choices - hpca, monaco_immune, 
                        immgen, dice, novershtern_hematopoietic, blueprint_encode
                      - If 'custom' is set for the '--annot_method' option,
                        provide custom reference Seurat object with normalized 
                        counts similar the project Seurat object (query object)
                        and cell labels column named 'cell_type' 
                      - If 'mapping' is set for '--annot_method' option provide a
                        TSV mapping file of annotations projects samples with the 
                        matching cells (barcodes) as 'cell_id' column and annotated 
                        cell labels as 'cell_type' column",
              metavar = "character")
)

# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project <- opt$project
assay <- opt$assay
annot_method <- opt$annot_method
ref_data <- opt$ref_data

# Check method options
stopifnot(assay %in% c("RNA", "SCT"))
stopifnot(annot_method %in% c("celldex", "custom", "mapping"))

# Set seed for reproducibility
set.seed(123)

# Establish base directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to module, data, and plots directories
module_dir <- file.path(root_dir, "analyses", "cell-typing")
objects_dir <- file.path(module_dir, "objects", project)
plots_dir <- file.path(module_dir, "plots", project)
results_dir <- file.path(module_dir, "results", project)

# Create objects folder if it doesn't exist
if (!dir.exists(objects_dir)) {
  dir.create(objects_dir , recursive = TRUE)
}

# Create plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Create results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Load samples Seurat object and join layers
sample_obj <- readRDS(file.path(objects_dir, glue::glue("{project}-clusters.rds")))
num_of_samples <-  sample_obj$sample_id %>% unique() %>% length()
if (num_of_samples > 1) {
  sample_obj<- SeuratObject::JoinLayers(sample_obj)
}

# Get normalized counts matrix
matrix_normized <- SeuratObject::LayerData(sample_obj, assay = assay, 
                                           layer = 'data')

if (annot_method == "celldex") {
  # check celldex reference dataset options
  stopifnot(ref_data %in% c("hpca", "monaco_immune", "immgen", "dice", 
                             "novershtern_hematopoietic", "blueprint_encode"))
  
  # get reference dataset from celldex
  version <- celldex::fetchLatestVersion(ref_data)
  reference <- celldex::fetchReference(ref_data, version)
  cell_types <- reference$label.main
  sample_obj <- annotate_cells(sample_obj, matrix_normized, reference, 
                               cell_types, project, annot_method)
  
} else if (annot_method == "custom") {
  # load SCE RDS reference object and gene cell types 
  reference <- readRDS(file.path(ref_data))
  expr_matrix <- SeuratObject::GetAssayData(reference, assay = assay, 
                                            layer = "data")
  cell_types <- reference$cell_type
  sample_obj <- annotate_cells(sample_obj, matrix_normized, expr_matrix, 
                               cell_types, project, annot_method)
} else {
  # load 'cell_id' to 'cell_type' file
  cell_types <- readr::read_tsv(file.path(ref_data)) %>% 
    tibble::column_to_rownames(var = "cell_id")
  
  # add cell type annotations to Seurat object
  sample_obj <-
    Seurat::AddMetaData(sample_obj, cell_types, col.name = "cell_type")
}

# Remove cells without cell type annotations
sample_obj <- subset(x = sample_obj, subset = cell_type != is.na(cell_type))

# Plot UMAP with cell types
pdf(file=file.path(plots_dir, glue::glue("{project}-{annot_method}-annotations-umap.pdf")), 
    width = 10)
scCustomize::DimPlot_scCustom(sample_obj, reduction = "umap", 
                              group.by = "cell_type")
dev.off()

# Save cell annotations
sample_obj@meta.data %>% dplyr::select(cell_id, sample_id, seurat_clusters, 
                                       cell_type) %>% 
  readr::write_tsv(file.path(results_dir, 
                             glue::glue("{project}-{annot_method}-annotations.tsv.gz")))

# Save cell annotated Seurat object
saveRDS(sample_obj, file = file.path(objects_dir, 
                                     glue::glue("{project}-{annot_method}-annotations.rds")))
