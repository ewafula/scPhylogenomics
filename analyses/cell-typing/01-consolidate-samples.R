# Project Sample Consolidation
 
# Eric Wafula
# 2025

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(scCustomize))

# Functions
create_seurat_object <- function(sample_data_dir, sample_name) {
  # create seurat object
  sample_data <- Seurat::Read10X(data.dir = sample_data_dir)
  seurat_object <- Seurat::CreateSeuratObject(counts = sample_data,
                                              project = sample_name)
  return(seurat_object)
}

# Set up optparse options
option_list <- list(
  make_option(opt_str = "--project", type = "character", default = NULL,
              help = "Project directory name with 10x Cell Ranger counts sample runs",
              metavar = "character"),
  make_option(opt_str = "--integrate", action = "store_true", default = TRUE, 
              help = "Integate project samples - [default \"%default\"]
                      If set to FALSE, project samples will not be integrated"),
  make_option(opt_str = "--normalize_method", type = "character", default = "LogNormalize",
              help = "Normalize the count data - [default \"%default\"]
                      - Available choices - LogNormalize, and SCT",
              metavar = "character"),
  make_option(opt_str = "--integrate_method", type = "character", default = "Harmony",
              help = "Method for integrating project samples - [default \"%default\"]
                      - Available choices - Harmony, CCA, FastMNN, and RPCA",
              metavar = "character"),
  make_option(opt_str = "--components", type = "integer", default = 30,
              help = "Number of components for PCA [default %default]",
              metavar = "number"),
  make_option(opt_str = "--resolution", type = "double", default = 0.8,
              help = "Seurat clustering resolution parameter - [default %default]",
              metavar = "number")
)

# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project <- opt$project
integrate <- opt$integrate
normalize_method <- opt$normalize_method
integrate_method <- opt$integrate_method
components <- opt$components
resolution <- opt$resolution

# Check method options
stopifnot(normalize_method %in% c("LogNormalize", "SCT"))
stopifnot(integrate_method %in% c("Harmony", "CCA", "FastMNN", "RPCA"))
integrate_method <- paste0(integrate_method, "Integration") 

# Set seed for reproducibility
set.seed(123)

# Establish base directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to module, data, and plots directories
project_dir <- file.path(root_dir, "data", "projects", project)
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

# Get samples names
samples <- list.files(project_dir)

# Create a list of sample Seurat objects to merge
sample_obj_list <- list()
sample_names <- vector()
for (sample in samples) {
  # create sample Seurat objects
  sample_obj <- create_seurat_object(file.path(project_dir, sample, 
                                               "miqc"), sample)
  sample_obj_list <- append(sample_obj_list, sample_obj)
  sample_names <- append(sample_names, sample)
}

# Merge samples 
if (length(samples) > 1) {
  # merge
  sample_obj <- merge(x = sample_obj_list[[1]], y = sample_obj_list[-1], 
                    add.cell.ids = sample_names, project = project)
  
  # remove sample object list to conserve memory
  rm(sample_obj_list, sample_names)
} else {
  sample_obj <- SeuratObject::RenameCells(object = sample_obj, 
                                          add.cell.id = sample_names)
}

# Add 'cell_id' and 'sample_id' to metadata
sample_obj$cell_id <- rownames(sample_obj@meta.data)
sample_obj$sample_id <- sample_obj$orig.ident

# Apply normalization and pca dimension reduction 
if (normalize_method == "SCT") {
  sample_obj <- Seurat::SCTransform(sample_obj)
  sample_obj <- Seurat::RunPCA(sample_obj, npcs = 50, verbose = F)
} else {
  sample_obj <- Seurat::NormalizeData(sample_obj)
  sample_obj <- Seurat::FindVariableFeatures(sample_obj)
  sample_obj <- Seurat::ScaleData(sample_obj)
}

# Run a PCA dimensionality reduction and create elbow plot
sample_obj <- Seurat::RunPCA(sample_obj, npcs = 50, verbose = F)
pdf(file=file.path(plots_dir, glue::glue("{project}-elbowplot.pdf")))
print(Seurat::ElbowPlot(sample_obj, ndims = 50))
dev.off()

# integrate samples 
if (integrate) {
  if (normalize_method == "SCT") {
    sample_obj <- Seurat::IntegrateLayers(object = sample_obj, 
                                          method = integrate_method, 
                                          orig.reduction = "pca", 
                                          new.reduction = "integrated",
                                          normalization.method = "SCT", 
                                          verbose = FALSE)
    sample_obj <- Seurat::FindNeighbors(sample_obj, dims = 1:components, 
                                        reduction = "integrated")
    sample_obj <- Seurat::FindClusters(sample_obj, resolution = resolution)
    sample_obj <- Seurat::RunUMAP(sample_obj, dims = 1:components,
                                  reduction = "integrated", 
                                  reduction.name = "umap")
  } else {
    sample_obj <- Seurat::IntegrateLayers(object = sample_obj,
                                          method = integrate_method,
                                          orig.reduction = "pca",
                                          new.reduction = "integrated",
                                          verbose = FALSE)       
    sample_obj <- Seurat::FindNeighbors(sample_obj, dims = 1:components, 
                                        reduction = "pca")
    sample_obj <- Seurat::FindClusters(sample_obj, resolution = resolution)
    sample_obj <- Seurat::RunUMAP(sample_obj, dims = 1:components,
                                  reduction = "pca", reduction.name = "umap")
  }				
} else {
  sample_obj <- Seurat::FindNeighbors(sample_obj, dims = 1:components,
                                      reduction = "pca")
  sample_obj <- Seurat::FindClusters(sample_obj, resolution = resolution)
  sample_obj <- Seurat::RunUMAP(sample_obj, dims = 1:components,
                                reduction = "pca", reduction.name = "umap")
}

# Plot UMAP with cell clusters
pdf(file=file.path(plots_dir, glue::glue("{project}-clusters-umap.pdf")))
scCustomize::DimPlot_scCustom(sample_obj, reduction = "umap",  label = TRUE, 
                              repel = TRUE, group.by = "seurat_clusters")
dev.off()

# Save Seurat object
saveRDS(sample_obj, file = file.path(objects_dir, 
                                     glue::glue("{project}-clusters.rds")))

