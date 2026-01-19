# Remove compromised cells from either Cell Ranger or DoubletFinder filtered 
# sample feature barcode matrices

# Eric Wafula
# 2025

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(flexmix))


# Set up optparse options
option_list <- list(
  make_option(opt_str = "--project", type = "character", default = NULL,
              help = "Project directory name with 10x Cell Ranger counts sample runs",
              metavar = "character"),
  make_option(opt_str = "--doubletfinder", action = "store_true", default = TRUE, 
              help = "DoubletFinder filtered feature barcode matrices - [default \"%default\"]
                      If set to FALSE, Cell Ranger filtered feature barcode matrices will be used"),
  make_option(opt_str = "--min_cells", type = "integer", default = 3,
              help = "Minimum features detected in at least this many cells [default %default]",
              metavar = "number"),
  make_option(opt_str = "--min_features", type = "integer", default = 200,
              help = "Minimum cells where at least this many features are detected [default %default]",
              metavar = "number"),
  make_option(opt_str = "--percentile_filter", type = "double", default = 0.95,
              help = "Filtering when there is no significant number of compromised cells 
                      (i.e., there is only one distribution detected by the mixture model) [default %default]",
              metavar = "percentile")
)

# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project <- opt$project
doubletfinder <- opt$doubletfinder
min_cells <- opt$min_cells
min_features <- opt$min_features
percentile_filter <- opt$percentile_filter

# Set seed for reproducibility
set.seed(123)

# Establish base directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to module, data, and plots directories
project_dir <- file.path(root_dir, "data", "projects", project)
module_dir <- file.path(root_dir, "analyses", "data-preprocessing")
plots_dir <- file.path(module_dir, "plots", project)
results_dir <- file.path(module_dir, "results", project)


# Create plots folder if it doesn't exist
if (!dir.exists(plots_dir )) {
  dir.create(plots_dir , recursive = TRUE)
}

# Create results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# get samples names
samples <- list.files(project_dir)

for (sample in samples) {
  # Set input and output sub-directories
  input_dir <- file.path(project_dir, sample, "doubletfinder")
  if (!doubletfinder) {
    input_dir <- file.path(project_dir, sample, "outs", "filtered_feature_bc_matrix")
  }
  miqc_dir <- file.path(project_dir, sample, "miqc")
  if (!dir.exists(miqc_dir)) {
    dir.create(miqc_dir)
  }
  
  # Create a Seurat object
  seuratObj <- Seurat::Read10X(data.dir = input_dir) %>% 
    Seurat::CreateSeuratObject(project = sample, min.cell = min_cells, 
                               min.features = min_features)
  
  # Calculate the percent mitochondrial reads in a cell w 
  seuratObj[["percent.mt"]] <- 
    Seurat::PercentageFeatureSet(seuratObj, pattern = "^MT[-\\.]")
  
  # Plot Seurat QC metrics for visual inspection
  Seurat::VlnPlot(seuratObj,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, pt.size=0, layer = "counts")
  ggplot2::ggsave(file.path(plots_dir, paste(sample,
                                    "metrics-distribution.pdf", sep = "-")),
         width = 15)
  p1 <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA",
                               feature2 = "percent.mt")
  p2 <- Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA",
                               feature2 = "nFeature_RNA")
  p1 + p2
  ggplot2::ggsave(file.path(plots_dir, paste(sample, 
                                    "metrics_correlation.pdf", sep = "-")),
         width = 15)
  
  # Automatically determine mitochondrial cutoffs using miQC and flexmix packages
  seuratObj <- SeuratWrappers::RunMiQC(seuratObj,
                    percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA",
                    posterior.cutoff = 0.75, model.slot = "flexmix_model",
                    backup.option = "percentile", backup.percentile = percentile_filter)
  
  if (!is.null(seuratObj@misc$flexmix_model)) {
     # Check MiQC parameters and posterior values
     flexmix::parameters(Misc(seuratObj, "flexmix_model"))
     head(flexmix::posterior(Misc(seuratObj, "flexmix_model")))
    
     # Plot the miQC model results
     SeuratWrappers::PlotMiQC(seuratObj, color.by = "miQC.probability") +
       ggplot2::scale_color_gradient(low = "grey", high = "purple")
     ggplot2::ggsave(file.path(plots_dir, paste(sample, "metrics-miqc-model.pdf", 
                                               sep = "-")))
    # Plot miQC filtering 
    SeuratWrappers::PlotMiQC(seuratObj, color.by = "miQC.keep")
    ggplot2::ggsave(file.path(plots_dir, paste(sample, 
                                               "metrics-miqc-filtering.pdf", 
                                               sep = "-")))
  } else {
    #  Plot percentile filtering 
    Seurat::FeatureScatter(seuratObj,feature1 = "nFeature_RNA", 
                           feature2 = "percent.mt", group.by = "miQC.keep")
    ggplot2::ggsave(file.path(plots_dir, paste(sample, 
                                               "metrics-miqc-filtering.pdf", 
                                               sep = "-")))
  }
  
  # Remove compromised cells from the Seurat object
  seuratObj_filtered <- subset(seuratObj, miQC.keep == "keep")
  seuratObj_filtered
  
  # Save compromised filtered results to file
  DropletUtils::write10xCounts(x = seuratObj_filtered@assays$RNA$counts, 
                               path = miqc_dir, overwrite = TRUE, version='3')
}
