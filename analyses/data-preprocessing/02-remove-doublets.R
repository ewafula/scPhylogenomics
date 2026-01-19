# Remove doublets and multiplets cells from either Cell Ranger or SoupX filtered
# sample feature barcode matrices

# Eric Wafula
# 2025

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(DropletUtils))

# Functions
create_seurat_object <- function(soupx_dir, sample, components) {
  # create seurat object
  sample_data <- Seurat::Read10X(data.dir = soupx_dir)
  seuratObj <- Seurat::CreateSeuratObject(counts = sample_data, 
                                              project = sample)
  # pre-process Seurat object (standard)
  seuratObj <- Seurat::NormalizeData(seuratObj)
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, selection.method = "vst",
                                            nfeatures = 2000)  
  seuratObj <- Seurat::ScaleData(seuratObj)
  seuratObj <- Seurat::RunPCA(seuratObj)
  seuratObj <- Seurat::FindNeighbors(seuratObj, dims = 1:components)
  seuratObj <- Seurat::FindClusters(seuratObj)
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = 1:components)
  return(seuratObj)
}

determine_doublet_rate <- function(seuratObj, num_cells, reagent_kit) {
  if (reagent_kit == "v3.1_DualIndex") {
    if (num_cells > 0 && num_cells <= 550)
      { doublet_rate <- 0.004 }
    if (num_cells > 550 && num_cells <= 1500)
      { doublet_rate <- 0.008 }
    if (num_cells > 1500 && num_cells <= 2500)
      { doublet_rate <- 0.016 }
    if (num_cells > 2500 && num_cells <= 3500)
      { doublet_rate <- 0.024 }
    if (num_cells > 3500 && num_cells <= 4500)
      { doublet_rate <- 0.032 }
    if (num_cells > 4500 && num_cells <= 5500)
      { doublet_rate <- 0.040 }
    if (num_cells > 5500 && num_cells <= 6500)
      { doublet_rate <- 0.048 }
    if (num_cells > 6500 && num_cells <= 7500)
      { doublet_rate <- 0.056 }
    if (num_cells > 7500 && num_cells <= 8500)
      { doublet_rate <- 0.064 }
    if (num_cells > 8500 && num_cells <= 9500)
      { doublet_rate <- 0.072 }
    if (num_cells > 9500)
      { doublet_rate <- 0.080 }
  }
  
  if (reagent_kit == "v3.1_Automated"){ 
    if (num_cells > 0 && num_cells <= 550)
      { doublet_rate <- 0.004 }
    if (num_cells > 550 && num_cells <= 1500)
      { doublet_rate <- 0.008 }
    if (num_cells > 1500 && num_cells <= 2500)
      { doublet_rate <- 0.016 }
    if (num_cells > 2500 && num_cells <= 3500)
      { doublet_rate <- 0.023 }
    if (num_cells > 3500 && num_cells <= 4500)
      { doublet_rate <- 0.031 }
    if (num_cells > 4500 && num_cells <= 5500)
      { doublet_rate <- 0.039 }
    if (num_cells > 5500 && num_cells <= 6500)
      { doublet_rate <- 0.046 }
    if (num_cells > 6500 && num_cells <= 7500)
      { doublet_rate <- 0.054 }
    if (num_cells > 7500 && num_cells <= 8500)
      { doublet_rate <- 0.061 }
    if (num_cells > 8500 && num_cells <= 9500)
      { doublet_rate <- 0.069 }
    if (num_cells > 9500)
      { doublet_rate <- 0.076 }    
  }
  
  if (reagent_kit == "v4") { 
    if (num_cells > 0 && num_cells <= 550)
      { doublet_rate <- 0.002 }
    if (num_cells > 550 && num_cells <= 1500)
      { doublet_rate <- 0.004 }
    if (num_cells > 1500 && num_cells <= 2500)
      { doublet_rate <- 0.008 }
    if (num_cells > 2500 && num_cells <= 3500)
      { doublet_rate <- 0.012 }
    if (num_cells > 3500 && num_cells <= 4500)
      { doublet_rate <- 0.016 }
    if (num_cells > 4500 && num_cells <= 5500)
      { doublet_rate <- 0.020 }
    if (num_cells > 5500 && num_cells <= 6500)
      { doublet_rate <- 0.024 }
    if (num_cells > 6500 && num_cells <= 7500)
      { doublet_rate <- 0.028 }
    if (num_cells > 7500 && num_cells <= 8500)
      { doublet_rate <- 0.032 }
    if (num_cells > 8500 && num_cells <= 9500)
      { doublet_rate <- 0.036 }
    if (num_cells > 9500 && num_cells <= 10500)
      { doublet_rate <- 0.040 }
    if (num_cells > 9500 && num_cells <= 10500)
      { doublet_rate <- 0.040 }
    if (num_cells > 10500 && num_cells <= 11500)
      { doublet_rate <- 0.044 }
    if (num_cells > 11500 && num_cells <= 12500)
      { doublet_rate <- 0.048 }
    if (num_cells > 12500 && num_cells <= 13500)
      { doublet_rate <- 0.052 }
    if (num_cells > 13500 && num_cells <= 14500)
      { doublet_rate <- 0.056 }
    if (num_cells > 14500 && num_cells <= 15500)
      { doublet_rate <- 0.060 }
    if (num_cells > 15500 && num_cells <= 16500)
      { doublet_rate <- 0.064 }
    if (num_cells > 16500 && num_cells <= 17500)
      { doublet_rate <- 0.068 }
    if (num_cells > 17500 && num_cells <= 18500)
      { doublet_rate <- 0.072 }
    if (num_cells > 18500 && num_cells <= 19500)
      { doublet_rate <- 0.076 }
    if (num_cells > 19500 )
      { doublet_rate <- 0.080 }
  }
  return(doublet_rate)
} 

get_pK <- function(seuratObj, sample, components, results_dir) {
  # calculate pk
  sweep.res.list <- DoubletFinder::paramSweep(seuratObj, 
                                                 PCs = 1:components, 
                                                 sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  # write pk values to file
  file_name <- file.path(results_dir, paste0(sample,'-pk-values.txt'))
  write.table(bcmvn, file_name, append = FALSE, sep = '\t', dec = '.', 
              quote=FALSE, row.names = FALSE, col.names = TRUE)
  # sweep plot
  pdf(file = file.path(plots_dir, paste0(sample, '-pk-sweep-plot.pdf')))
  pK=as.numeric(as.character(bcmvn$pK))
  BCmetric=bcmvn$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  par(mar=c(5,4,4,1)+.1,cex.main=1.2,font.main=2)
  plot(x = pK, y = BCmetric, pch = 16, type="b",
       col = "blue",lty=1, cex = 1.2)
  abline(v=pK_choose,lwd=2,col='red',lty=2)
  title("The BCmvn distributions")
  text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 2,col = "red")
  dev.off()
  # get the pK values 
  pK <- readr::read_tsv(file_name) %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>% 
    dplyr::pull(pK)
  return(pK)
}

run_doubletfinder <- function(seuratObj, components, doublet_rate, pN, pK) {
  doublet_rate <- as.numeric(doublet_rate)
  pK <- as.numeric(pK)
  annotations <- seuratObj@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate*nrow(seuratObj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seuratObj <- 
    DoubletFinder::doubletFinder(seuratObj, PCs = 1:components, pN = pN, 
                                 pK = pK, nExp = nExp_poi.adj, 
                                 reuse.pANN = NULL, sct = FALSE) 
  # umap plot
  pdf(file = file.path(plots_dir, paste0(sample, '-umap-doublets.pdf')))
  metadata <- names(seuratObj@meta.data)
  metadata <- metadata[grep("^DF\\.classifications_", metadata)]
  print(Seurat::DimPlot(seuratObj, reduction = "umap", group.by = metadata))
  dev.off()
  return(seuratObj)
  seurat_object
}

# Set up optparse options
option_list <- list(
  make_option(opt_str = "--project", type = "character", default = NULL,
              help = "Project directory name with 10x Cell Ranger counts sample runs",
              metavar = "character"),
  make_option(opt_str = "--soupx", action = "store_true", default = TRUE, 
              help = "SoupX filtered feature barcode matrices - [default \"%default\"]
                      If set to FALSE, Cell Ranger filtered feature barcode matrices will be used"),  
  make_option(opt_str = "--components", type = "integer", default = 20,
              help = "Number of components for PCA [default %default]",
              metavar = "number"),
  make_option(opt_str = "--reagent_kit", type = "character", default = "v3.1_Automated",
              help = "10x Chromium Single Cell 3' Reagent Kit used for sample preparation - [default \"%default\"]
                      - For determining sample doublet rate from recovered cells (v3.1_DualIndex, v3.1_Automated or v4)",
              metavar = "character")
)

# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project <- opt$project
sample <- opt$sample
soupx <- opt$soupx
components <- opt$components
reagent_kit <- opt$reagent_kit

# Assign variables
pN <- 0.25 # default recommended for DoubletFinder

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
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Create results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# get samples names
samples <- list.files(project_dir)

for (sample in samples) {
  # Set input and output sub-directories
  input_dir <- file.path(project_dir, sample, "soupx")
  if (!soupx) {
    input_dir <- file.path(project_dir, sample, "outs", "filtered_feature_bc_matrix")
  }
  doubletfinder_dir <- file.path(project_dir, sample, "doubletfinder")
  if (!dir.exists(doubletfinder_dir)) {
    dir.create(doubletfinder_dir)
  }
  
  # create Seurat object and determine doublet rate
  seuratObj <- create_seurat_object(input_dir, sample, components)
  num_cells <- seuratObj %>%  Seurat::Cells() %>% length() 
  doublet_rate <- determine_doublet_rate(seuratObj, num_cells, reagent_kit)
  
  # Run doublet removal analysis 
  pK <- get_pK(seuratObj, sample, components, results_dir)
  seuratObj_doublets <- run_doubletfinder(seuratObj, components, doublet_rate, 
                                          pN, pK)
  
  # Remove doublets from filtered Seurat object
  metadata <- names(seuratObj_doublets@meta.data)
  metadata <- metadata[grep("^DF\\.classifications_", metadata)]
  doublets <- 
    seuratObj_doublets[, seuratObj_doublets@meta.data[, metadata] == "Doublet"]
  doublet_ids <- unique(colnames(doublets))
  seuratObj_doublets <- 
    seuratObj_doublets[,!colnames(seuratObj_doublets) %in% doublet_ids]
  
  # Save doublet filtered results
  DropletUtils::write10xCounts(x = seuratObj_doublets@assays$RNA$counts, 
                               path = doubletfinder_dir, overwrite = TRUE, 
                               version='3')
  
  # Removed pK sweep plot created by doubletfinder. A better pK sweep plot  
  # already created using the pK distribution values output file.
  unlink("Rplots.pdf")
}
