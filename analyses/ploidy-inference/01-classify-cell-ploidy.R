# Classify Sample into Tumor and Normal Cells 

# Eric Wafula
# 2025

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(copykat))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(scCustomize))

# Functions
run_copykat <- function(count_matrix, project, plot_genes, nthreads, plots_dir) {
  # predict cell ploidy
  cell_ploidy <- copykat::copykat(rawmat = count_matrix, id.type = "Symbol",
                                  ngene.chr = 5, win.size = 25, KS.cut = 0.1, 
                                  sam.name = project, distance = "euclidean", 
                                  norm.cell.names = "", output.seg = "FALSE", 
                                  plot.genes = plot_genes, genome = "hg20",
                                  n.cores = nthreads)
  
  # save cell-by-gene copy number heatmap 
  src <- glue::glue("{project}_copykat_heatmap.jpeg")
  dest <- file.path(plots_dir, basename(src))
  
  # remove CopyKAT output files not required
  file.copy(src, dest, overwrite = TRUE)
  files <- list.files(pattern = project)
  file.remove(files)
  
  return(cell_ploidy)
}

# Set up optparse options
option_list <- list(
  make_option(opt_str = "--project", type = "character", default = NULL,
              help = "A valid scPhylogenomics project name",
              metavar = "character"),
  make_option(opt_str = "--annot_object", type = "character", default = NULL,
              help = "Project cell type annotated Seurat object RDS file name",
              metavar = "character"),
  make_option(opt_str = "--plot_genes", action = "store_true", default = TRUE,
              help = "Plot copy number heatmap with gene by cell matrix - [default \"%default\"]
              If set to FALSE, copy number heatmap not be generated"),
  make_option(opt_str = "--nthreads", type = "integer", default = 4,
              help = "Number of threads (CPUs) [default = %default]",
              metavar = "character")
)

# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project <- opt$project
annot_object <- opt$annot_object
plot_genes <-opt$plot_genes
nthreads <- opt$nthreads

# Set seed for reproducibility
set.seed(123)

# Establish base directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set path to module, data, and plots directories
module_dir <- file.path(root_dir, "analyses", "ploidy-inference")
annot_dir <- file.path(root_dir, "analyses", "cell-typing", "objects")
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

# Load cell annotation object and extract count matrix
sample_obj <- readRDS(file.path(annot_dir, project, annot_object))
count_matrix <- SeuratObject::GetAssayData(sample_obj, assay = "RNA", 
                                           layer = "counts")
# Predict cell ploidy
cell_ploidy <- run_copykat(count_matrix, project, plot_genes, nthreads, plots_dir)

# Classify cells ploidy and to Seurat object metadata
cell_class <- cell_ploidy$prediction %>%  tibble::as_tibble() %>% 
  tibble::column_to_rownames(var = "cell.names")  %>% 
  dplyr::rename(cell_ploidy = `copykat.pred`) %>% 
  dplyr::mutate(cell_class = case_when(cell_ploidy == "diploid" ~ "normal",
                                       cell_ploidy == "aneuploid" ~ "tumor",
                                       cell_ploidy == "not.defined" ~ "unknown"))
sample_obj <- Seurat::AddMetaData(sample_obj, cell_class)

# Remove cells with cell type annotations
sample_obj <- subset(x = sample_obj, subset = cell_class != is.na(cell_class))


# Plot UMAP with cell types
pdf(file=file.path(plots_dir, glue::glue("{project}-cell-ploidy-prediction-umap.pdf")), 
    width = 17)
scCustomize::DimPlot_scCustom(sample_obj, reduction = "umap",
                              group.by = c("cell_type", "cell_class"))
dev.off()

# Save cell annotations
sample_obj@meta.data %>% dplyr::select(cell_id, sample_id, seurat_clusters, 
                                       cell_type, cell_type, cell_ploidy, cell_class) %>% 
  readr::write_tsv(file.path(results_dir, 
                             glue::glue("{project}-cell-ploidy-prediction.tsv.gz")))

# Save cell annotated Seurat object
saveRDS(sample_obj, file = file.path(objects_dir, 
                                     glue::glue("{project}-cell-ploidy-prediction.rds")))
