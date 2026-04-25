# Create Custome Cell Typing RDS

# Eric Wafula
# 2025

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))

# Set up optparse options
option_list <- list(
  make_option(opt_str = "--expr_matrix", type = "character", default = NULL,
              help = "Cell types gene expression raw counts matrix RDS file 
                      with genes as rows and cells (barcodes) as columns
                     - convert to a sparse matrix if too large",
              metavar = "character"),
  make_option(opt_str = "--cell_types", type = "character", default = NULL, 
              help = "A cell type annotation TSV file with the counts matrix
                      cells (barcodes) as 'cell_id' column and annotated cell 
                      labels as 'cell_type' column",
              metavar = "character"),
  make_option(opt_str = "--normalize_method", type = "character", default = "LogNormalize",
              help = "Normalize the count data - [default \"%default\"]
                      - Available choices - LogNormalize, and SCT",
              metavar = "character"),
  make_option(opt_str = "--rds_file", type = "character", default = "reference-data.rds",
              help = "RDS output file name for annotation object - [default \"%default\"]",
              metavar = "character")
)


# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
expr_matrix <- opt$expr_matrix
cell_types <- opt$cell_types
normalize_method <- opt$normalize_method
rds_file <- opt$rds_file

# Check method options
stopifnot(normalize_method %in% c("LogNormalize", "SCT"))

# Set seed for reproducibility
set.seed(123)

# Load expression matrix
expr_matrix <- readRDS(file.path(expr_matrix))

# Load cell type annotation table
cell_types <- readr::read_tsv(file.path(cell_types)) %>% 
  tibble::column_to_rownames(var = "cell_id")

# Create Seurat object
seurat_obj <- Seurat::CreateSeuratObject(counts = expr_matrix)

# Apply normalization  
if (normalize_method == "SCT") {
  seurat_obj <- Seurat::SCTransform(seurat_obj)
} else {
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
}

# Add cell type annotations to Seurat object 
seurat_obj <-
  Seurat::AddMetaData(seurat_obj, cell_types, col.name = "cell_type")

# Save cell types annotated Seurat object
saveRDS(seurat_obj, file = file.path(rds_file))
