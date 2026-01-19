# Get sample cell type tumor cell barcodes for SNV calling

# Eric Wafula
# 2025

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(glue))

# Set up optparse options
option_list <- list(
  make_option(opt_str = "--project", type = "character", default = NULL,
              help = "A valid scPhylogenomics project name",
              metavar = "character"),
  make_option(opt_str = "--cell_types", type = "character", default = "aneuploid",
              help = "Cell types to used SNV calling - [default \"%default\"]
                      - Either comma-separated list of cell type names or a single name,
                      - cell type names in quotations if includes white space(s)
                      - if set to `all`, all cells will be used
                      - if set to `aneuploid`, cancer cells determined using copy 
                      - number changes will be used - dafault",
              metavar = "character")
)

# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project <- opt$project
cell_types <- opt$cell_types

 # Set seed for reproducibility
set.seed(123)

# Establish base directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set data directories
module_dir <- file.path(root_dir, "analyses", "snv-calling")
ploidy_dir <- file.path(root_dir, "analyses", "ploidy-inference", "results", project)
inputs_dir <- file.path(module_dir, "inputs",  project)

# Create inputs folder if it doesn't exist
if (!dir.exists(inputs_dir)) {
  dir.create(inputs_dir, recursive = TRUE)
}

# Load ploidy classification table and subset for cancer cell types of interest 
# i.e., either all or some cell types that have aneuploid cells.
ploidy_results <- file.path(ploidy_dir, 
                            glue::glue("{project}-cell-ploidy-prediction.tsv.gz"))

# Get cell types of interest
if (cell_types == "aneuploid"){
  ploidy_prediction <- readr::read_tsv(ploidy_results) %>% 
    dplyr::filter(cell_ploidy %in% "aneuploid") %>% 
    dplyr::mutate(cell_type = "aneuploid")
} else if (cell_types == "all"){
  ploidy_prediction <- readr::read_tsv(ploidy_results) %>% 
    dplyr::mutate(cell_type = gsub("\\s+", "-", cell_type)) %>% 
    dplyr::mutate(cell_type = "all-cell-types")
} else {
  cells <- unlist(strsplit(cell_types, ","))
  ploidy_prediction <- readr::read_tsv(ploidy_results)
  if (length(cells) > 1) {
    ploidy_prediction <- ploidy_prediction %>%
      dplyr::filter(cell_type %in% cells) %>% 
      dplyr::mutate(cell_type = gsub("\\s+", "-", cell_type)) %>% 
      dplyr::mutate(cell_type = "select-cell-types")
  } else {
    ploidy_prediction <- ploidy_prediction %>% 
      dplyr::filter(cell_type %in% cells) %>% 
      dplyr::mutate(cell_type = gsub("\\s+", "-", cell_type))
  }
}

# Subset by sample to use in SNV calling
for (sample in unique(ploidy_prediction$sample_id)) {
  ploidy_prediction %>% dplyr::filter(sample_id == sample) %>% 
    dplyr::select(cell_id, cell_type) %>% 
    dplyr::rename(Index = cell_id, Cell_type = cell_type) %>% 
    dplyr::mutate(Index = sub(".*_", "", Index)) %>% 
    readr::write_tsv(file.path(inputs_dir, 
                               glue::glue("{sample}-cancer-cells-barcodes.tsv.gz")))
}
