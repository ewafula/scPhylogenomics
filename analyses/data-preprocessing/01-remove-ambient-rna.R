# Remove background RNA (contaminants) from unfiltered Cell Ranger count output
# using SoupX

# Eric Wafula
# 2025

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(DropletUtils))

# Functions
run_soupx <- function(cellranger_dir, results_dir)
{
  sc <- SoupX::load10X(cellranger_dir)
  sc <- SoupX::autoEstCont(sc, soupQuantile = 0.75, tfidfMin = 0.5)
  results <- SoupX::adjustCounts(sc)
  DropletUtils::write10xCounts(results_dir, results, version='3', overwrite=TRUE)
}

# Set up optparse options
option_list <- list(
  make_option(opt_str = "--project", type = "character", default = NULL,
              help = "Project directory name with output from 10x Cell Ranger count for all samples",
              metavar = "character")
)

# Parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project <- opt$project

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

# Get samples names
samples <- list.files(project_dir)

# Remove contaminants RNA - background (ambient) RNA
for (sample in samples) {
  # set input and output directories
  cellranger_dir <- file.path(project_dir, sample, "outs")
  soupx_dir <- file.path(project_dir, sample, "soupx")
  
  # create SoupX output directory if it doesn't exist
  if (!dir.exists(soupx_dir)) {
    dir.create(soupx_dir)
  }
  
  # remove sample contaminants with SoupX
  pdf(file = file.path(plots_dir, paste0(sample,"-soupx-contaminant.pdf")), 
      height = 8, width = 11)
  run_soupx(cellranger_dir, soupx_dir)
  dev.off()
}
