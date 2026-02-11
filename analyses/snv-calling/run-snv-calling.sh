#!/bin/bash
# Eric Wafula
# 2025

# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

printf '\nStart TNBC single sample SNV calling...\n'

# Single sample patient data with selected tumor cell types - TNBC5
printf '\n-- Get sample cell type tumor cell barcodes for SNV calling...\n'
Rscript --vanilla 01-get-cell-type-barcodes.R \
  --project TNBC \
  --cell_types 'luminal epithelial cell of mammary gland,basal cell'

printf '\n-- Call sample SNVs using the cell-snp-lite pipeline...\n'
python 02-snv-calling.py \
  --project TNBC \
  --num_threads 16

printf '\n-- Create sample filtered SNP MSA and binary matrix ...\n'
python 03-generate-snp-msa.py \
  --project TNBC \
  --min_cells_per_snp 0.065 \
  --min_snps_per_cell 0.065
 
  
printf '\nStart SPECTRUM multiple sample SNV calling...\n'  
  
# Muiltiple sample patient data with cancer cell types
printf '\n-- Get sample cell type tumor cell barcodes for SNV calling...\n'
Rscript --vanilla 01-get-cell-type-barcodes.R \
  --project SPECTRUM \
  --cell_types Ovarian-cancer

printf '\n-- Call sample SNVs using the cell-snp-lite pipeline...\n'
python 02-snv-calling.py \
  --project SPECTRUM \
  --num_threads 16

printf '\n-- Create sample filtered SNP MSA and binary matrix ...\n'
python 03-generate-snp-msa.py \
  --project SPECTRUM \
  --min_cells_per_snp 0.03 \
  --min_snps_per_cell 0.03
  
printf '\nStart AML single sample SNV calling...\n'  
  
# Single sample patient data with all cell types - LE1
printf '\n-- Get sample cell type tumor cell barcodes for SNV calling...\n'
Rscript --vanilla 01-get-cell-type-barcodes.R \
  --project AML \
  --cell_types all

printf '\n-- Call sample SNVs using the cell-snp-lite pipeline...\n'
python 02-snv-calling.py \
  --project AML \
  --num_threads 16

printf '\n-- Create sample filtered SNP MSA and binary matrix ...\n'
python 03-generate-snp-msa.py \
  --project AML \
  --min_cells_per_snp 0.015 \
  --min_snps_per_cell 0.05

printf '\nSNV calling Done...\n'


