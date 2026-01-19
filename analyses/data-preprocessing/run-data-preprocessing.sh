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

printf '\nStart data preprocessing...\n'

# TNBC dataset
printf '\n-- Preprocessing TNBC samples...\n'

printf '\n---- Removing background RNA...\n'
Rscript --vanilla 01-remove-ambient-rna.R \
  --project TNBC

printf '\n---- Removing doublets...\n'
Rscript --vanilla  02-remove-doublets.R \
  --project TNBC \
  --soupx TRUE \
  --components 20 \
  --reagent_kit v3.1_Automated

printf '\n---- Removing compromised cells...\n'
Rscript --vanilla 03-remove-compromised.R \
  --project TNBC \
  --doubletfinder TRUE \
  --min_cells 3 \
  --min_features 200 \
  --percentile_filter 0.95
  
# AML dataset
printf '\n-- Preprocessing AML samples...\n'

printf '\n---- Removing background RNA...\n'
Rscript --vanilla 01-remove-ambient-rna.R \
  --project AML

printf '\n---- Removing doublets...\n'
Rscript --vanilla  02-remove-doublets.R \
  --project AML \
  --soupx TRUE \
  --components 20 \
  --reagent_kit v3.1_Automated

printf '\n---- Removing compromised cells...\n'
Rscript --vanilla 03-remove-compromised.R \
  --project AML \
  --doubletfinder TRUE \
  --min_cells 3 \
  --min_features 200 \
  --percentile_filter 0.95

printf '\nData preprocessing Done...\n'
