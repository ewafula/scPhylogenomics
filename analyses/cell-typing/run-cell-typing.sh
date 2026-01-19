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

printf '\nStart cell typing...\n'
# TNBC dataset individual sample
printf '\n-- Cell typing of TNBC individual sample...\n'
Rscript --vanilla 01-consolidate-samples.R \
  --project TNBC \
  --integrate FALSE \
  --normalize_method LogNormalize \
  --components 20 \
  --resolution 0.4
  
printf '\n---- Annotating cell types using custom HBCA dataset...\n'  
Rscript --vanilla  02-annotate-cell-types.R \
  --project TNBC \
  --assay RNA \
  --annot_method custom \
  --ref_data inputs/HBCA-snRNA-seq-all-cells.rds
  
# AML dataset individual sample
printf '\n-- Cell typing of AML individual sample...\n' 
Rscript --vanilla 01-consolidate-samples.R \
  --project AML \
  --integrate FALSE \
  --normalize_method LogNormalize \
  --components 20 \
  --resolution 1.4
  
printf '\n---- Annotating cell types using the Leukemia Cell Atlas...\n'
pushd ../../scripts/ > /dev/nul
Rscript -e "rmarkdown::render('create-leukemia-reference-mapping.Rmd', clean = TRUE)"
popd > /dev/null
  
printf '\n---- Assigning cell types using mappings from Leukemia Cell Atlas...\n'
Rscript --vanilla  02-annotate-cell-types.R \
  --project AML \
  --assay RNA \
  --annot_method mapping \
  --ref_data inputs/AML-LE1-cell-annotation.tsv.gz

printf '\nCell typing Done...\n'


