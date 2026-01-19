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

printf '\nStart ploidy prediction and classification...\n'

# TNBC dataset individual sample
printf '\n-- Cell ploidy prediction and classification of TNBC individual sample...\n'
  
Rscript --vanilla 01-classify-cell-ploidy.R \
  --project TNBC \
  --annot_object TNBC-custom-annotations.rds \
  --nthreads 30

# AML dataset individual sample
printf '\n-- Cell ploidy prediction and classification of AML individual sample...\n'

Rscript --vanilla 01-classify-cell-ploidy.R \
  --project AML \
  --annot_object AML-mapping-annotations.rds \
  --nthreads 30

printf '\nCell ploidy prediction and classification Done...\n'
