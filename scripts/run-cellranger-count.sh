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


# Check if required arguments are provided
if [ $# -ne 4 ]; then
    printf "Usage: $0 <PROJECT_NAME> <SAMPLE SHEET> <REFERENCE DATA DIR> <THREADS>\n"
    exit 1
fi

# Set up directories to store project sample feature, barcode and matrix data 
DATA_DIR="../data"
PROJECTS_DIR="$DATA_DIR/projects"
PROJECT_NAME="$1"
SAMPLE_SHEET="$2"
REFERENCE_DIR="$3"
THREADS="$4"

# Crete project data sub-directory
PROJECT_DIR="$PROJECTS_DIR/$PROJECT_NAME"
mkdir -p $PROJECT_DIR

printf "\nStart creating sample feature, barcode and matrix data for project scRNA-Seq fastq ...\n\n"

# Create sample matrix data and copy move out to project directory
while IFS=$'\t' read -r sample fastq_dir _; do
  
  printf "Running Cell Ranger count on $sample sample fastq...\n\n"
  
  # run Cell Ranger count
  cellranger count \
      --id=$sample \
      --fastqs=$fastq_dir \
      --sample=$sample \
      --transcriptome=$REFERENCE_DIR \
      --create-bam=true \
      --localcores=$THREADS
  
  # move sample feature, barcode and matrix data directory to the project directory
  SAMPLE_DIR="$sample"
  mv $SAMPLE_DIR $PROJECT_DIR
  
done < "$SAMPLE_SHEET"

printf "\nDone running cellranger count...\n\n"
