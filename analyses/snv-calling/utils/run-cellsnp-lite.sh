#!/bin/bash
# Eric Wafula
# 2025

set -e
set -o pipefail

###############################################
# Validate arguments
###############################################
# [UPDATED] Usage message now reflects that defaults are NO filtering.
if [ "$#" -lt 5 ] || [ "$#" -gt 9 ]; then
    echo "Usage: $0 <project> <sample> <barcodes file> <cell type> [num threads] [min maf] [min count] [editing file] [pon file]"
    echo "Defaults: num threads=4, min maf=0.1, min count=100"
    echo "Defaults: editing file=None, pon file=None (No filtering performed unless paths provided)"
    exit 1
fi

###############################################
# Global paths
###############################################
UTILS=utils
RESULTS="results"
DATA=../../data

###############################################
# Required parameters
###############################################
PROJECT=$1
SAMPLE=$2
BARCODES_FILE=$3
CELL_TYPE=$4

###############################################
# Optional parameters with defaults
###############################################
NPROC=${5:-4}
MINMAF=${6:-0.1}
MINCOUNT=${7:-100}

# Defaults to empty strings ("") unless user provides file paths.
# This means if the user stops typing after arg 7, these variables are empty.
EDITING_FILE=${8:-""}
PON_FILE=${9:-""}

###############################################
# Environment setup
###############################################
CONDA_ENV=$HOME/miniconda3
SCRATCH=../../scratch/cellsnp/$PROJECT/$SAMPLE
mkdir -p $SCRATCH

# Conda activation
source $CONDA_ENV/bin/activate cellsnp_env

printf "\n Working on $PROJECT $SAMPLE $CELL_TYPE cells...\n"

###############################################
# Step 0: Subset barcodes
###############################################
printf '\n Getting sample cell type barcodes...\n'

zcat "$BARCODES_FILE" | \
    awk -v var="$CELL_TYPE" 'NR==1 || $2==var { print }' \
    > "$SCRATCH/$CELL_TYPE.barcode.tsv"

###############################################
# Step 1: Split BAM by cell type
###############################################
printf "\nSplitting alignment file in cell type specific bams...\n"

bam_file="$DATA/projects/$PROJECT/$SAMPLE/outs/possorted_genome_bam.bam"
sample=$SAMPLE
cell_barcodes="$SCRATCH/$CELL_TYPE.barcode.tsv"

output_dir1="$SCRATCH/Step1_BamCellTypes"
mkdir -p "$output_dir1"

python $UTILS/split-bam.py \
    --bam "$bam_file" \
    --meta "$cell_barcodes" \
    --id "$sample" \
    --n_trim 5 \
    --max_nM 5 \
    --max_NH 1 \
    --outdir "$output_dir1"

###############################################
# Step 2: Variant calling
###############################################
printf "\nDetecting somatic mutations...\n"

cell_type="$CELL_TYPE"
output_dir2="$SCRATCH/Step2_VariantCalling"
mkdir -p "$output_dir2"

# [UPDATED] Build command dynamically based on presence of file paths
CMD=(python $UTILS/variant-calling.py)
CMD+=(--sam "$output_dir1/${sample}.${cell_type}.bam")
CMD+=(--barcode "$output_dir1/${sample}.${cell_type}.barcodes.tsv")
CMD+=(--outdir "$output_dir2")
CMD+=(--nproc "$NPROC")
CMD+=(--minMAF "$MINMAF")
CMD+=(--minCOUNT "$MINCOUNT")
CMD+=(--id "$sample")
CMD+=(--cell_type "$cell_type")

# Only add --editing if the user provided a non-empty string AND the file exists
if [[ -n "$EDITING_FILE" ]]; then
    if [[ -f "$EDITING_FILE" ]]; then
        printf " >> Using Editing sites file: $EDITING_FILE\n"
        CMD+=(--editing "$EDITING_FILE")
    else
        echo "Error: Editing file path provided ('$EDITING_FILE') but file not found."
        exit 1
    fi
else
    printf " >> No Editing file provided. Skipping Editing filtering.\n"
fi

# Only add --pon if the user provided a non-empty string AND the file exists
if [[ -n "$PON_FILE" ]]; then
    if [[ -f "$PON_FILE" ]]; then
        printf " >> Using PoN file: $PON_FILE\n"
        CMD+=(--pon "$PON_FILE")
    else
        echo "Error: PoN file path provided ('$PON_FILE') but file not found."
        exit 1
    fi
else
    printf " >> No PoN file provided. Skipping PoN filtering.\n"
fi

# Execute command
"${CMD[@]}"

###############################################
# Step 3: Copy results to module result dir
###############################################
printf "\n Copying final filtered cellsnp-lite results...\n"

results="results/$PROJECT"
mkdir -p "$results"

cp "$output_dir2/${sample}.${cell_type}.cellSNP.base.vcf.gz" \
   "$results/${sample}.${cell_type}.cellSNP.base.vcf.gz"

gzip -c "$output_dir2/${sample}.${cell_type}.cellSNP.samples.tsv" \
   > "$results/${sample}.${cell_type}.cellSNP.samples.tsv.gz"

gzip -c "$output_dir2/${sample}.${cell_type}.cellSNP.tag.AD.mtx" \
   > "$results/${sample}.${cell_type}.cellSNP.tag.AD.mtx.gz"

gzip -c "$output_dir2/${sample}.${cell_type}.cellSNP.tag.DP.mtx" \
   > "$results/${sample}.${cell_type}.cellSNP.tag.DP.mtx.gz"

gzip -c "$output_dir2/${sample}.${cell_type}.cellSNP.tag.OTH.mtx" \
   > "$results/${sample}.${cell_type}.cellSNP.tag.OTH.mtx.gz"

printf "\nAnalysis Done...\n"

###############################################
# Environment cleanup
###############################################
conda deactivate || true
conda deactivate || true
