#!/bin/bash
# Updated: 2026 - Automatic Sample Detection, Merging, and Results Aggregation
# Author: Eric Wafula

set -e
set -o pipefail

###############################################
# Validate arguments
###############################################
if [ "$#" -lt 2 ] || [ "$#" -gt 7 ]; then
    echo "Usage: $0 <project> <cell type> [num threads] [min maf] [min count] [editing file] [pon file]"
    echo "Note: <sample> and <barcodes file> are automatically detected from the 'inputs' directory."
    exit 1
fi

###############################################
# Global paths
###############################################
UTILS=utils
RESULTS="results"
DATA=../../data
INPUTS="inputs"

###############################################
# Required parameters
###############################################
PROJECT=$1
CELL_TYPE=$2

###############################################
# Optional parameters with defaults
###############################################
NPROC=${3:-4}
MINMAF=${4:-0.1}
MINCOUNT=${5:-100}
EDITING_FILE=${6:-""}
PON_FILE=${7:-""}

###############################################
# Environment setup
###############################################
if [ -d "/opt/conda" ]; then
    CONDA_BASE="/opt/conda"
elif [ -d "$HOME/miniconda3" ]; then
    CONDA_BASE="$HOME/miniconda3"
else
    echo "Error: Could not find miniconda3 in /opt/conda or $HOME/miniconda3"
    exit 1
fi

# Conda activation
source "$CONDA_BASE/bin/activate" cellsnp_env

# Automatically determine samples based on the flat file structure
# Extracts everything before "-cancer-cells-barcodes.tsv.gz" as the SAMPLE name
SAMPLES=($(ls $INPUTS/$PROJECT/*-cancer-cells-barcodes.tsv.gz | xargs -n 1 basename | sed 's/-cancer-cells-barcodes.tsv.gz//'))
NUM_SAMPLES=${#SAMPLES[@]}

if [ "$NUM_SAMPLES" -eq 0 ]; then
    echo "Error: No barcode files found in $INPUTS/$PROJECT matching pattern *-cancer-cells-barcodes.tsv.gz"
    exit 1
fi

printf "\nDetected $NUM_SAMPLES samples in project: $PROJECT\n"

# Arrays to keep track of outputs for merging
BAMS_TO_MERGE=()
BARCODES_TO_MERGE=()

###############################################
# Step 1: Iterate through samples and split BAMs
###############################################
for SAMPLE in "${SAMPLES[@]}"; do
    printf "\n--- Processing Sample: $SAMPLE ---\n"
    
    # NEW: Updated path for flat structure
    BARCODES_FILE="$INPUTS/$PROJECT/$SAMPLE-cancer-cells-barcodes.tsv.gz"
    
    if [ ! -f "$BARCODES_FILE" ]; then
        echo "Warning: Barcode file $BARCODES_FILE not found. Skipping $SAMPLE."
        continue
    fi

    SCRATCH_SAMPLE="../../scratch/cellsnp/$PROJECT/$SAMPLE"
    mkdir -p "$SCRATCH_SAMPLE"

    # Subset barcodes for current sample/cell type
    printf "Filtering barcodes for $CELL_TYPE...\n"
    zcat "$BARCODES_FILE" | \
        awk -v var="$CELL_TYPE" 'NR==1 || $2==var { print }' \
        > "$SCRATCH_SAMPLE/$CELL_TYPE.barcode.tsv"

    # BAM files are still expected in the standard data/project/sample path
    bam_file="$DATA/projects/$PROJECT/$SAMPLE/outs/possorted_genome_bam.bam"
    output_dir1="$SCRATCH_SAMPLE/Step1_BamCellTypes"
    mkdir -p "$output_dir1"

    if [ ! -f "$bam_file" ]; then
        echo "Error: BAM file not found at $bam_file"
        exit 1
    fi

    printf "Splitting BAM...\n"
    python $UTILS/split-bam.py \
        --bam "$bam_file" \
        --meta "$SCRATCH_SAMPLE/$CELL_TYPE.barcode.tsv" \
        --id "$SAMPLE" \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 \
        --outdir "$output_dir1"

    # Save paths for Step 2
    BAMS_TO_MERGE+=("$output_dir1/${SAMPLE}.${CELL_TYPE}.bam")
    BARCODES_TO_MERGE+=("$output_dir1/${SAMPLE}.${CELL_TYPE}.barcodes.tsv")
done

###############################################
# Step 2: Merge logic
###############################################
FINAL_SCRATCH="../../scratch/cellsnp/$PROJECT"
mkdir -p "$FINAL_SCRATCH"

if [ "$NUM_SAMPLES" -gt 1 ]; then
    printf "\nMultiple samples detected. Merging into 'MERGED' BAM...\n"
    FINAL_ID="MERGED"
    FINAL_BAM="$FINAL_SCRATCH/$FINAL_ID.$CELL_TYPE.bam"
    FINAL_BARCODES="$FINAL_SCRATCH/$FINAL_ID.$CELL_TYPE.barcodes.tsv"
    
    # Use helper script to merge, sort, and index
    python $UTILS/merge-bams.py --inputs "${BAMS_TO_MERGE[@]}" --output "$FINAL_BAM"
    
    # Combine unique barcodes across all samples
    cat "${BARCODES_TO_MERGE[@]}" | sort | uniq > "$FINAL_BARCODES"
else
    printf "\nSingle sample detected. Proceeding without merge...\n"
    FINAL_ID="${SAMPLES[0]}"
    FINAL_BAM="${BAMS_TO_MERGE[0]}"
    FINAL_BARCODES="${BARCODES_TO_MERGE[0]}"
fi

###############################################
# Step 3: Variant calling
###############################################
printf "\nDetecting somatic mutations for $FINAL_ID ($CELL_TYPE)...\n"

output_dir2="$FINAL_SCRATCH/Step2_VariantCalling"
mkdir -p "$output_dir2"

CMD=(python $UTILS/variant-calling.py)
CMD+=(--sam "$FINAL_BAM")
CMD+=(--barcode "$FINAL_BARCODES")
CMD+=(--outdir "$output_dir2")
CMD+=(--nproc "$NPROC")
CMD+=(--minMAF "$MINMAF")
CMD+=(--minCOUNT "$MINCOUNT")
CMD+=(--id "$FINAL_ID")
CMD+=(--cell_type "$CELL_TYPE")

[[ -f "$EDITING_FILE" ]] && CMD+=(--editing "$EDITING_FILE")
[[ -f "$PON_FILE" ]] && CMD+=(--pon "$PON_FILE")

# Execute variant calling
"${CMD[@]}"

###############################################
# Step 4: Refined Copy results (with Compression)
###############################################
printf "\nCopying and compressing final results to $RESULTS/$PROJECT...\n"

DEST_DIR="$RESULTS/$PROJECT"
mkdir -p "$DEST_DIR"

SRC_PREFIX="$output_dir2/${FINAL_ID}.${CELL_TYPE}.cellSNP"

# 1. VCF
if [ -f "${SRC_PREFIX}.base.vcf.gz" ]; then
    cp "${SRC_PREFIX}.base.vcf.gz" "$DEST_DIR/"
elif [ -f "${SRC_PREFIX}.base.vcf" ]; then
    gzip -c "${SRC_PREFIX}.base.vcf" > "$DEST_DIR/${FINAL_ID}.${CELL_TYPE}.cellSNP.base.vcf.gz"
fi

# 2. Samples TSV
if [ -f "${SRC_PREFIX}.samples.tsv" ]; then
    gzip -c "${SRC_PREFIX}.samples.tsv" > "$DEST_DIR/${FINAL_ID}.${CELL_TYPE}.cellSNP.samples.tsv.gz"
fi

# 3. Matrix Files (AD, DP, OTH)
for tag in AD DP OTH; do
    if [ -f "${SRC_PREFIX}.tag.${tag}.mtx" ]; then
        gzip -c "${SRC_PREFIX}.tag.${tag}.mtx" > "$DEST_DIR/${FINAL_ID}.${CELL_TYPE}.cellSNP.tag.${tag}.mtx.gz"
    fi
done

printf "\nAnalysis Done. Results available in $DEST_DIR\n"

###############################################
# Environment cleanup
###############################################
conda deactivate || true
