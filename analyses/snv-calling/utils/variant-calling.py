import argparse
import os
import subprocess
import sys
import pandas as pd
import numpy as np
import shutil
import csp_utils

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run cellsnp-lite and filter output based on PoN and Editing sites.")

    # cellsnp-lite arguments
    parser.add_argument("-s", "--sam", required=True, help="Indexed BAM sample file.")
    parser.add_argument("-b", "--barcode", required=True, help="A plain file listing all effective cell barcodes in the BAM File.")
    parser.add_argument("-O", "--outdir", required=True, help="Output directory.")
    parser.add_argument("-p", "--nproc", type=int, default=4, help="Number of threads to use. Default is 4.")
    parser.add_argument("--minMAF", type=float, default=0.1, help="Minimum MAF. Default 0.1.")
    parser.add_argument("--minCOUNT", type=int, default=100, help="Minimum COUNT. Default 100.")

    # Metadata and Filtering arguments
    parser.add_argument("--id", required=True, help="Sample name string.")
    parser.add_argument("--cell_type", required=True, help="Cell types being analyzed.")
    parser.add_argument("--editing", required=False, help="Path to RNA editing sites file (columns: chr, pos).")
    parser.add_argument("--pon", required=False, help="Path to Panel of Normals (PoN) file.")

    return parser.parse_args()

def normalize_chrom(chrom):
    """
    Normalizes chromosome names to handle mismatches like 'chr1' vs '1'.
    Removes 'chr' prefix and converts to string.

    NOTE: This is used ONLY for comparison/lookup keys. It does not alter 
    the original VCF data content to ensure output format stability.
    """
    return str(chrom).replace("chr", "")

def load_exclusion_list(filepath, file_type):
    """
    Loads PoN or Editing sites into a set of (normalized_chrom, pos) tuples.
    """
    print(f"Loading exclusion list: {filepath}")
    blacklist = set()

    try:
        if file_type == "pon":
            # PoN has a header: #CHROM POS ...
            df = pd.read_csv(filepath, sep="\t")
            # Check common column names
            chrom_col = '#CHROM' if '#CHROM' in df.columns else 'CHROM'
            pos_col = 'POS'
        elif file_type == "editing":
            # Editing file example showed no header: chr10 100001056
            # We assume no header, columns 0 and 1
            df = pd.read_csv(filepath, sep="\t", header=None)
            chrom_col = 0
            pos_col = 1

        for _, row in df.iterrows():
            # CONVERSION: We convert 'chr1' -> '1' here for the lookup set.
            chrom = normalize_chrom(row[chrom_col])
            pos = int(row[pos_col])
            blacklist.add((chrom, pos))

    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        sys.exit(1)

    print(f"Loaded {len(blacklist)} sites from {file_type}.")
    return blacklist

def run_cellsnp(args):
    """
    Constructs and executes the cellsnp-lite command.
    """
    # Create a temporary subdirectory for raw cellsnp output to avoid naming conflicts
    temp_out_dir = os.path.join(args.outdir, "temp_cellsnp_raw")
    if not os.path.exists(temp_out_dir):
        os.makedirs(temp_out_dir, exist_ok=True)

    cmd = [
        "cellsnp-lite",
        "-s", args.sam,
        "-b", args.barcode,
        "-O", temp_out_dir,
        "-p", str(args.nproc),
        "--minMAF", str(args.minMAF),
        "--minCOUNT", str(args.minCOUNT),
        "--gzip" # Usually recommended for VCF output
    ]

    print(f"Running command: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print("cellsnp-lite finished successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running cellsnp-lite:\n{e.stderr}")
        sys.exit(1)

    return temp_out_dir

def filter_and_save(temp_dir, final_dir, args, blacklist):
    """
    Loads data, filters based on blacklist, and saves with new names.
    """
    print("Loading cellSNP output for filtering...")

    # Load data using csp_utils
    # Assuming output is gzipped because we passed --gzip to cellsnp-lite
    try:
        adata = csp_utils.csp_load_data(temp_dir, is_genotype=False, is_gzip=True)
    except Exception:
        # Fallback if gzip failed or wasn't generated
        print("Gzip load failed, trying uncompressed...")
        adata = csp_utils.csp_load_data(temp_dir, is_genotype=False, is_gzip=False)

    original_count = adata.n_obs
    print(f"Original variant count: {original_count}")

    # Create Filter Mask
    keep_mask = []
    vcf_df = adata.obs

    # Determine column names in loaded VCF (csp_utils renames #CHROM to CHROM usually)
    chrom_key = 'CHROM' if 'CHROM' in vcf_df.columns else '#CHROM'

    filtered_counter = 0
    for _, row in vcf_df.iterrows():
        # NORMALIZE FOR CHECK ONLY:
        # We convert the VCF chromosome to standard string (e.g., "1") to match the blacklist.
        # This DOES NOT change the value in the dataframe, so the output format is preserved.
        chrom = normalize_chrom(row[chrom_key])
        pos = int(row['POS'])

        if (chrom, pos) in blacklist:
            keep_mask.append(False)
            filtered_counter += 1
        else:
            keep_mask.append(True)

    print(f"Filtered out {filtered_counter} variants based on PoN and Editing lists.")

    # Apply Filter (AnnData slicing handles X, layers, and obs simultaneously)
    adata_filtered = adata[keep_mask, :]
    print(f"Final variant count: {adata_filtered.n_obs}")

    if adata_filtered.n_obs == 0:
        print("WARNING: All variants were filtered out!")

    # Save the filtered data back to the temp dir first using standard naming
    # We use csp_utils to ensure matrix integrity.
    # Since we never modified adata.obs['CHROM'] in place, it saves as "1", "2", etc.
    csp_utils.csp_save_data(adata_filtered, temp_dir)

    # Now Rename and Move to Final Directory
    prefix = f"{args.id}.{args.cell_type}"

    # Map standard cellsnp names to desired output names
    files_map = {
        "cellSNP.tag.AD.mtx": f"{prefix}.cellSNP.tag.AD.mtx",
        "cellSNP.tag.DP.mtx": f"{prefix}.cellSNP.tag.DP.mtx",
        "cellSNP.tag.OTH.mtx": f"{prefix}.cellSNP.tag.OTH.mtx",
        "cellSNP.base.vcf.gz": f"{prefix}.cellSNP.base.vcf.gz",
        "cellSNP.samples.tsv": f"{prefix}.cellSNP.samples.tsv" # Usually good to keep samples too
    }

    print(f"Renaming and moving files to {final_dir}...")
    if not os.path.exists(final_dir):
        os.makedirs(final_dir, exist_ok=True)

    for old_name, new_name in files_map.items():
        src = os.path.join(temp_dir, old_name)
        dst = os.path.join(final_dir, new_name)

        if os.path.exists(src):
            shutil.move(src, dst)
            print(f"Saved: {new_name}")
        else:
            # Handle case where VCF might not be GZ depending on saved state
            if old_name.endswith(".gz") and os.path.exists(src[:-3]):
                shutil.move(src[:-3], dst[:-3]) # Move unzipped version
                print(f"Saved: {new_name[:-3]}")
            else:
                print(f"Warning: Expected output file {old_name} not found.")

    # Clean up temp dir
    shutil.rmtree(temp_dir)
    print("Processing complete.")

def main():
    args = parse_arguments()

    # 1. Load Exclusion Lists (PoN and Editing)
    # Logic to check if arguments exist before loading
    
    # Handle PoN
    if args.pon:
        pon_blacklist = load_exclusion_list(args.pon, "pon")
    else:
        print("No PoN file provided. Skipping PoN filtering.")
        pon_blacklist = set()

    # Handle RNA Editing
    if args.editing:
        editing_blacklist = load_exclusion_list(args.editing, "editing")
    else:
        print("No Editing sites file provided. Skipping Editing site filtering.")
        editing_blacklist = set()

    # Combine blacklists (works even if sets are empty)
    full_blacklist = pon_blacklist.union(editing_blacklist)

    # 2. Run cellsnp-lite
    temp_dir = run_cellsnp(args)

    # 3. Filter Output and Save with new names
    filter_and_save(temp_dir, args.outdir, args, full_blacklist)

if __name__ == "__main__":
    main()
