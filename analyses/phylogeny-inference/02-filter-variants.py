import argparse
import os
import sys
import pandas as pd

# Add the parent directory to the Python path
sys.path.append(os.path.dirname(__file__))
import utils.csp_utils as csp_utils

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter pre-computed cellsnp-lite results using inclusion lists.")

    # New required arguments
    parser.add_argument("--project", required=True, help="A valid scPhylogenomics project name string (e.g., TNBC).")
    parser.add_argument("--sample", required=True, help="Sample name string (e.g., TNBC5).")
    parser.add_argument("--cell_type", required=True, help="Cell types being analyzed (e.g., select-cell-types).")

    return parser.parse_args()

def load_filter_lists(data_dir, prefix):
    """
    Loads the inclusion lists for barcodes and SNP sites.

    Expected files:
      - {prefix}cellSNP.barcodes.filtered.tsv.gz (Col: Index)
      - {prefix}cellSNP.sites.filtered.tsv.gz (Cols: CHROM, POS)
    """
    print(f"Loading filter lists from: {data_dir} with prefix: {prefix}")

    # 1. Load Barcode Whitelist
    barcode_file = os.path.join(data_dir, f"{prefix}cellSNP.barcodes.filtered.tsv.gz")
    if not os.path.exists(barcode_file):
        print(f"Error: Barcode filter file not found: {barcode_file}")
        sys.exit(1)

    try:
        bc_df = pd.read_csv(barcode_file, sep="\t")
        if 'Index' not in bc_df.columns:
            print(f"Error: Column 'Index' not found in {barcode_file}")
            sys.exit(1)
        valid_barcodes = set(bc_df['Index'])
        print(f"Loaded {len(valid_barcodes)} valid barcodes.")
    except Exception as e:
        print(f"Error reading barcode file: {e}")
        sys.exit(1)

    # 2. Load Site Whitelist
    sites_file = os.path.join(data_dir, f"{prefix}cellSNP.sites.filtered.tsv.gz")
    if not os.path.exists(sites_file):
        print(f"Error: Sites filter file not found: {sites_file}")
        sys.exit(1)

    try:
        site_df = pd.read_csv(sites_file, sep="\t")
        required_cols = {'CHROM', 'POS'}
        if not required_cols.issubset(site_df.columns):
            print(f"Error: Columns 'CHROM' and 'POS' required in {sites_file}")
            sys.exit(1)

        # Create a set of tuples (CHROM, POS).
        # Note: CHROM is kept as string for consistent comparison, even if input is int.
        valid_sites = set(zip(site_df['CHROM'].astype(str), site_df['POS']))
        print(f"Loaded {len(valid_sites)} valid SNP sites.")
    except Exception as e:
        print(f"Error reading sites file: {e}")
        sys.exit(1)

    return valid_barcodes, valid_sites

def main():
    args = parse_arguments()

    # Construct paths and prefixes
    # Input Dir: ../snv-calling/results/<project>/
    input_base = "../snv-calling/results"
    input_dir = os.path.join(input_base, args.project)

    # Output Dir: ./results/<project>/
    output_dir = os.path.join("./results", args.project)

    # Prefix: <sample>.<cell_type>.
    prefix = f"{args.sample}.{args.cell_type}."

    # --- UPDATED SECTION: Directory Check/Creation ---
    # Ensure the output directory structure exists (e.g., ./results/TNBC/)
    if not os.path.exists(output_dir):
        print(f"Output directory does not exist. Creating: {output_dir}")
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            print(f"Error creating directory {output_dir}: {e}")
            sys.exit(1)
    else:
        print(f"Output directory exists: {output_dir}")
    # -------------------------------------------------

    if not os.path.exists(input_dir):
        print(f"Error: Input directory does not exist: {input_dir}")
        sys.exit(1)

    # 1. Load Filter Lists (Whitelist)
    valid_barcodes, valid_sites = load_filter_lists(input_dir, prefix)

    # 2. Load Existing cellSNP Data
    print("Loading existing cellSNP output...")
    try:
        # csp_utils loads the .gz files with the correct prefix
        adata = csp_utils.csp_load_data(input_dir, prefix=prefix, is_gzip=True)
    except Exception as e:
        print(f"Error loading cellSNP data: {e}")
        sys.exit(1)

    original_vars = adata.n_obs
    original_cells = adata.n_vars # csp_utils loads samples into .var
    print(f"Original data: {original_vars} variants x {original_cells} cells")

    # 3. Apply Filtering

    # A. Filter Variants (Rows - adata.obs)
    # Ensure CHROM column is treated as string for comparison with valid_sites set
    current_chroms = adata.obs['CHROM'].astype(str)
    current_pos = adata.obs['POS']

    variant_mask = [
        (c, p) in valid_sites
        for c, p in zip(current_chroms, current_pos)
    ]

    # B. Filter Cells (Columns - adata.var)
    # csp_utils loads the samples.tsv into adata.var with column name "cell"
    cell_mask = adata.var['cell'].isin(valid_barcodes)

    # Apply both masks simultaneously (Automatic Concordance)
    adata_filtered = adata[variant_mask, cell_mask]

    n_vars_final = adata_filtered.n_obs
    n_cells_final = adata_filtered.n_vars
    print(f"Filtered data: {n_vars_final} variants x {n_cells_final} cells")

    if n_vars_final == 0:
        print("WARNING: All variants were filtered out.")
    if n_cells_final == 0:
        print("WARNING: All cells were filtered out.")

    # 4. Save Filtered Output
    print(f"Saving filtered results to: {output_dir}")
    try:
        csp_utils.csp_save_data(adata_filtered, output_dir, prefix=prefix)
        print("Save complete.")
    except Exception as e:
        print(f"Error saving data: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
