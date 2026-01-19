import pandas as pd
import numpy as np
import gzip
import os
import argparse
import glob
import sys
import math

# --- Logic for File Discovery and Processing ---

def parse_arguments():
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(description="Generate SNP MSA and Binary Matrices from cellSNP results.")

    # Required Argument
    parser.add_argument(
        "--project",
        required=True,
        help="A valid scPhylogenomics project name. Script looks in results/<project>/"
    )

    # Base Calling Argument
    parser.add_argument(
        "--af_threshold",
        type=float,
        default=0.5,
        help="Allele Frequency (AF) threshold for calling a variant (AD/DP). Default: 0.5"
    )

    # New Filtering Arguments
    parser.add_argument(
        "--min_cells_per_snp",
        type=float,
        default=0.05,
        help="Filter SNPs: Keep SNP if variant is present in >= X fraction of cells. Default: 0.05 (5%)"
    )

    parser.add_argument(
        "--min_snps_per_cell",
        type=float,
        default=0.05,
        help="Filter Cells: Keep Cell if it has variants in >= X fraction of the remaining SNPs. Default: 0.05 (5%)"
    )

    return parser.parse_args()

def get_file_prefix(filename):
    """
    Extracts the 'sample_id.cell_type' prefix from a VCF filename.
    Expected format: <prefix>.cellSNP.base.vcf.gz
    """
    if ".cellSNP." in filename:
        return filename.split(".cellSNP.")[0]
    return None

def read_vcf_metadata(vcf_path):
    """Reads SNP metadata (REF, ALT) and column order from the VCF file."""
    print(f"  -> Reading SNP metadata from {os.path.basename(vcf_path)}...")

    try:
        vcf_df = pd.read_csv(
            vcf_path,
            sep='\t',
            comment='#',
            header=None,
            compression='gzip',
            usecols=[0, 1, 3, 4],  # CHROM, POS, REF, ALT
            names=['CHROM', 'POS', 'REF', 'ALT']
        )
    except Exception as e:
        print(f"Error reading VCF {vcf_path}: {e}")
        return None

    # Create a unique SNP identifier
    vcf_df['SNP_ID'] = vcf_df['CHROM'].astype(str) + '_' + vcf_df['POS'].astype(str) + '_' + vcf_df['REF'] + '_' + vcf_df['ALT']

    # The order in the VCF must match the 1-based index (row index) in the MTX files
    vcf_df['SNP_Index'] = vcf_df.index + 1

    return vcf_df.set_index('SNP_Index')

def read_samples(samples_path):
    """Reads ordered cell barcodes."""
    print(f"  -> Reading cell barcodes from {os.path.basename(samples_path)}...")
    if not os.path.exists(samples_path):
        print(f"Warning: Samples file not found: {samples_path}")
        return None

    barcodes_df = pd.read_csv(
        samples_path,
        sep='\t',
        header=None,
        names=['Barcode']
    )
    # Cell_Index maps to the MTX column index
    barcodes_df['Cell_Index'] = barcodes_df.index + 1
    return barcodes_df.set_index('Cell_Index')

def read_mtx_to_sparse_df(mtx_path, value_name):
    """Reads a sparse Matrix Market file into a Pandas DataFrame."""
    print(f"  -> Reading sparse matrix {os.path.basename(mtx_path)}...")
    if not os.path.exists(mtx_path):
        print(f"Warning: Matrix file not found: {mtx_path}")
        return None

    # Use delim_whitespace=True for Matrix Market files (spaces, not tabs)
    df = pd.read_csv(
        mtx_path,
        delim_whitespace=True,
        header=None,
        skiprows=3,
        names=['SNP_Index', 'Cell_Index', value_name]
    )

    # Explicitly convert indices to numeric to prevent string errors
    df['SNP_Index'] = pd.to_numeric(df['SNP_Index'], errors='coerce').astype('Int64')
    df['Cell_Index'] = pd.to_numeric(df['Cell_Index'], errors='coerce').astype('Int64')

    return df.dropna(subset=['SNP_Index', 'Cell_Index'])

def call_bases_and_create_matrices(vcf_df, samples_df, ad_df, dp_df, args):
    """
    1. Identifies SNPs where AF >= threshold in AT LEAST one cell.
    2. Generates the dense matrices.
    3. Filters SNPs (Columns) based on % of cells having the variant.
    4. Filters Cells (Rows) based on % of remaining SNPs having the variant.
    """
    af_threshold = args.af_threshold
    min_cells_fraction = args.min_cells_per_snp
    min_snps_fraction = args.min_snps_per_cell

    print("  -> Combining data and identifying informative SNPs...")

    # 1. Merge AD and DP
    combined_df = pd.merge(ad_df, dp_df, on=['SNP_Index', 'Cell_Index'], how='outer').fillna(0)

    # 2. Calculate Allele Frequency (AF)
    combined_df['AF'] = np.where(combined_df['DP'] > 0, combined_df['AD'] / combined_df['DP'], 0)

    # --- PRE-FILTER: Keep site if AF >= threshold in at least 1 cell ---
    informative_snps_mask = combined_df['AF'] >= af_threshold
    informative_snp_indices = combined_df.loc[informative_snps_mask, 'SNP_Index'].unique()
    informative_snp_indices = sorted(informative_snp_indices)

    num_total_snps = len(vcf_df)
    num_kept_snps = len(informative_snp_indices)

    if num_kept_snps == 0:
        print(f"     WARNING: No variants found with AF >= {af_threshold}. Output will be empty.")
        return None, None, None, None

    print(f"     Step 1 (Pre-filter): Keeping {num_kept_snps} sites with at least 1 variant call (out of {num_total_snps}).")

    # Filter VCF and Sparse Data
    filtered_vcf_df = vcf_df.loc[informative_snp_indices].copy()
    combined_df = combined_df[combined_df['SNP_Index'].isin(informative_snp_indices)]
    snp_index_map = {old_idx: new_idx for new_idx, old_idx in enumerate(informative_snp_indices)}

    # --- MATRIX GENERATION ---
    num_cells_total = len(samples_df)
    binary_matrix = np.full((num_cells_total, num_kept_snps), 3, dtype=int) # 3=Missing
    msa_matrix = np.full((num_cells_total, num_kept_snps), '-', dtype=str)

    print("  -> Populating initial dense matrices...")

    for row in combined_df.itertuples():
        old_snp_idx = int(row.SNP_Index)
        cell_idx = int(row.Cell_Index)
        af = row.AF
        dp = row.DP
        col_idx = snp_index_map[old_snp_idx]
        row_idx = cell_idx - 1

        if row_idx >= num_cells_total: continue
        if dp == 0: continue

        try:
            ref_base = filtered_vcf_df.at[old_snp_idx, 'REF']
            alt_base = filtered_vcf_df.at[old_snp_idx, 'ALT']
        except KeyError: continue

        if af >= af_threshold:
            binary_matrix[row_idx, col_idx] = 1
            msa_matrix[row_idx, col_idx] = alt_base
        else:
            binary_matrix[row_idx, col_idx] = 0
            msa_matrix[row_idx, col_idx] = ref_base

    # --- STEP 2: SNP FILTERING (Frequency in Population) ---
    print(f"  -> Filtering SNPs: Must be present in >= {min_cells_fraction*100}% of cells...")

    # Calculate how many cells have the variant (1) for each SNP (column)
    variant_counts_per_snp = np.sum(binary_matrix == 1, axis=0)

    # Calculate required count (rounding up)
    required_cell_count = math.ceil(num_cells_total * min_cells_fraction)

    # Create mask
    keep_snp_mask = variant_counts_per_snp >= required_cell_count

    num_snps_before = binary_matrix.shape[1]
    num_snps_after = np.sum(keep_snp_mask)
    print(f"     SNP Filter: Dropped {num_snps_before - num_snps_after} sites. Remaining: {num_snps_after} SNPs.")

    if num_snps_after == 0:
        print("     WARNING: All SNPs removed by prevalence filter.")
        return None, None, None, None

    # Apply SNP mask
    binary_matrix = binary_matrix[:, keep_snp_mask]
    msa_matrix = msa_matrix[:, keep_snp_mask]
    filtered_vcf_df = filtered_vcf_df.iloc[keep_snp_mask]

    # --- STEP 3: CELL FILTERING (Coverage of Variances) ---
    print(f"  -> Filtering Cells: Must have variants in >= {min_snps_fraction*100}% of remaining SNPs...")

    # Calculate how many variants (1) each cell (row) has
    variant_counts_per_cell = np.sum(binary_matrix == 1, axis=1)

    # Calculate required count based on REMAINING SNPs
    required_snp_count = math.ceil(num_snps_after * min_snps_fraction)

    # Create mask (Also ensures we remove cells with 0 variants)
    keep_cell_mask = variant_counts_per_cell >= required_snp_count

    num_cells_before = binary_matrix.shape[0]
    num_cells_after = np.sum(keep_cell_mask)
    print(f"     Cell Filter: Dropped {num_cells_before - num_cells_after} cells. Remaining: {num_cells_after} Cells.")

    if num_cells_after == 0:
        print("     WARNING: All Cells removed by variant coverage filter.")
        return None, None, None, None

    # Apply Cell mask
    binary_matrix = binary_matrix[keep_cell_mask, :]
    msa_matrix = msa_matrix[keep_cell_mask, :]
    filtered_samples_df = samples_df.iloc[keep_cell_mask].copy()

    return binary_matrix, msa_matrix, filtered_vcf_df, filtered_samples_df

def write_outputs(prefix, output_dir, vcf_df, samples_df, binary_matrix, msa_matrix):
    """Writes the FASTA, Binary Matrix, Filtered Sites, and Filtered Barcodes."""

    # Define Filenames
    fasta_filename = f"{prefix}.cellSNP.msa.fasta.gz"
    matrix_filename = f"{prefix}.cellSNP.snp.mtx.gz"
    sites_filename = f"{prefix}.cellSNP.sites.filtered.tsv.gz"
    barcodes_filename = f"{prefix}.cellSNP.barcodes.filtered.tsv.gz"

    # Define Paths
    fasta_path = os.path.join(output_dir, fasta_filename)
    matrix_path = os.path.join(output_dir, matrix_filename)
    sites_path = os.path.join(output_dir, sites_filename)
    barcodes_path = os.path.join(output_dir, barcodes_filename)

    # 1. Write MSA (FASTA)
    print(f"  -> Writing MSA to {fasta_filename}...")
    ref_sequence = "".join(vcf_df['REF'].tolist())
    cell_barcodes = samples_df['Barcode'].tolist()

    with gzip.open(fasta_path, 'wt') as f:
        f.write(">REFERENCE\n")
        f.write(ref_sequence + "\n")
        for i, barcode in enumerate(cell_barcodes):
            cell_sequence = "".join(msa_matrix[i, :])
            f.write(f">{barcode}\n")
            f.write(cell_sequence + "\n")

    # 2. Write Binary Matrix
    print(f"  -> Writing Binary Matrix to {matrix_filename}...")
    snp_ids = vcf_df['SNP_ID'].tolist()
    bin_df = pd.DataFrame(binary_matrix, index=cell_barcodes, columns=snp_ids)
    bin_df.to_csv(matrix_path, sep='\t', compression='gzip')

    # 3. Write Filtered Sites (CHROM, POS)
    print(f"  -> Writing Filtered SNP Sites to {sites_filename}...")
    # Select only CHROM and POS from the filtered VCF dataframe
    vcf_df[['CHROM', 'POS']].to_csv(
        sites_path,
        sep='\t',
        index=False,
        compression='gzip'
    )

    # 4. Write Filtered Barcodes (Index, Cell_type)
    print(f"  -> Writing Filtered Barcodes to {barcodes_filename}...")

    # Determine Cell_type from the prefix (value after the first dot)
    # Example: TNBC5.select-cell-types -> select-cell-types
    if '.' in prefix:
        cell_type_val = prefix.split('.', 1)[1]
    else:
        cell_type_val = prefix # Fallback if no dot

    # Create DataFrame for Barcodes output
    barcodes_out_df = pd.DataFrame(index=cell_barcodes)
    barcodes_out_df.index.name = 'Index'
    barcodes_out_df['Cell_type'] = cell_type_val

    barcodes_out_df.to_csv(
        barcodes_path,
        sep='\t',
        compression='gzip'
    )

def process_single_run(prefix, path_map, args, output_dir):
    """Orchestrates the processing for a single Sample+CellType combination."""
    print(f"\nProcessing group: {prefix}")

    vcf_df = read_vcf_metadata(path_map['vcf'])
    samples_df = read_samples(path_map['samples'])
    ad_df = read_mtx_to_sparse_df(path_map['ad'], 'AD')
    dp_df = read_mtx_to_sparse_df(path_map['dp'], 'DP')

    if any(x is None for x in [vcf_df, samples_df, ad_df, dp_df]):
        print(f"Skipping {prefix} due to missing input files.")
        return

    print(f"  -> Dimensions (Raw): {len(vcf_df)} SNPs, {len(samples_df)} Cells")

    # Generate Matrices with ROW and COLUMN FILTERING
    binary_mtx, msa_mtx, filtered_vcf, filtered_samples = call_bases_and_create_matrices(
        vcf_df, samples_df, ad_df, dp_df, args
    )

    if binary_mtx is None:
        print(f"Skipping writing for {prefix} (no informative data left).")
        return

    # Write Outputs using the FILTERED VCF and FILTERED SAMPLES
    write_outputs(prefix, output_dir, filtered_vcf, filtered_samples, binary_mtx, msa_mtx)

def main():
    args = parse_arguments()

    project_dir = os.path.join("results", args.project)

    if not os.path.exists(project_dir):
        print(f"Error: Project directory not found: {project_dir}")
        sys.exit(1)

    print(f"--- Starting Analysis for Project: {args.project} ---")
    print(f"--- AF Threshold: {args.af_threshold} ---")
    print(f"--- Filter: Min Cells per SNP: {args.min_cells_per_snp * 100}% ---")
    print(f"--- Filter: Min SNPs per Cell: {args.min_snps_per_cell * 100}% ---")

    vcf_files = glob.glob(os.path.join(project_dir, "*.cellSNP.base.vcf.gz"))

    if not vcf_files:
        print("No VCF files found in the project directory.")
        sys.exit(0)

    for vcf_path in vcf_files:
        filename = os.path.basename(vcf_path)
        prefix = get_file_prefix(filename)

        if not prefix: continue

        samples_path_gz = os.path.join(project_dir, f"{prefix}.cellSNP.samples.tsv.gz")
        samples_path_txt = os.path.join(project_dir, f"{prefix}.cellSNP.samples.tsv")
        samples_path = samples_path_gz if os.path.exists(samples_path_gz) else samples_path_txt

        path_map = {
            'vcf': vcf_path,
            'samples': samples_path,
            'ad': os.path.join(project_dir, f"{prefix}.cellSNP.tag.AD.mtx.gz"),
            'dp': os.path.join(project_dir, f"{prefix}.cellSNP.tag.DP.mtx.gz"),
        }

        if not os.path.exists(path_map['ad']): path_map['ad'] = path_map['ad'][:-3]
        if not os.path.exists(path_map['dp']): path_map['dp'] = path_map['dp'][:-3]

        process_single_run(prefix, path_map, args, project_dir)

    print("\nAll processing finished.")

if __name__ == "__main__":
    main()
