#!/usr/bin/env python3
# Eric Wafula
# 2025

import os
import argparse
import subprocess
import logging
import shutil
import gzip
from pathlib import Path

# Define script location relative to execution
home = Path(__file__).resolve().parent

# =============================================================================
# CONSTANTS & CONFIGURATION
# =============================================================================

# Relative paths
SNV_CALLING_DIR = Path("../snv-calling/results")
MODULE_RESULTS_DIR = Path("results")
SCRATCH_DIR = Path("../../scratch")

def setup_logger(log_file):
    """Setup logging to console and file."""
    # Ensure log directory exists
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def setup_directories(path):
    if not os.path.exists(path):
        os.makedirs(path)

def decompress_file(input_path, output_path):
    """Decompresses a gzipped file to a target output path."""
    with gzip.open(input_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def infer_tree(input_file, output_prefix, method, threads):
    """
    Infer phylogenetic tree using FastTree or IQ-TREE.
    Always decompresses gzipped input to a temporary file for tool consistency.
    """
    temp_input = None
    actual_input = input_file

    try:
        # --- Standardized Decompression ---
        # If input is gzipped, decompress it first regardless of the method.
        if input_file.suffix == '.gz':
            logging.info(f"Input is gzipped. Decompressing...")
            
            # Define temp file path (original name without .gz)
            # Located in the scratch directory alongside the output
            temp_input = output_prefix.parent / input_file.stem 
            
            decompress_file(input_file, temp_input)
            actual_input = temp_input
            logging.info(f"Temporary decompressed file created: {actual_input}")

        # --- Construct Command ---
        cmd = []
        if method == "fasttree":
            tree_file = f"{output_prefix}.tree"
            cmd = [
                "FastTree", 
                "-nt", 
                "-gtr", 
                "-gamma",
                "-out", tree_file,
                str(actual_input)
            ]
            
        elif method == "iqtree":
            cmd = [
                "iqtree3",
                "-s", str(actual_input),
                "-st", "DNA",
                "-fast",
                "-mset", "JC,K80,HKY,GTR",
                "-m", "MF+G+I",
                "-nt", str(threads),
                "-pre", str(output_prefix),
                "-redo"
            ]
        else:
            raise ValueError(f"Unsupported method: {method}")

        # --- Execute ---
        logging.info(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        logging.info(f"Inferred tree using {method}: {output_prefix}")

    finally:
        # --- Cleanup ---
        # Remove the temporary decompressed file to save space
        if temp_input and temp_input.exists():
            try:
                temp_input.unlink()
                logging.info(f"Removed temporary decompressed file: {temp_input}")
            except Exception as e:
                logging.warning(f"Failed to remove temp file {temp_input}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Phylogeny Inference Module (scPhylogenomics)")

    # Required Arguments
    parser.add_argument("--project", required=True, help="Project ID (e.g., AML)")
    parser.add_argument("--sample", required=True, help="Sample ID (e.g., LE1)")
    parser.add_argument("--cell_type", required=True, help="Cell type string (e.g., all-cell-types)")

    # Optional Arguments
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of parallel threads (default = 4).")
    parser.add_argument("--method", choices=["iqtree", "fasttree"], default="iqtree", help="Tree inference method (default = iqtree).")
    parser.add_argument("--log", default="phylogeny_inference.log", help="Log file name.")

    args = parser.parse_args()

    # 1. Define Directories
    # Input: ../snv-calling/results/{project}/
    input_dir = SNV_CALLING_DIR / args.project
    
    # Scratch Output: ../../scratch/phylogeny/{project}/{sample}/
    scratch_output_dir = SCRATCH_DIR / "phylogeny" / args.project / args.sample
    
    # Module Results: results/{project}/
    module_results_dir = MODULE_RESULTS_DIR / args.project

    setup_directories(scratch_output_dir)
    setup_directories(module_results_dir)

    # 2. Setup Logging
    setup_logger(scratch_output_dir / args.log)

    # 3. Locate Input File
    msa_filename = f"{args.sample}.{args.cell_type}.cellSNP.msa.fasta.gz"
    input_msa = input_dir / msa_filename

    if not input_msa.exists():
        logging.error(f"Input MSA file not found: {input_msa}")
        return

    logging.info(f"Processing Project: {args.project} | Sample: {args.sample} | Cell Type: {args.cell_type}")
    logging.info(f"Input File: {input_msa}")
    logging.info(f"Output Scratch Directory: {scratch_output_dir}")

    # 4. Define Output Prefix
    output_base_name = f"{args.sample}.{args.cell_type}.{args.method}"
    output_prefix = scratch_output_dir / output_base_name

    # 5. Run Tree Inference
    try:
        infer_tree(input_msa, output_prefix, args.method, args.threads)
    except subprocess.CalledProcessError as e:
        logging.error(f"Tree inference failed for {msa_filename}: {e}")
        return
    except Exception as ex:
        logging.exception(f"An unexpected error occurred during inference.")
        return

    # 6. Copy Final Tree to Results Directory
    if args.method == "iqtree":
        # IQ-TREE usually creates .treefile
        source_tree = scratch_output_dir / f"{output_base_name}.treefile"
        dest_tree_name = f"{args.sample}.{args.cell_type}.iqtree.tree"
    else:
        # FastTree creates .tree directly
        source_tree = scratch_output_dir / f"{output_base_name}.tree"
        dest_tree_name = f"{args.sample}.{args.cell_type}.fasttree.tree"

    dest_tree_path = module_results_dir / dest_tree_name

    if source_tree.exists():
        logging.info(f"Copying final tree to results directory...")
        try:
            shutil.copy(source_tree, dest_tree_path)
            logging.info(f"Successfully created: {dest_tree_path}")
        except Exception as e:
            logging.error(f"Failed to copy tree file: {e}")
    else:
        logging.error(f"Expected output tree file not found at: {source_tree}")

if __name__ == "__main__":
    main()
    
