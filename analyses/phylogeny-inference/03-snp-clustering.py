#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
import subprocess
from pathlib import Path

# Add the directory containing 'utils' to sys.path
sys.path.append(os.path.dirname(__file__))

# Import logic from separate modules
# We wrap this in try/except because manifold_utils might fail to import
# if the main environment lacks torch/SNPmanifold, which is expected.
try:
    from utils import hierarchical_utils
except ImportError as e:
    print(f"Error importing utility modules: {e}")
    sys.exit(1)

# --- CONFIGURATION ---

# Detect environment location (Docker vs Local)
if os.path.exists("/opt/conda/envs/snpmanifold_env/bin/python"):
    CONDA_ENV_PYTHON = "/opt/conda/envs/snpmanifold_env/bin/python"
elif os.path.exists(os.path.expanduser("~/miniconda3/envs/snpmanifold_env/bin/python")):
    CONDA_ENV_PYTHON = os.path.expanduser("~/miniconda3/envs/snpmanifold_env/bin/python")
else:
    print("Error: Could not locate snpmanifold_env python executable.")
    sys.exit(1)

# Paths relative to analyses/phylogeny-inference/
SCRATCH_DIR = Path("../../scratch")
MODULE_RESULTS_DIR = Path("results")
MODULE_PLOTS_DIR = Path("plots")
SNV_CALLING_DIR = Path("../snv-calling/results")

def setup_directories(path):
    if not os.path.exists(path):
        os.makedirs(path)

def run_manifold_subprocess(args, scratch_output_dir):
    """
    Constructs a command to run utils/manifold_utils.py using the
    snpmanifold_env Python executable.
    """
    # Path to the worker script (now named manifold_utils.py)
    task_script = Path(os.path.dirname(__file__)) / "utils" / "manifold_utils.py"

    if not task_script.exists():
        raise FileNotFoundError(f"Could not find worker script: {task_script}")

    if not os.path.exists(CONDA_ENV_PYTHON):
        raise FileNotFoundError(f"Could not find Conda Python at: {CONDA_ENV_PYTHON}\nPlease check the path in 03-snp-clustering.py")

    # Build command: [python_path, script_path, --arg, val, ...]
    cmd = [
        str(CONDA_ENV_PYTHON),
        str(task_script),
        "--project", args.project,
        "--sample", args.sample,
        "--cell_type", args.cell_type,
        "--max_cluster", str(args.max_cluster),
        "--scratch_dir", str(SCRATCH_DIR),
        # Pass the pre-calculated specific output dir to ensure consistency
        "--output_dir", str(scratch_output_dir)
    ]

    if args.final:
        cmd.append("--final")

    print(f"--- Dispatching Manifold Task to Environment: {CONDA_ENV_PYTHON} ---")
    print(f"Command: {' '.join(cmd)}")

    # Run the subprocess and wait for it to finish
    try:
        subprocess.check_call(cmd)
        print("--- Manifold Task Completed Successfully ---")
    except subprocess.CalledProcessError as e:
        print(f"--- Manifold Task Failed with Exit Code {e.returncode} ---")
        sys.exit(e.returncode)

def main():
    parser = argparse.ArgumentParser(
        description="Unified SNP Clustering for scPhylogenomics (Manifold or Hierarchical)"
    )

    # Required Arguments
    parser.add_argument("--project", required=True, help="Project ID (e.g., AML)")
    parser.add_argument("--sample", required=True, help="Sample ID (e.g., LE1)")
    parser.add_argument("--cell_type", required=True, help="Cell type string (e.g., all-cell-types)")

    # Core Optional Arguments
    parser.add_argument("--method", choices=['manifold', 'hierarchical'], default='manifold',
                        help="Clustering method to use.")
    parser.add_argument("--max_cluster", type=int, default=10,
                        help="Maximum number of clusters (Manifold) or Max k for narrow search (Hierarchical).")
    parser.add_argument("--final", action="store_true",
                        help="Flag to indicate final run; copies outputs to module results/ and plots/.")

    # Hierarchical-Specific Tuning Arguments
    parser.add_argument("--force_k", type=int, help="[Hierarchical] Force specific narrow clusters per clade.")
    parser.add_argument("--broad_k", type=int, default=4, help="[Hierarchical] Number of broad hierarchical clades (default = 4).")
    parser.add_argument("--min_k", type=int, default=2, help="[Hierarchical] Min k for narrow cluster search (default = 2).")

    args = parser.parse_args()

    # Determine Output Directory
    stage = "final" if args.final else "initial"
    scratch_sub = "snpmanifold" if args.method == "manifold" else "hierarchical"

    # Path: ../../scratch/{method}/{project}/{sample}/{stage}/
    scratch_output_dir = SCRATCH_DIR / scratch_sub / args.project / args.sample / stage
    setup_directories(scratch_output_dir)

    print(f"Output directory: {scratch_output_dir}")

    # Dispatch logic
    if args.method == 'manifold':
        run_manifold_subprocess(args, scratch_output_dir)
    elif args.method == 'hierarchical':
        hierarchical_utils.run_hierarchical(
            args,
            scratch_output_dir,
            SNV_CALLING_DIR,
            MODULE_RESULTS_DIR,
            MODULE_PLOTS_DIR
        )

if __name__ == "__main__":
    main()
