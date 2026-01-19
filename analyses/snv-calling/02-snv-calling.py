#!/usr/bin/env python3
# Eric Wafula
# 2025

import os
import sys
import glob
import gzip
import subprocess
import argparse
from datetime import datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

home = Path(__file__).resolve().parent

# Argument parsing
parser = argparse.ArgumentParser(
    description="scPhylogenomics SNV Calling",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument('--project', required=True, help='A valid scPhylogenomics project name')
parser.add_argument('--num_threads', type=int, default=4, help='Number of threads (default: 4)')
parser.add_argument("--minMAF", type=float, default=0.1, help="Minimum MAF. Default 0.1.")
parser.add_argument("--minCOUNT", type=int, default=100, help="Minimum COUNT. Default 100.")
parser.add_argument("--editing", required=False, help="Path to RNA editing sites file (columns: chr, pos). Optional.")
parser.add_argument("--pon", required=False, help="Path to Panel of Normals (PoN) file. Optional.")

args = parser.parse_args()

project = args.project
num_threads = args.num_threads
minMAF = args.minMAF
minCOUNT = args.minCOUNT

# Initialize variables with empty strings if arguments are not provided.
# This ensures we pass *something* to the bash script (an empty string) to keep positional arguments aligned.
editing_file = args.editing if args.editing else ""
pon_file = args.pon if args.pon else ""

print()
print(f"{datetime.now()} - Starting SNV calling for {project} project\n")

inputs_dir = home / "inputs" / project
commands = []

# Collect all commands to be run
for cell_barcodes in glob.glob(f"{inputs_dir}/*.gz"):
    sample = os.path.basename(cell_barcodes)
    if "-cancer-cells-barcodes.tsv.gz" in sample:
        sample = sample.split("-cancer-cells-barcodes.tsv.gz")[0]
    else:
        continue

    cell_types = set()
    with gzip.open(cell_barcodes, 'rt') as f:
        for line in f:
            if line.startswith("Index"):
                continue
            fields = line.strip().split('\t')
            if len(fields) > 1:
                cell_types.add(fields[1])

    # Parallelization for cell type currently not utilized in scPhylogenomics, 
    # but would if multi-cell-specific SComatic analysis performed.
    for cell_type in sorted(cell_types): 
        cmd = [
            "bash", "utils/run-cellsnp-lite.sh",
            project, 
            sample, 
            cell_barcodes, 
            cell_type, 
            str(num_threads), 
            str(minMAF), 
            str(minCOUNT),
            editing_file,
            pon_file
        ]
        commands.append(cmd)

# Run commands in parallel
def run_command(cmd):
    # Pass the command list directly; subprocess handles empty strings correctly as arguments
    result = subprocess.run(cmd, capture_output=True, text=True)
    return {
        "cmd": " ".join(cmd),
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr
    }

with ThreadPoolExecutor(max_workers=num_threads) as executor:
    futures = {executor.submit(run_command, cmd): cmd for cmd in commands}

    for future in as_completed(futures):
        result = future.result()
        if result["returncode"] != 0:
            print(f"[ERROR] Command failed: {result['cmd']}")
            print(result['stderr'])
        else:
            print(f"[OK] Finished: {result['cmd']}")

print(f"\n{datetime.now()} - Completed SNV calling for {project} project\n")
