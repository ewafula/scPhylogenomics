#!/usr/bin/env python3
# Eric Wafula

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
parser.add_argument("--editing", required=False, help="Path to RNA editing sites file. Optional.")
parser.add_argument("--pon", required=False, help="Path to Panel of Normals (PoN) file. Optional.")

args = parser.parse_args()

project = args.project
num_threads = args.num_threads
minMAF = args.minMAF
minCOUNT = args.minCOUNT

editing_file = args.editing if args.editing else ""
pon_file = args.pon if args.pon else ""

print()
print(f"{datetime.now()} - Starting SNV calling for {project} project")
print("Scanning project directory for cell types...")

# Match your actual path: inputs/SPECTRUM/
inputs_dir = home / "inputs" / project

# 1. Identify all unique cell types from the flat file list
all_cell_types = set()
# Search for the specific naming pattern in the flat directory
barcode_files = glob.glob(f"{inputs_dir}/*-cancer-cells-barcodes.tsv.gz")

if not barcode_files:
    print(f"[ERROR] No barcode files found in {inputs_dir}")
    print(f"Expected files matching pattern: *-cancer-cells-barcodes.tsv.gz")
    sys.exit(1)

for bc_file in barcode_files:
    try:
        with gzip.open(bc_file, 'rt') as f:
            for line in f:
                if line.startswith("Index"):
                    continue
                fields = line.strip().split('\t')
                if len(fields) > 1:
                    all_cell_types.add(fields[1])
    except Exception as e:
        print(f"[WARNING] Could not read {bc_file}: {e}")

print(f"Found {len(all_cell_types)} unique cell types: {', '.join(sorted(all_cell_types))}\n")

# 2. Prepare commands: One call per Cell Type
commands = []
for cell_type in sorted(all_cell_types):
    cmd = [
        "bash", "utils/run-cellsnp-lite.sh",
        project,
        cell_type,
        str(num_threads),
        str(minMAF),
        str(minCOUNT),
        editing_file,
        pon_file
    ]
    commands.append(cmd)

# 3. Parallelization
def run_command(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    return {
        "cell_type": cmd[2], 
        "cmd": " ".join(cmd),
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr
    }

with ThreadPoolExecutor(max_workers=max(1, num_threads // 2)) as executor:
    futures = {executor.submit(run_command, cmd): cmd for cmd in commands}

    for future in as_completed(futures):
        result = future.result()
        if result["returncode"] != 0:
            print(f"[ERROR] Workflow failed for cell type: {result['cell_type']}")
            print(result['stderr'])
        else:
            print(f"[OK] Finished analysis for cell type: {result['cell_type']}")

print(f"\n{datetime.now()} - Completed project {project} SNV calling pipeline.\n")
