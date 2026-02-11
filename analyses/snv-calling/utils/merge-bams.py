import pysam
import argparse
import os

def merge_bams(input_bams, output_bam):
    """
    Merges multiple BAM files, sorts the result, and creates an index.
    """
    print(f"Merging {len(input_bams)} BAM files into {output_bam}...")

    # Merge
    pysam.merge("-f", output_bam, *input_bams)

    # Sort (pysam.sort creates a new file, so we overwrite the merged one)
    tmp_sort = output_bam + ".tmp"
    pysam.sort("-o", tmp_sort, output_bam)
    os.replace(tmp_sort, output_bam)

    # Index
    print(f"Indexing merged BAM: {output_bam}")
    pysam.index(output_bam)

def main():
    parser = argparse.ArgumentParser(description='Merge and index multiple BAM files')
    parser.add_argument('--inputs', nargs='+', required=True, help='List of input BAM files')
    parser.add_argument('--output', type=str, required=True, help='Output merged BAM path')
    args = parser.parse_args()

    merge_bams(args.inputs, args.output)

if __name__ == '__main__':
    main()
