import pysam
import pandas as pd
import argparse
import timeit
import sys
import numpy as np

def meta_to_dict(txt, tissue):
    """
    Reads the metadata file, cleans the cell barcodes, and prepares a dictionary
    mapping clean barcodes to cell types.

    Args:
        txt (str): Path to the metadata file.
        tissue (str or None): Optional tissue ID to prefix to the cell type name.

    Returns:
        tuple: (DICT, ALL_CELL_TYPES) where DICT is the barcode-to-cell_type mapping
            and ALL_CELL_TYPES is a list of unique cell types.
    """
    metadata = pd.read_csv(txt, delimiter="\t")

    # Clean index column (remove the '-1' suffix from 10x barcodes)
    metadata['Index_clean'] = metadata['Index'].str.replace('-.*$', '', regex=True)

    # If tissue provided, append tissue id to cell type to be recognised in downstream analysis
    if tissue is None:
        metadata['Cell_type_clean'] = metadata['Cell_type'].str.replace(' ', '_', regex=True)
    else:
        # Prefix cell type with tissue ID for uniqueness
        tissue = tissue.replace(" ", "_")
        metadata['Cell_type_clean'] = metadata['Cell_type'].str.replace(' ', '_', regex=True)
        metadata['Cell_type_clean'] = str(tissue) + '__' + metadata['Cell_type_clean'].astype(str)

    # Create dictionary with cell types and cell barcodes (using the cleaned barcode as key)
    DICT = metadata.set_index('Index_clean')['Cell_type_clean'].to_dict()
    ALL_CELL_TYPES = metadata['Cell_type_clean'].unique()

    del metadata

    return (DICT, ALL_CELL_TYPES)


def split_bam(bam, txt, outdir, donor, tissue, max_NM, max_NH, min_MAPQ, n_trim):
    """
    Splits the input BAM file into cell-type specific BAMs, applies filtering,
    prefixes cell barcodes with the donor/sample ID, and saves the new barcodes.

    Args:
        bam (str): Path to the input BAM file.
        txt (str): Path to the metadata file.
        outdir (str): Output directory path.
        donor (str): Sample ID (used as barcode prefix).
        tissue (str or None): Tissue ID (optional, used in meta_to_dict).
        max_NM (int or None): Maximum number of mismatches.
        max_NH (int or None): Maximum number of alignment hits.
        min_MAPQ (int): Minimum mapping quality.
        n_trim (int): Number of bases to trim (set quality to 0).
    """

    start = timeit.default_timer()

    # 1. Transform table to dictionary
    DICT, ALL_CELL_TYPES = meta_to_dict(txt, tissue)

    if len(DICT) < 1:
        print('Warning: No cell barcodes found in the --meta file')
        sys.exit()

    # 2. Open infile
    infile = pysam.Samfile(bam, "rb")

    # 3. Create and open out bam files and initialize barcode storage
    DICT_files = {}
    # 💡 NEW: Dictionary to store unique, prefixed barcodes for each cell type
    DICT_barcodes = {}
    for cell_type in ALL_CELL_TYPES:
        outfile = "{}/{}.{}.bam".format(outdir, donor, cell_type)
        outfile_wb = pysam.AlignmentFile(outfile, "wb", template=infile)
        DICT_files[cell_type] = outfile_wb
        # Use a set to automatically handle unique barcodes for each cell type
        DICT_barcodes[cell_type] = set()

    # 4. Start read counts
    total_reads = 0
    FILTER_dict = {'Total_reads': 0, 'Pass_reads': 0, 'CB_not_found': 0, 'CB_not_matched': 0} # To store filter reasons
    
    # 5. Check reads and split them in bam files
    for read in infile.fetch():
        FILTER_dict['Total_reads'] += 1 # To count the total number of reads analysed
        total_reads += 1

        if (total_reads % 5000000 == 0):
            print('Number of reads already processed: ' + str(total_reads))

        # Check if CB tag is present
        try:
            barcode_with_suffix = read.opt("CB") # e.g., AAACCCAAGACCAACG-1
        except:
            FILTER2 = 'CB_not_found'
            FILTER_dict[FILTER2] += 1
            continue

        # Extract the core barcode for lookup (e.g., AAACCCAAGACCAACG)
        barcode_no_suffix = barcode_with_suffix.split("-")[0]
        
        # Check if CB code matches with an annotated cell type
        try:
            CELL_TYPE = DICT[barcode_no_suffix]
        except:
            FILTER2 = 'CB_not_matched'
            FILTER_dict[FILTER2] += 1
            continue

        # NEW: Create the new prefixed barcode (e.g., Le1_AAACCCAAGACCAACG-1)
        new_barcode = f"{donor}_{barcode_with_suffix}"

        # Final filters
        FILTER = []
        # Check number of mismatches in the read
        if (max_NM != None):
            try:
                if (read.opt("nM") > max_NM):
                    FILTER.append('nM')
            except KeyError as e:
                FILTER.append('nM_not_found')
        # Check number of hits
        if (max_NH != None):
            try:
                if (read.opt("NH") > max_NH):
                    FILTER.append('NH')
            except KeyError as e:
                FILTER.append('NH_not_found')

        # Check mapping quality
        if (min_MAPQ > 0):
            try:
                if (read.mapq < min_MAPQ):
                    FILTER.append('MAPQ')
            except:
                FILTER.append('MAPQ_not_found')

        # Making a decision about the read (filtered or not)
        if (len(FILTER) > 0): # If there are reasons to filter read
            FILTER2 = ';'.join(FILTER)
            
            FILTER_dict[FILTER2] = FILTER_dict.get(FILTER2, 0) + 1
            continue
        else:
            FILTER_dict['Pass_reads'] += 1

        # Only for PASS reads
        # Trim last and first bases of the read (reduce quality) if specified
        if (n_trim > 0):
            CIGAR_tuple = read.cigartuples
            if (CIGAR_tuple != None and len(CIGAR_tuple) > 1):
                # Check the number of bases to trim at the beginning of the read
                if (CIGAR_tuple[0][0] == 4):
                    # Conservative trimming logic
                    if (CIGAR_tuple[0][1] >= 20 and CIGAR_tuple[0][1] < 30):
                        trim_start = 30 + n_trim
                    else:
                        trim_start = CIGAR_tuple[0][1] + n_trim
                else:
                    trim_start = n_trim

                # Check the number of bases to trim at the end of the read
                if (CIGAR_tuple[-1][0] == 4):   
                    if (CIGAR_tuple[-1][1] >= 20 and CIGAR_tuple[-1][1] < 30):
                        trim_end = 30 + n_trim
                    else:
                        trim_end = CIGAR_tuple[-1][1] + n_trim
                else:
                    trim_end = n_trim
            else:
                trim_start = n_trim
                trim_end = n_trim

            # Set base qualities to 0 at the beginning and the end of the read if specified
            # Get first and last qualities
            Q_indexes = list(np.arange(trim_start)) + list((np.arange(trim_end) + 1)*-1)
            Q = read.query_qualities
            for Qi in Q_indexes:
                Q[Qi] = 0

            # Substitute qualities
            read.query_qualities = Q
        
        # MODIFIED: Update the 'CB' tag in the read with the new prefixed barcode
        read.set_tag('CB', new_barcode, value_type='Z', replace=True)
        
        # NEW: Store the unique prefixed barcode for later saving
        DICT_barcodes[CELL_TYPE].add(new_barcode)

        # Print passed read to the cell-type specific BAM file
        DICT_files[CELL_TYPE].write(read)

    # 6. Close opened files
    # (Indent is 4 spaces - outside the read loop)
    infile.close()
    for key in DICT_files.keys():
        DICT_files[key].close()
        
    # Write Barcode files for each cell type
    print('Writing barcode files...')
    for cell_type, barcodes in DICT_barcodes.items():
        # Format the output file name as requested: {outdir}/{donor}.{cell_type}.barcodes.tsv
        outfile_barcodes = "{}/{}.{}.barcodes.tsv".format(outdir, donor, cell_type)
        with open(outfile_barcodes, 'w') as fout:
            # Write each unique, prefixed barcode on a new line
            for barcode in sorted(list(barcodes)):
                fout.write(f"{barcode}\n")
    print('Barcode files successfully written.')

    # 7. Get and create report
    outfile_report = "{}/{}.report.txt".format(outdir, donor)
    stop = timeit.default_timer()
    endtime = round((stop - start), 2)
    FILTER_dict['Total_time'] = endtime

    data_df = pd.DataFrame([FILTER_dict])
    data_df.to_csv(outfile_report, index=False, sep='\t')

    # 8. Index bam files
    print('Indexing BAM files...')
    for cell_type in ALL_CELL_TYPES:
        bam = "{}/{}.{}.bam".format(outdir, donor, cell_type)
        pysam.index(bam)
    print('BAM indexing complete.')


def initialize_parser():
    parser = argparse.ArgumentParser(description='Split alignment file into cell type specific BAMs')
    parser.add_argument('--bam', type=str, default=1, help='BAM file to be analysed (Sorted by coordinate)', required=True)
    parser.add_argument('--meta', type=str, default=1, help='Metadata file mapping cell barcodes to cell type information', required=True)
    parser.add_argument('--id', type=str, default='Sample', help='Sample ID', required=False)
    parser.add_argument('--max_nM', type=int, default=None, help='Maximum number of mismatches permitted. [Default: Switched off]', required=False)
    parser.add_argument('--max_NH', type=int, default=None, help='Maximum number of alignment hits permitted. [Default: Switched off]', required=False)
    parser.add_argument('--min_MQ', type=int, default=255, help='Minimum mapping quality required. [Default: 255]', required=False)
    parser.add_argument('--n_trim', type=int, default=0, help='Number of bases trimmed by setting the base quality to 0 [Default: 0]', required=False)
    parser.add_argument('--outdir', default='.', help='Out directory', required=False)
    return (parser)


def main():
    # 1. Arguments
    parser = initialize_parser()
    args = parser.parse_args()

    bam = args.bam
    outdir = args.outdir
    txt = args.meta
    donor = args.id
    tissue = None
    max_NM = args.max_nM
    min_MAPQ = args.min_MQ
    n_trim = args.n_trim
    max_NH = args.max_NH

    # 2. Split bam file
    split_bam(bam, txt, outdir, donor, tissue, max_NM, max_NH, min_MAPQ, n_trim)


#-------------
# Execute code
#-------------

if __name__ == '__main__':
    main()
