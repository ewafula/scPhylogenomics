# SNV Calling

## Purpose
This module calls and characterizes somatic single nucleotide variants (SNVs) across cell populations for tumor samples at single-cell resolution. The workflow begins by extracting barcode identifiers for cell types of interest to ensure only relevant cells are analyzed. It then automates sample- and cell type-specific SNV calling using [cellsnp-lite](https://academic.oup.com/bioinformatics/article/37/23/4569/6272512), efficiently processing sequencing data to detect SNVs and generate generate genotype matrices. The pipeline includes optional filtering using RNA Editing sites and Panels of Normals (PoN) if available. Finally, the results are processed to generate mutation matrices and pseudo-multiple sequence alignments (MSAs) suitable for phylogenomic and clonal clustering analyses.

## Analysis
This module performs SNV calling to produce mutation matrices and MSAs suitable for phylogenomic analyses.

#### Module directory structure
Illustration of module directory structure based on the current implementation.

```
.
|-- 01-get-cell-type-barcodes.R
|-- 02-snv-calling.py
|-- 03-generate-snp-msa.py
|-- run-snv-calling.sh
|-- README.md
|-- inputs
|   |-- AllEditingSites.hg38.txt
|   |-- PoN.scRNAseq.hg38.tsv
|   |-- AML
|   |   `-- LE1-cancer-cells-barcodes.tsv.gz
|   `-- TNBC
|       `-- TNBC5-cancer-cells-barcodes.tsv.gz
|-- results
|   |-- AML
|   |   |-- LE1.all-cell-types.cellSNP.barcodes.filtered.tsv.gz
|   |   |-- LE1.all-cell-types.cellSNP.base.vcf.gz
|   |   |-- LE1.all-cell-types.cellSNP.msa.fasta.gz
|   |   |-- LE1.all-cell-types.cellSNP.samples.tsv.gz
|   |   |-- LE1.all-cell-types.cellSNP.sites.filtered.tsv.gz
|   |   |-- LE1.all-cell-types.cellSNP.snp.mtx.gz
|   |   |-- LE1.all-cell-types.cellSNP.tag.AD.mtx.gz
|   |   |-- LE1.all-cell-types.cellSNP.tag.DP.mtx.gz
|   |   `-- LE1.all-cell-types.cellSNP.tag.OTH.mtx.gz
|   `-- TNBC
|       |-- TNBC5.select-cell-types.cellSNP.barcodes.filtered.tsv.gz
|       |-- TNBC5.select-cell-types.ellSNP.base.vcf.gz
|       |-- TNBC5.select-cell-types.cellSNP.msa.fasta.gz
|       |-- TNBC5.select-cell-types.cellSNP.samples.tsv.gz
|       |-- TNBC5.select-cell-types.cellSNP.sites.filtered.tsv.gz
|       |-- TNBC5.select-cell-types.cellSNP.snp.mtx.gz
|       |-- TNBC5.select-cell-types.cellSNP.tag.AD.mtx.gz
|       |-- TNBC5.select-cell-types.cellSNP.tag.DP.mtx.gz
|       `-- TNBC5.select-cell-types.cellSNP.tag.OTH.mtx.gz
`-- utils
    |-- csp_utils.py
    |-- run-cellsnp-lite.sh
    |-- split-bam.py
    `-- variant-calling.py
```

#### `run-snv-calling.sh`
This script serves as a module wrapper script to automate project SNV calling tasks and can be adapted for any project. All specified paths in this script are relative to the module directory. Therefore, it should always be executed relative to this [module directory](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/snv-calling). 

###### Example usage:
```bash
bash run-snv-calling.sh
```

#### `01-get-cell-type-barcodes.R`
This script retrieves all barcode identifiers for selected cell types and/or those classified as tumor cells (aneuploid) by the [ploidy-inference module](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/ploidy-inference).  These barcodes are then saved in a suitable format for sample-specific SNV calling using the [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/main/manual.html), facilitating efficient detection of somatic mutations across various cell populations.

###### Example usage:
```
Rscript --vanilla 01-get-cell-type-barcodes.R \
  --project AML \
  --cell_types all-cell-types
```

###### Argument descriptions:
```
Usage: 01-get-cell-type-barcodes.R [options]


Options:
	--project=CHARACTER
		A valid scPhylogenomics project name

	--cell_types=CHARACTER
		Cell types to used SNV calling - [default "aneuploid"]
                      - Either comma-separated list of cell type names or a single name,
                      - cell type names in quotations if includes white space(s)
                      - if set to `all`, all cells will be used
                      - if set to `aneuploid`, cancer cells determined using copy 
                      - number changes will be used - dafault

	-h, --help
		Show this help message and exit
```

#### `02-snv-calling.py`
This script automates SNV calling for scPhylogenomics projects using the [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/main/manual.html) bash workflow ([utils/run-cellsnp-lite.sh](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/snv-calling/utils/run-cellsnp-lite.sh)). It scans the project input directories for barcode files associated with seleted cells, identifies unique cell types, and launches jobs in parallel. Each job runs the `cellsnp-lite workflow`, which activates the `cellsnp-lite conda environment` and processes data through multiple steps, including, splitting BAM files by selcted cell types, calling variants using cellsnp-lite, and pptionally filtering variants using a Panel of Normals (PoN) and RNA Editing sites if provided. Results are stored in a structured output directory for the subsequent creation of SNV psuedo-multiple sequence alignments (MSA). 

If the analysis is conducted outside the scPhylogenomics Docker environment, the location of the cellsnp-lite conda environment installation must be updated in the bash workflow script ([utils/run-cellsnp-lite.sh](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/snv-calling/utils/run-cellsnp-lite.sh)).

###### Example usage:
```
python 02-snv-calling.py \
  --project AML\
  --num_threads 16
```

###### Argument descriptions:
```
usage: 02-snv-calling.py [-h] --project PROJECT [--num_threads NUM_THREADS] [--minMAF MINMAF] [--minCOUNT MINCOUNT] [--editing EDITING]
                         [--pon PON]

scPhylogenomics SNV Calling

optional arguments:
  -h, --help            show this help message and exit
  --project PROJECT     A valid scPhylogenomics project name
  --num_threads NUM_THREADS
                        Number of threads (default: 4)
  --minMAF MINMAF       Minimum MAF. Default 0.1.
  --minCOUNT MINCOUNT   Minimum COUNT. Default 100.
  --editing EDITING     Path to RNA editing sites file (columns: chr, pos). Optional.
  --pon PON             Path to Panel of Normals (PoN) file. Optional.
```

#### `03-generate-snp-msa.py`
This script processes the raw or filtered cellsnp-lite outputs to generate the final SNP Multiple Sequence Alignment (MSA) and a binary mutation profile matrix. It applies strict filtering to ensure only informative SNPs and high-quality cells are retained for downstream phylogenomic and clonal clustering analyses. It produces SNP psuedo-multiple sequence alignments (MSA) per cell, a binary mutation profile matrix coded with 0/1/3(Ref/Alt/Missing, and filtered VCF, sites, cell barcode for list of SNPs retained after filtering.


###### Example usage:
```
python 03-generate-snp-msa.py \
  --project SPECTRUM \
  --af_threshold 0.5 \
  --min_cells_per_snp 0.05 \
  --min_snps_per_cell 0.05
```

###### Argument descriptions:

```
usage: 03-generate-snp-msa.py [-h] --project PROJECT [--af_threshold AF_THRESHOLD]
                              [--min_cells_per_snp MIN_CELLS_PER_SNP]
                              [--min_snps_per_cell MIN_SNPS_PER_CELL]

Generate SNP MSA and Binary Matrices from cellSNP results.

options:
  -h, --help            show this help message and exit
  --project PROJECT     A valid scPhylogenomics project name.
  --af_threshold AF_THRESHOLD
                        Allele Frequency (AF) threshold for calling a variant. Default: 0.5
  --min_cells_per_snp MIN_CELLS_PER_SNP
                        Filter SNPs: Keep SNP if variant is present in >= X fraction of cells. Default: 0.05 (5%)
  --min_snps_per_cell MIN_SNPS_PER_CELL
                        Filter Cells: Keep Cell if it has variants in >= X fraction of the remaining SNPs. Default: 0.05 (5%)
```
