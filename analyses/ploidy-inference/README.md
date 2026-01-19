# Ploidy Inference

## Purpose

This module infers ploidy for annotated single cells from the [cell-typing module](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/cell-typing) and classifies cells by malignancy status. Cancer cells often exhibit genome-wide aneuploidy, characterized by abnormal chromosome numbers, while normal stromal and immune cells typically remain diploid. The ploidy inference in scPhylogenomics for single-cell RNA sequencing (scRNA-Seq) data aims to detect aneuploidy in individual cells by leveraging genome-wide cellular copy number variations (CNVs) inferred from gene expression patterns. The goal is to identify malignant cells, which can help determine tumor subclones and genetic heterogeneity within tumors, utilizing the SNP-based phylogenomics approach implemented in the subsequent modules of the scPhylogenomics workflow. 


## Analysis
This module infers cell ploidy for annotated scRNA-Seq data produced by the [cell-typing module](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/cell-typing) using the [copykat R package](https://github.com/navinlabcode/copykat). It classifies cells as either `tumor (malignant)`, `normal` or `unknown` status to facilitate tumor cell-type-specific cancer mutational profiling using other [scPhylogenomics analyses modules](https://github.com/ewafula/scPhylogenomics/tree/main/analyses).

#### Module directory structure
Illustration of module directory structure using 10x scRNA-Seq dataset from the scPhylogenomics analysis workflow study. 
```
.
|-- 01-classify-cell-ploidy.R
|-- run-ploidy-inference.sh
|-- README.md
|-- objects
|   |-- AML
|   |   `-- AML-cell-ploidy-prediction.rds
|   `-- TNBC
|       `-- TNBC-cell-ploidy-prediction.rds
|-- plots
|   |-- AML
|   |   |-- AML-cell-ploidy-prediction-umap.pdf
|   |   `-- AML_copykat_heatmap.jpeg
|   `-- TNBC
|       |-- TNBC-cell-ploidy-prediction-umap.pdf
|       `-- TNBC_copykat_heatmap.jpeg
|-- results
    |-- AML
    |   `-- AML-cell-ploidy-prediction.tsv.gz
    `-- TNBC
       `-- TNBC-cell-ploidy-prediction.tsv.gz
```
**Note**: Seurat objects (RDS files) are not uploaded on the GitHub repository because of file size limitations. 

## General usage of scripts

#### `run-ploidy-inference.sh`
This script serves as a module wrapper script to automate project ploidy inference tasks and can be adapted for any project. All specified paths in this script are relative to the module directory. Therefore, it should always be executed relative to this [module directory](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/ploidy-inference). 

###### Example usage:
```
bash run-ploidy-inference.sh
```

#### `01-classify-cell-ploidy.R`
This script utilizes the [copykat R package](https://github.com/navinlabcode/copykat) to infer cell ploidy from expression counts matrix of an annotated Seurat object with either `integrated samples`, `merged samples` or `a single sample` of a project. It then classifies the ploidy cell predictions based on cell malignancy status. The results of the cell ploidy predictions and malignancy classifications are saved back into the object's metadata.

The input annotated RDS Seurat objects are sourced from the [object folder of cell-typing module object](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/cell-typing/objects) (`{PROJECT_NAME}-{ANNOTATION-METHOD}-annotations.rds`). This script will generate an RDS Seurat object, named `{PROJECT_NAME}-cell-ploidy-prediction.rds`, which will contain updated metadata saved in the `objects/<PROJECT_NAME>/` sub-directory, a summary classification tables, named `{PROJECT_NAME}-cell-ploidy-prediction.tsv.gz` saved in the `results/<PROJECT_NAME>/` sub-directory, and corresponding UMAP plots, named `{PROJECT_NAME}-cell-ploidy-prediction-umap.pdf` saved in the `plots/<PROJECT_NAME>/` sub-directory of this module.

###### Example usage:
```
Rscript --vanilla 01-classify-cell-ploidy.R \
  --project TNBC \
  --annot_object TNBC-custom-annotations.rds \
  --nthreads 8
```

###### Argument descriptions:
```
Usage: 01-classify-cell-ploidy.R [options]


Options:
	--project=CHARACTER
		A valid scPhylogenomics project name

	--annot_object=CHARACTER
		Project cell type annotated Seurat object RDS file name
    
    --plot_genes
		Plot copy number heatmap with gene by cell matrix - [default "TRUE"]
              If set to FALSE, copy number heatmap not be generated

	--nthreads=CHARACTER
		Number of threads (CPUs) [default = 4]

	-h, --help
		Show this help message and exit
```