# Cell Typing

## Purpose

This module consolidates preprocessed 10x single-cell RNA-Seq (scRNA-Seq) data produced by the [data-processing module](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/data-preprocessing) by either merging or integrating samples and annotating cell types to facilitate cancer mutational profiling.

## Analysis
This module consolidates filtered 10x scRNA-Seq samples produced by the [data-processing module](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/data-preprocessing) for a project into either a `merged` or an `integrated` dataset following the [Seurat](https://www.nature.com/articles/s41587-023-01767-y) scRNA-Seq analysis protocol. The resulting dataset is then annotated with cell types using [SingleR](https://www.nature.com/articles/s41590-018-0276-y) against either publicly available scRNA-Seq reference databases or user-provided custom references to facilitate cell-type-specific cancer mutational profiling using other [scPhylogenomics analyses modules](https://github.com/ewafula/scPhylogenomics/tree/main/analyses).

#### Module directory structure
Illustration of module directory structure using 10x scRNA-Seq dataset from the scPhylogenomics analysis workflow study. 
```  
.
|-- 01-consolidate-samples.R
|-- 02-annotate-cell-types.R
|-- run-cell-typing.sh
|-- README.md
|-- inputs
|   |-- AML-LE1-cell-annotation.tsv.gz
|   |-- BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds
|   |-- BoneMarrowMap_SymphonyReference.rds
|   |-- BoneMarrowMap_uwot_model.uwot
|   `-- HBCA-snRNA-seq-all-cells.rds
|-- objects
|   |-- AML
|   |   |-- AML-clusters.rds
|   |   `-- AML-mapping-annotations.rds
|   `-- TNBC
|       |-- TNBC-clusters.rds
|       `-- TNBC-custom-annotations.rds
|-- plots
|   |-- AML
|   |   |-- AML-clusters-umap.pdf
|   |   |-- AML-elbowplot.pdf
|   |   `-- AML-mapping-annotations-umap.pdf
|   `-- TNBC
|       |-- TNBC-clusters-umap.pdf
|       |-- TNBC-custom-annotations-umap.pdf
|       |-- TNBC-custom-delta-heatmap.pdf
|       `-- TNBC-elbowplot.pdf
|-- results
    |-- AML
    |   `-- AML-mapping-annotations.tsv.gz
    `-- TNBC
        `-- TNBC-custom-annotations.tsv.gz      
```
**Note**: Seurat objects (RDS files) are not uploaded on the GitHub repository because of file size limitations. 

## General usage of scripts

#### `run-cell-typing.sh`
This script serves as a module wrapper script to automate project cell type annotation tasks and can be adapted for any project. All specified paths in this script are relative to the module directory. Therefore, it should always be executed relative to this [module directory](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/cell-typing). 

###### Example usage:
```
bash run-cell-typing.sh
```

#### `01-consolidate-samples.R`
This script utilizes the [Seurat](https://satijalab.org/seurat/articles/seurat5_integration) to `merge`, `normalize`, optionally `integrate`, and `cluster` filtered scRNA-Seq samples of a project into a consolidated dataset. Only `normalization` and `clustering` are performed when there is a single sample to be analyzed.

A project's filtered sample data are in the `scPhylogenomics/data/projects/<PROJECT_NAME>/cellranger/<SAMPLE_NAME>/miqc/` subdirectory. The script will create an RDS Seurat object (`{PROJECT_NAME}-clusters.rds`) with consolidated project samples in the `objects/<PROJECT_NAME>/` sub-directory of this module.

###### Example usage:
```
Rscript --vanilla 01-consolidate-samples.R 
  --project TNBC \
  --integrate TRUE \
  --normalize_method LogNormalize \
  --integrate_method Harmony \
  --components 20 \
  --resolution 0.4
```

###### Argument descriptions:
```
Usage: 01-consolidate-samples.R [options]


Options:
    --project=CHARACTER
        Project directory name with 10x Cell Ranger counts sample runs

    --integrate
        Integate project samples - [default "TRUE"]
                      If set to FALSE, project samples will not be integrated

    --normalize_method=CHARACTER
        Normalize the count data - [default "LogNormalize"]
                      - Available choices - LogNormalize, and SCT

    --integrate_method=CHARACTER
        Method for integrating project samples - [default "Harmony"]
                      - Available choices - Harmony, CCA, FastMNN, and RPCA

    --components=NUMBER
        Number of components for PCA [default 30]

    --resolution=NUMBER
        Seurat clustering resolution parameter - [default 0.8]

    -h, --help
        Show this help message and exit
```

#### `02-annotate-cell-types.R`
This script utilizes the [SingleR](https://bioconductor.org/packages/3.21/bioc/vignettes/SingleR/inst/doc/SingleR.html) Bioconductor R package to annotate cell types in the consolidated project dataset (Seurat object). We illustrate the utility of this module using the [scPhylogenomics workflow scRNA-Seq datasets](https://github.com/ewafula/scPhylogenomics/tree/main?tab=readme-ov-file#data-description). The `CopyKAT triple-negative breast cancer (TNBC)` dataset is annotated against the `Human Primary Cell Atlas (HPCA)` reference dataset provided in the [celldex](https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html) Bioconductor R package and the [Human Breast Cell Atlas (HBCA)](https://explore.data.humancellatlas.org/projects/1ffa2223-28a6-4133-a5a4-badd00faf4bc) reference dataset to illustrate the module's capability to annotate cells against built-in scRNA-Seq Bioconductor-celldex references and custom references respectively. In addition, the `MSK SPECTRUM High-grade serous ovarian cancer (HGSOC)` dataset is assigned actual cell types from the original study to illustrate the module's capability to incorporate annotations inferred outside of scPhylogenomics workflow. 

The module also produces a diagnostic delta plot that enable users to inspect the confidence of the predicted cell types across the dataset. As described in the [SingleR vignette](https://bioconductor.org/packages/3.21/bioc/vignettes/SingleR/inst/doc/SingleR.html), each cell, represented as a column in the heatmap, each cell, represented as a column in the heatmap, should ideally have one score that is significantly higher than the others.  This indicates a clear and unambiguous assignment to a cell type. Conversely, if a cell displays a range of similar scores, it suggests that the assignment is uncertain. While this uncertainty is generally not desirable, it may be acceptable if the scores are distributed among similar cell types that are difficult to distinguish.

The script will create an RDS Seurat object (`{PROJECT_NAME}-{ANNOTATION-METHOD}-annotations.rds`) with consolidated and annotated project samples in the `objects/<PROJECT_NAME>/` sub-directory of this module. In addition, to generate summary plots in the `plots/<PROJECT_NAME>/` module sub-directory and cell type annotation tables in the `results/<PROJECT_NAME>/` module sub-directory.

###### Example usage:
1) Using the celldex Bioconductor R package reference dataset - HPCA):
```
Rscript --vanilla  02-annotate-cell-types.R \
  --project TNBC \
  --assay RNA \
  --annot_method celldex \
  --ref_data hpca
```

2) Using a custom reference dataset (Seurat object) - `HBCA`:
```
Rscript --vanilla  02-annotate-cell-types.R \
  --project TNBC \
  --assay RNA \
  --annot_method custom \
  --ref_data inputs/HBCA-snRNA-seq-all-cells.rds
```

3) Using externally inferred cell types (mapping file) - `Leukemia Cell Atlas`
```
# Create cell type annotations mapping file from the Leukemia Cell Atlas
pushd ../../scripts/ > /dev/nul
Rscript -e "rmarkdown::render('create-leukemia-reference-mapping.Rmd', clean = TRUE)"
popd > /dev/null
  
# Assign cell types using mappings from Leukemia Cell Atlas...\n'
Rscript --vanilla  02-annotate-cell-types.R \
  --project AML \
  --assay RNA \
  --annot_method mapping \
  --ref_data inputs/AML-LE1-cell-annotation.tsv.gz
```

###### Argument descriptions:
```
Usage: 02-annotate-cell-types.R [options]


Options:
    --project=CHARACTER
        Project for samples in the Seurat  object

    --assay=CHARACTER
        Project Seurat object assay with normalized counts - [default "RNA"]
                      - Available choices - RNA and SCT

    --annot_method=CHARACTER
        Automatic cell type annotation method - [default "celldex"]
                      - Available choices - celldex, custom, and mapping

    --ref_data=CHARACTER
        Reference dataset - [default "hpca"]
                      - Available celldex reference choices - hpca, monaco_immune, 
                        immgen, dice, novershtern_hematopoietic, blueprint_encode
                      - If 'custom' is set for the '--annot_method' option,
                        provide custom reference Seurat object with normalized 
                        counts similar the project Seurat object (query object)
                        and cell labels column named 'cell_type' 
                      - If 'mapping' is set for '--annot_method' option provide a
                        TSV mapping file of annotations projects samples with the 
                        matching cells (barcodes) as 'cell_id' column and annotated
                        cell labels as 'cell_type' column"

    -h, --help
        Show this help message and exit
```

Suppose a reference dataset obtained from other scRNA-Seq reference databases such as [Azimuth](https://azimuth.hubmapconsortium.org/), [PangaloDB](https://panglaodb.se/index.html), [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home), [Single Cell Protal](https://singlecell.broadinstitute.org/single_cell), and [CZCellxGeneDiscover](https://cellxgene.cziscience.com/) is not available as a normalized Seurat object. In that case, We provide a [script](https://github.com/ewafula/scPhylogenomics/blob/main/scripts/create-reference-dataset.R) for creating a lightweight custom reference Seurat object dataset using a raw counts gene expression matrix and corresponding cell type annotations.

###### Example usage:
```
Rscript --vanilla create-reference-dataset.R \
  --expr_matrix expr_matrix.rds \
  --cell_types cell_types.tsv.gz
  ```
  
 ###### Argument descriptions:
 ```
 Usage: create-reference-dataset.R [options]


Options:
    --expr_matrix=CHARACTER
        Cell types gene expression raw counts matrix RDS file 
                      with genes as rows and cells (barcodes) as columns
                     - convert to a sparse matrix if too large

    --cell_types=CHARACTER
        A cell type annotation TSV file with the counts matrix
                      cells (barcodes) as 'cell_id' column and annotated cell 
                      labels as 'cell_type' column

    --normalize_method=CHARACTER
        Normalize the count data - [default "LogNormalize"]
                      - Available choices - LogNormalize, and SCT

    --rds_file=CHARACTER
        RDS output file name for annotation object - [default "reference-data.rds"]

    -h, --help
        Show this help message and exit
 
```
Here, we provide an example of how to convert a reference data object from [Azimuth single-cell reference database](https://azimuth.hubmapconsortium.org/) for use in annotating cell types in the scPhylogenmics workflow.

This is a [Human PBMC reference dataset](https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC), which consists of 24 samples processed with a CITE-seq panel of 228 TotalSeq A antibodies to generate single-cell RNA and ADT data.

Adopted from [Seurat Azimuth annotation](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html)
```
# load libraries
library(Seurat)
library(SeuratData)

# check for available reference data sets on Azimuth
> available_data <- SeuratData::AvailableData()
> available_data[grep("Azimuth", available_data[, 3]), 1:3]
adiposeref.SeuratData         adiposeref   1.0.0     Azimuth Reference: adipose
bonemarrowref.SeuratData   bonemarrowref   1.0.0  Azimuth Reference: bonemarrow
fetusref.SeuratData             fetusref   1.0.0       Azimuth Reference: fetus
heartref.SeuratData             heartref   1.0.0       Azimuth Reference: heart
humancortexref.SeuratData humancortexref   1.0.0 Azimuth Reference: humancortex
kidneyref.SeuratData           kidneyref   1.0.2      Azimuth Reference: kidney
lungref.SeuratData               lungref   2.0.0        Azimuth Reference: lung
mousecortexref.SeuratData mousecortexref   1.0.0 Azimuth Reference: mousecortex
pancreasref.SeuratData       pancreasref   1.0.0    Azimuth Reference: pancreas
pbmcref.SeuratData               pbmcref   1.0.0        Azimuth Reference: pbmc
tonsilref.SeuratData           tonsilref   2.0.0      Azimuth Reference: tonsil

# install pbmc reference dataset
> InstallData("pbmcref")

# get pbmc reference seurat object
> pbmcref <- LoadData("pbmcref", type="azimuth")
> pbmcref
$map
An object of class Seurat 
5228 features across 36433 samples within 2 assays 
Active assay: refAssay (5000 features, 0 variable features)
 1 layer present: data
 1 other assay present: ADT
 2 dimensional reductions calculated: refUMAP, refDR

$plot
An object of class Seurat 
1 features across 24760 samples within 1 assay 
Active assay: RNA (1 features, 0 variable features)
 2 layers present: counts, data
 1 dimensional reduction calculated: refUMAP

# reset defaut assay ()
DefaultAssay(pbmcref$map) <- "ADT"

# rename ADT to RNA - sPhylogenomics requires either RNA or SCT assay with data layer
> pbmcref <- SeuratObject::RenameAssays(object = pbmcref$map,  ADT = 'RNA')
>  pbmcref
An object of class Seurat 
5228 features across 36433 samples within 2 assays 
Active assay: ADT (228 features, 224 variable features)
 1 layer present: data
 1 other assay present: refAssay
 2 dimensional reductions calculated: refUMAP, refDR

# create cell_type metadata column required by scPhylogenomics module
pbmcref$cell_type <- pbmcref@meta.data$celltype.l2
> pbmcref@meta.data |> head()
                    celltype.l1    celltype.l2           celltype.l3 ori.index
L1_AAACGAATCCTCACCA     other T            gdT                 gdT_3        27
L1_AAACGCTAGAGCATTA       CD8 T        CD8 TEM             CD8 TEM_2        30
L1_AAACGCTCAACGATCT       CD8 T        CD8 TCM             CD8 TCM_1        35
L1_AAACGCTGTGCTCGTG     other T            dnT                 dnT_2        40
L1_AAACGCTTCTTGGTCC           B B intermediate B intermediate lambda        42
L1_AAAGAACCAAGCGGAT       CD4 T        CD4 TCM             CD4 TCM_3        46
                    nCount_refAssay nFeature_refAssay      cell_type nCount_RNA
L1_AAACGAATCCTCACCA               0                 0            gdT          0
L1_AAACGCTAGAGCATTA               0                 0        CD8 TEM          0
L1_AAACGCTCAACGATCT               0                 0        CD8 TCM          0
L1_AAACGCTGTGCTCGTG               0                 0            dnT          0
L1_AAACGCTTCTTGGTCC               0                 0 B intermediate          0
L1_AAAGAACCAAGCGGAT               0                 0        CD4 TCM          0
                    nFeature_RNA
L1_AAACGAATCCTCACCA            0
L1_AAACGCTAGAGCATTA            0
L1_AAACGCTCAACGATCT            0
L1_AAACGCTGTGCTCGTG            0
L1_AAACGCTTCTTGGTCC            0
L1_AAAGAACCAAGCGGAT            0

# save the object to RDS file to use as reference annotation object in scPylogenomics
saveRDS(pbmcref, "pbmcref.rds")

# use the converted reference object to annotate cell types
Rscript --vanilla  02-annotate-cell-types.R \
  --project BLOOD_PRJ \
  --assay RNA \
  --annot_method custom \
  --ref_data pbmcref.rds
```

We provide another example of how to convert a reference `Anndata h5ad` data object from [CZCellxGeneDiscover](https://cellxgene.cziscience.com/) for use in annotating cell types in the phylogenomics workflow.

This is a [Human blood and bone marrow reference dataset](https://cellxgene.cziscience.com/collections/93eebe82-d8c3-41bc-a906-63b5b5f24a9d), which consists of surface markers to cellular identities and biological processes across all main hematopoietic cell types in from two samples of a healthy young donor processed with BD Rhapsody Targeted mRNA 

Adopted from [Single Cell Genomics, Transcriptomics & Proteomics Tutorials](https://www.youtube.com/@Collection_of_online_tutorials)
```
# download Anndata object in h5ad format
wget https://datasets.cellxgene.cziscience.com/7dd9be25-93ed-452e-acb1-addda8565ca2.h5ad

# mkdir Blood_BM_normal/

# extract raw data from the h5ad object using python
>>> import scanpy as sc
>>> import numpy as np
>>> import pandas as pd
>>> from scipy import io
>>> from scipy.sparse import coo_matrix, csr_matrix
>>> import os

>>> adata = sc.read_h5ad("7dd9be25-93ed-452e-acb1-addda8565ca2.h5ad")

>>> adata = adata.raw.to_adata()

>>> io.mmwrite("Blood_BM_normal/matrix", adata.X.T)

>>> with open("Blood_BM_normal/barcodes.tsv", "w") as f:
      for item in adata.obs_names:
        f.write(item + "\n")

>>> with open("Blood_BM_normal/features.tsv", "w") as f:
      for item in adata.var_names:
        f.write(item + "\n")

>>> with open("Blood_BM_normal/features.tsv", "w") as f:
      for item in adata.var_names:
        f.write(item + "\n")

# extract metadata
>>> adata.obs.to_csv("Blood_BM_normal/metadata.csv")

# cell type annotation column in metadata named as required in scPhylogenomics (`cell_type`)
>>> adata
AnnData object with n_obs × n_vars = 15502 × 460
    obs: 'Cluster_ID', 'Sample', 'Cell_label', 'is_primary_data', 'organism_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'assay_ontology_term_id', 'tissue_ontology_term_id', 'Genotype', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'disease_ontology_term_id', 'donor_id', 'cell_type_ontology_term_id', 'suspension_type', 'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
    var: 'name', 'feature_is_filtered', 'Unnamed: 0', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
    uns: 'cell_type_ontology_term_id_colors', 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'title'
    obsm: 'X_projected', 'X_projectedmean'

#  switch to R and load libraries
> library(Seurat)
> library(tidyverse)
> library(R.utils)

# compress the raw data
> gzip("Blood_BM_normal/barcodes.tsv")
> gzip("Blood_BM_normal/features.tsv")
> gzip("Blood_BM_normal/matrix.mtx")

# load raw data in Seurat
> Raw_data <- Seurat::Read10X(data.dir = "Blood_BM_normal/", gene.colum=1)
> metadata <- read.csv("Blood_BM_normal/metadata.csv")
> rownames(metadata) <- metadata[["X"]]
> metadata[["X"]] <- NULL 

# create Seurat object
> Seurat_obj <- SeuratObject::CreateSeuratObject(count = Raw_data, meta.data = metadata)
> Seurat_obj
An object of class Seurat 
460 features across 15502 samples within 1 assay 
Active assay: RNA (460 features, 0 variable features)
 1 layer present: counts

# Normalize data
> Seurat_obj <- Seurat::NormalizeData(Seurat_obj)
> Seurat_obj
An object of class Seurat 
460 features across 15502 samples within 1 assay 
Active assay: RNA (460 features, 0 variable features)
 2 layers present: counts, data

# save the object to RDS file to use as reference annotation object in scPylogenomics
> saveRDS(Seurat_obj, "blood_and_bm.rds")

# use the converted reference object to annotate cell types
Rscript --vanilla  02-annotate-cell-types.R \
  --project BLOOD_PRJ \
  --assay RNA \
  --annot_method custom \
  --ref_data blood_and_bm.rds
```

In addition, we provide an example [RMarkdown notebook](https://github.com/ewafula/scPhylogenomics/tree/main/scripts/create-leukemia-reference-mapping.nb.html) that demonstrates how to annotate leukemic scRNA-seq datasets using the [Single-cell Transcriptional Atlas of Human Hematopoiesis](https://aacrjournals.org/bloodcancerdiscov/article/6/4/307/763153/Single-cell-Transcriptional-Atlas-of-Human) as a reference within the [cell-typing module](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/cell-typing). This atlas encompasses **263,159 human bone marrow hematopoietic cells**, with balanced representation of CD34⁺ stem/progenitor and differentiated populations. HSPC annotations were rigorously validated against decades of functional studies, with the transcriptional HSC state showing strong concordance with bona fide LT-HSCs. Using the provided notebook, users can project leukemic scRNA-seq datasets (e.g., AML) for **cell type prediction directly from raw count matrices**, enabling accurate contextualization of hematopoietic states.