# Phylogeny Inference

## Purpose
This module conducts single-cell clonal phylogenetic analysis by integrating mutation-based phylogenetic inference with advanced clonal population clustering. The goal is to elucidate the evolutionary relationships among cells within tumor samples.

The analysis begins with construction of phylogenetic trees from SNP-based multiple sequence alignments, supporting both [FastTree](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490) and [IQ-TREE](https://ecoevorxiv.org/repository/view/8916/) for robust model selection and tree inference.

To identify distinct clonal populations, the module utilizes a unified clustering interface supporting two powerful approaches:

1. Manifold Clustering: A deep learning approach using a Variational Autoencoder (VAE) via [SNPmanifold](https://pmc.ncbi.nlm.nih.gov/articles/PMC12465888/).

2. Hierarchical Clustering: A hybrid algorithm employing Ward's linkage and weighted non-negative matrix factorization (WNMF) for broad and fine-grained subclonal resolution.

Finally, clonal phylogeny visualizations are generated, overlaying clonal identity onto the topology of the mutation-based tree. This enables a clear interpretation of subclonal evolution and lineage divergence. The analysis process is wrapped in reproducible scripts that provide a cohesive pipeline for reconstructing and visualizing tumor clonal evolutionary dynamics at a single-cell resolution.

## Analysis
This module provides a reproducible pipeline for single-cell clonal phylogenetic analysis by integrating mutation-based tree inference with deep learning or WNMF-based clonal clustering.

#### Module directory structure
Illustration of module directory structure using 10x scRNA-Seq dataset from the scPhylogenomics analysis workflow study.
```
.
|-- 01-phylogeny-inference.py
|-- 02-filter-variants.py
|-- 03-snp-clustering.py
|-- 04-infer-clonal-phylogeny.R
|-- run-phylogeny-inference.sh
|-- README.md
|-- inputs
|   |-- AML-LE1-annotation_categories.tsv.gz
|   |-- AML-LE1-annotation_ploidy.tsv.gz
|   `-- TNBC-TNBC5-annotation_categories.tsv.gz
|-- utils
|   |-- csp_utils.py
|   |-- hierarchical_utils.py
|   `-- manifold_utils.py
|-- plots
|   |-- AML
|   |   |-- LE1.all-cell-types.clusters.hierarchical.png
|   |   |-- LE1.all-cell-types.clusters.manifold.png
|   |   |-- LE1.all-cell-types.fasttree.hierarchical.clonal_tree.png
|   |   |-- LE1.all-cell-types.fasttree.hierarchical.lineage_tree.png
|   |   |-- LE1.all-cell-types.fasttree.hierarchical.ploidy_tree.png
|   |   |-- LE1.all-cell-types.fasttree.manifold.clonal_tree.png
|   |   |-- LE1.all-cell-types.fasttree.manifold.lineage_tree.png
|   |   |-- LE1.all-cell-types.fasttree.manifold.ploidy_tree.png
|   |   |-- LE1.all-cell-types.iqtree.hierarchical.clonal_tree.png
|   |   |-- LE1.all-cell-types.iqtree.hierarchical.lineage_tree.png
|   |   |-- LE1.all-cell-types.iqtree.hierarchical.ploidy_tree.png
|   |   |-- LE1.all-cell-types.iqtree.manifold.clonal_tree.png
|   |   |-- LE1.all-cell-types.iqtree.manifold.lineage_tree.png
|   |   |-- LE1.all-cell-types.iqtree.manifold.ploidy_tree.png
|   |   |-- LE1.all-cell-types.scatter.2d.manifold.png
|   |   `-- LE1.all-cell-types.scatter.3d.manifold.png
|   `-- TNBC
|       |-- TNBC5.select-cell-types.clusters.hierarchical.png
|       |-- TNBC5.select-cell-types.clusters.manifold.png
|       |-- TNBC5.select-cell-types.fasttree.hierarchical.clonal_tree.png
|       |-- TNBC5.select-cell-types.fasttree.manifold.clonal_tree.png
|       |-- TNBC5.select-cell-types.iqtree.hierarchical.clonal_tree.png
|       |-- TNBC5.select-cell-types.iqtree.manifold.clonal_tree.png
|       `-- TNBC5.select-cell-types.scatter.2d.manifold.png
|-- results
    |-- AML
    |   |-- LE1.all-cell-types.cellSNP.base.vcf.gz
    |   |-- LE1.all-cell-types.cellSNP.samples.tsv.gz
    |   |-- LE1.all-cell-types.cellSNP.tag.AD.mtx.gz
    |   |-- LE1.all-cell-types.cellSNP.tag.DP.mtx.gz
    |   |-- LE1.all-cell-types.cellSNP.tag.OTH.mtx.gz
    |   |-- LE1.all-cell-types.clones.hierarchical.tsv.gz
    |   |-- LE1.all-cell-types.clones.manifold.tsv.gz
    |   |-- LE1.all-cell-types.fasttree.hierarchical.clonal.treefile
    |   |-- LE1.all-cell-types.fasttree.hierarchical.lineage.treefile
    |   |-- LE1.all-cell-types.fasttree.hierarchical.ploidy.treefile
    |   |-- LE1.all-cell-types.fasttree.manifold.clonal.treefile
    |   |-- LE1.all-cell-types.fasttree.manifold.lineage.treefile
    |   |-- LE1.all-cell-types.fasttree.manifold.ploidy.treefile
    |   |-- LE1.all-cell-types.fasttree.tree
    |   |-- LE1.all-cell-types.iqtree.hierarchical.clonal.treefile
    |   |-- LE1.all-cell-types.iqtree.hierarchical.lineage.treefile
    |   |-- LE1.all-cell-types.iqtree.hierarchical.ploidy.treefile
    |   |-- LE1.all-cell-types.iqtree.manifold.clonal.treefile
    |   |-- LE1.all-cell-types.iqtree.manifold.lineage.treefile
    |   |-- LE1.all-cell-types.iqtree.manifold.ploidy.treefile
    |   `-- LE1.all-cell-types.iqtree.tree
    `-- TNBC
        |-- cell_barcodes_SNPmanifold_clusters.tsv.gz
        |-- TNBC5.select-cell-types.cellSNP.base.vcf.gz
        |-- TNBC5.select-cell-types.cellSNP.samples.tsv.gz
        |-- TNBC5.select-cell-types.cellSNP.tag.AD.mtx.gz
        |-- TNBC5.select-cell-types.cellSNP.tag.DP.mtx.gz
        |-- TNBC5.select-cell-types.cellSNP.tag.OTH.mtx.gz
        |-- TNBC5.select-cell-types.clones.hierarchical.tsv.gz
        |-- TNBC5.select-cell-types.clones.manifold.tsv.gz
        |-- TNBC5.select-cell-types.fasttree.hierarchical.clonal.treefile
        |-- TNBC5.select-cell-types.fasttree.manifold.clonal.treefile
        |-- TNBC5.select-cell-types.fasttree.tree
        |-- TNBC5.select-cell-types.iqtree.hierarchical.clonal.treefile
        |-- TNBC5.select-cell-types.iqtree.manifold.clonal.treefile
        `-- TNBC5.select-cell-types.iqtree.tree
``` 

## General usage of scripts

#### `run-phylogeny-inference.sh`
This script serves as a module wrapper script to automate project ploidy inference tasks and can be adapted for any project. All specified paths in this script are relative to the module directory. Therefore, it should always be executed relative to this [module directory](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/phylogeny-inference). 

###### Example usage:
```
bash run-snv-calling.sh
```

#### `01-phylogeny-inference.py`
This script infers a phylogenetic tree from SNP-based DNA multiple sequence alignment (MSA) files estimated by the upstream [snv-calling module](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/snv-calling). Designed for scalability, it automates the process of phylogenetic tree inference supporting two widely used programs:

* [FastTree](https://morgannprice.github.io/fasttree/): Constructs trees using the GTR substitution model and accounts for site rate variation with a gamma distribution.

* [IQ-TREE](https://iqtree.github.io/): Leverages advanced model selection (JC, K80, HKY, GTR) with various rate heterogeneity options (MF+G+I), automatically determining the best-fitting model.
 

###### Example usage:
```
python3 01-phylogeny-inference.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types \
  --method iqtree \
  --threads 16
```

###### Argument descriptions:
```
usage: 01-phylogeny-inference.py [-h] --project PROJECT --sample SAMPLE --cell_type CELL_TYPE
                                 [-t THREADS] [--method {iqtree,fasttree}] [--log LOG]

Phylogeny Inference Module (scPhylogenomics)

options:
  -h, --help            show this help message and exit
  --project PROJECT     Project ID (e.g., AML)
  --sample SAMPLE       Sample ID (e.g., LE1)
  --cell_type CELL_TYPE
                        Cell type string (e.g., all-cell-types)
  -t THREADS, --threads THREADS
                        Number of parallel threads (default = 4).
  --method {iqtree,fasttree}
                        Tree inference method (default = iqtree).
  --log LOG             Log file name.         Keep intermediate trimmed alignments.
```

#### `02-filter-variants.py`
This script filters the pre-computed `cellSNP` results (`variants` and `cells`) generated by the upstream [snv-calling module](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/snv-calling). It uses the retained `SNPs` and `cell barcodes` from the trimmed MSA as inclusion lists. This preprocessing step ensures that the matrices used for `manifold clonal clustering` contain only high-quality, relevant variants and cells that correspond to the MSA and the binary profile matrix utilized in `phylogeny inference` and `hierarchical clonal clustering`.

The script loads cellSNP outputs (`VCFs` and `matrices`) and subsets them based on:

* Barcode Whitelist: Defined in {prefix}cellSNP.barcodes.filtered.tsv.gz.
* Site Whitelist: Defined in {prefix}cellSNP.sites.filtered.tsv.gz.

The filtered matrices (AD, DP) and VCFs are saved to the project results directory, for subsequent clonal clustering.

###### Example usage:
```
python3 02-filter-variants.py \
  --project TNBC \
  --sample TNBC5 \
  --cell_type select-cell-types
```

###### Argument descriptions:
```
usage: 02-filter-variants.py [-h] --project PROJECT --sample SAMPLE --cell_type CELL_TYPE

Filter pre-computed cellsnp-lite results using inclusion lists.

options:
  -h, --help            show this help message and exit
  --project PROJECT     A valid scPhylogenomics project name string (e.g., TNBC).
  --sample SAMPLE       Sample name string (e.g., TNBC5).
  --cell_type CELL_TYPE
                        Cell types being analyzed (e.g., select-cell-types).
```

#### `03-snp-clustering.py`
This script serves as a unified entry point for SNP-based clonal clustering. It supports two distinct methodologies: `Manifold` (Deep Learning) and `Hierarchical` (NMF). It delegates the core logic to specialized utility modules (`utils/manifold_utils.py` and `utils/hierarchical_utils.py`).

**Method 1: Manifold Clustering** (`--method manifold`)
Uses a `Variational Autoencoder` (VAE) implemented in PyTorch ([SNPmanifold](https://github.com/StatBiomed/SNPmanifold/tree/main)).
* **Environment:** Automatically dispatches execution to a specific conda environment (snpmanifold_env) to handle deep learning dependencies.
* **Process:** Performs training, dimensionality reduction (UMAP), and clustering (K-Means/GMM) on the latent space.
* **Outputs:** 3D scatter plots and clonal assignment tables.

**Method 2: Hierarchical Clustering** (`--method hierarchical`)
Uses a hybrid two-stage approach:
* **Broad Clade Clustering:** Hierarchical clustering (Ward's method) identifies major evolutionary clades (`--broad_k`).
* **Narrow Subclone Clustering:** Weighted Non-negative Matrix Factorization (WNMF) decomposes profiles within clades to find subtle subclones (`--min_k` to `--max_cluster`).

We recommend using a two-step clustering approach for either method. First, run the analysis without utilizing the `--final` option and with a large value of `k` to help determine the optimal number of clusters. After identifying the optimal clusters, rerun the analysis with the `--final` option set. This will determine whether the clustering results and plots are saved in the workflow [scratch directory](https://github.com/ewafula/scPhylogenomics/tree/main/scratch) for evaluation or the module [results directory](https://github.com/ewafula/scPhylogenomics/tree/main/analyses/phylogeny-inference/results).

###### Example usage (Manifold):
```
python3 03-snp-clustering.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types \
  --method manifold \
  --max_cluster 3 \
  --final 
```

###### Example usage (Hierarchical):
```
python3 03-snp-clustering.py \
  --project TNBC \
  --sample TNBC5 \
  --cell_type select-cell-types \
  --method hierarchical \
  --broad_k 3 \
  --final
```

###### Argument descriptions:
```
usage: 03-snp-clustering.py [-h] --project PROJECT --sample SAMPLE --cell_type CELL_TYPE
                            [--method {manifold,hierarchical}] [--max_cluster MAX_CLUSTER]
                            [--final] [--force_k FORCE_K] [--broad_k BROAD_K] [--min_k MIN_K]

Unified SNP Clustering for scPhylogenomics (Manifold or Hierarchical)

options:
  -h, --help            show this help message and exit
  --project PROJECT     Project ID (e.g., AML)
  --sample SAMPLE       Sample ID (e.g., LE1)
  --cell_type CELL_TYPE
                        Cell type string (e.g., all-cell-types)
  --method {manifold,hierarchical}
                        Clustering method to use.
  --max_cluster MAX_CLUSTER
                        Maximum number of clusters (Manifold) or Max k for narrow search (Hierarchical).
  --final               Flag to indicate final run; copies outputs to module results/ and plots/.
  --force_k FORCE_K     [Hierarchical] Force specific narrow clusters per clade.
  --broad_k BROAD_K     [Hierarchical] Number of broad hierarchical clades (default = 4).
  --min_k MIN_K         [Hierarchical] Min k for narrow cluster search (default = 2).
```

### `04-infer-clonal-phylogeny.R`
This script integrates the phylogenetic tree (from 01) with the clonal clusters (from 03). It visualizes the evolutionary relationships by overlaying clonal identity onto the tree topology.

The workflow involves:
1. **Rooting:** Roots the tree via a reference sequence or midpoint rooting.
2. **Pruning:** Removes outlier tips with excessively long branches (`--percentile_outlier`).
3. **Annotation:** Maps clonal IDs (`Manifold` or `Hierarchical`) or sample IDs to the tree tips.

4, **Visualization:** Generate rectangular-layout trees with tips colored by clone/sample and clades labeled accordingly.

The script automatically detects which clustering outputs (`Manifold` or `Hierarchical`) are available in the results directory and generates plots for each.


###### Example usage:
```
Rscript 04-infer-clonal-phylogeny.R \
  --project AML \
  --refseq \
  --annot_type Clone \
  --rescale
```

###### Argument descriptions:
```
Usage: 04-infer-clonal-phylogeny.R [options]

Options:
  --project=CHARACTER
    A valid scPhylogenomics project name

  --refseq
    Whether to include reference sequence [default: FALSE]

  --rescale
    Whether to rescale branch lengths to compress extremes [default: FALSE]

  --percentile_outlier=PERCENTILE
    Determine likely outlier super long branches to exclude [default: 0.99]

  --annot_type=CHARACTER
    Clone metadata annotation - [default: Clone]
    - Available choices: Clone, Lineage, and Sample
    - Sample choice only works with merged samples alignment
    
  --cell_category=CHARACTER
    Optional cell annotation categories file (if not provided, skipped)

  -h, --help
    Show this help message and exit
```

