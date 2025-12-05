# scPhylogenomics

## Overview
scPhylogenomics is a fully modular, scalable, and reproducible computational workflow for `single-cell RNA sequencing (scRNA-Seq)-based phylogenomics`. The workflow addresses key analytical stages, from high-throughput data preprocessing (ambient RNA removal, doublet/multiplet exclusion, and mitochondrial quality control) through sample integration and reference-assisted cell-type annotation, to ploidy inference that stratifies malignant from non-malignant cells. Leveraging a streamlined somatic SNV calling approach restricted to tumor cell barcodes, the workflow summarizes variant profiles. It automates the construction of mutation matrices and multiple sequence alignments suitable for downstream evolutionary inference. Finally, robust phylogenetic and clonal population structure reconstruction is achieved by integrating tree-building (FastTree/IQ-TREE) with weighted non-negative matrix factorization (WNMF)-based clonal clustering, producing visualizations and metrics for tumor heterogeneity and clonal evolution. The presented pipeline emphasizes automation, reproducibility, multi-project compatibility, and direct adaptability to high-performance and cloud computing environments, thus addressing key needs in modern single-cell cancer phylogenomics and computational oncology.

## Data description 

### Study data 

##### CopyKAT
Single-cell transcriptomic data for human tumors from the [CopyKAT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148673) (Copy number Karyotyping of Aneuploid Tumors) methodology development study to distinguish normal cell types in the tumor microenvironment from malignant cells and to resolve clonal substructure within the tumor utilizing cell copy number profiles. 

- TNBC5 sample from an individual with `triple-negative breast cancer (TNBC)`

##### Acute Myeloid Leukemia (AML) dataset
**Pending:**
 Need brief description of the AML dataset for sample `Le1` 

 - Le1 sample from an individual with `acute myeloid leukemia (AML)`
 
### Processing raw sequence data
We provide a [script](scripts/run-cellranger-count.sh) to process scRNA-Seq data using the [Cell Ranger Count pipeline](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct), developed by 10x Genomics to transforms raw sequencing reads (FASTQ files) into a structured gene expression matrix suitable for downstream analysis. The pipeline aligns sample sequencing reads in FASTQ files to a reference transcriptome, and performs filtering, barcode counting, and UMI counting to generate raw and filtered feature-barcode matrix and position-sorted read alignments. The resulting sample output data is stored in the scPhylogenomics projects data directory (see directory tree structure below).

The single cell reference transcriptomes datasets releases are avilable to download from the [10x Genomics Cell Ranger Download center](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads) and moved to the `data/refdata directory`.

##### Example usage:
```
bash run-cellranger-count.sh
Usage: run-cellranger-count.sh <PROJECT_NAME> <SAMPLE SHEET> <REFERENCE DATA DIR>

cd sPhylogenomics/scripts/
bash  run-cellranger-count.sh CopyKAT sample_sheet.tsv ../data/refdata
```

##### Example sample sheet (tab-separated):
|sample_name|path_to_sample_fastq_dir    |
|:----------|:---------------------------|
|Le1      | /path/to/sample/fastq/dir/ |
|TNBC5      | /path/to/sample/fastq/dir/ |

##### Example ouput: 
```
.
├── data
    ├── projects
    |   └── CopyKAT (project name)
    |   |   ├── TNBC5 (sample name)
    |   |        ├── outs
    |   |            ├── analysis
    |   |            ├── cloupe.cloupe
    |   |            ├── filtered_feature_bc_matrix
    |   |            ├── filtered_feature_bc_matrix.h5
    |   |            ├── metrics_summary.csv
    |   |            ├── molecule_info.h5
    |   |            ├── possorted_genome_bam.bam
    |   |            ├── possorted_genome_bam.bam.bai
    |   |            ├── raw_feature_bc_matrix
    |   |            ├── raw_feature_bc_matrix.h5
    |   |            └── web_summary.html`  
    |   |
    |   └── AML (project name)
    |       ├── Le1 (sample name)
    |           ├── outs
    |               ├── analysis
    |               ├── cloupe.cloupe
    |               ├── filtered_feature_bc_matrix
    |               ├── filtered_feature_bc_matrix.h5
    |               ├── metrics_summary.csv
    |               ├── molecule_info.h5
    |               ├── possorted_genome_bam.bam
    |               ├── possorted_genome_bam.bam.bai
    |               ├── raw_feature_bc_matrix
    |               ├── raw_feature_bc_matrix.h5
    |               └── web_summary.html`       
    |               
    └── refdata
        ├── cellrange_reference.info
        ├── fasta
        |   ├── genome.fa
        |   └── genome.fa.fai
        ├── genes
        |   └── genes.gtf.gz
        ├── refdata-gex-GRCh38-2024-A
        ├── reference.json
        └── star
            ├── Genome
            ├── SA
            ├── SAindex
            ├── chrLength.txt
            ├── chrName.txt
            ├──chrNameLength.txt
            ├── chrStart.txt
            ├── exonGeTrInfo.tab
            ├── exonInfo.tab
            ├── geneInfo.tab
            ├── genomeParameters.txt
            ├── sjdbInfo.txt
            ├── sjdbList.fromGTF.out.tab
            ├── sjdbList.out.tab
            └──  transcriptInfo.tab             
```


## Software dependencies
We recommend performing analyses using the [scPhylogenomics Docker image](https://hub.docker.com/repository/docker/ewafula/scphylogenomics/general) due to the large number of software dependencies required to work through the complete workflow, including all analysis modules. While we encourage development within the Docker container, it is also possible to conduct analyses without it, if preferred. If you choose to analyze without Docker, be sure to install all the necessary software dependencies listed in the [Dockerfile](https://github.com/ewafula/scPhylogenomics/blob/main/Dockerfile) and the module scripts in your local environment. 


### Analysis in the Docker container
The scPhlogenomics Docker image is based on a a versioned [`tidyverse`](https://hub.docker.com/r/rocker/tidyverse) image from the [Rocker Project](https://www.rocker-project.org/) (`R version 4.4.0`) which allows development via RStudio in the Docker container (instance of the image) if users to a customize the scPhlogenomics module code. Otherwise, module analyses are conducted locally on the terminal in a Docker container.

##### 1). Clone the scPhylogenomics repository from GitHub
Clone repository and change into scPhylogenomics local repository directory.
```bash
git clone https://github.com/ewafula/scPhylogenomics.git
cd scPhylogenomics/
```

##### 2). Download the latest scPhylogenomics Docker image from DockerHub
Download the most latest scPhylogenomics Docker image from DockerHub and check if the Docker image was sucessfull downloaded
```bash
docker pull ewafula/scphylogenomics:latest
docker image ls 

REPOSITORY                                       TAG       IMAGE ID       CREATED        SIZE
ewafula/scphylogenomics                          latest    c1b310fdfbbe   11 hours ago   10GB
```

##### 3). Create a Docker container
Set the volume to point to the local scPhylogenomics repository directory, with a container name of your choice (e.g., `scphylogenomics`) and a local RStudio password of your choice (e.g., `pass`) and check if the `scphylogenomics` continer was successfully created.
```bash
docker run --name scphylogenomics -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/scPhylogenomics ewafula/scphylogenomics:latest
docker container ls 

CONTAINER ID   IMAGE                            COMMAND   CREATED              STATUS              PORTS                    NAMES
a94205e27ae8   ewafula/scphylogenomics:latest   "/init"   About a minute ago   Up About a minute   0.0.0.0:8787->8787/tcp   scphylogenomics
```

##### 4). Run scPhlogenomics module analyses
The scPhlogenomics module analyses can peformed on terminal in docker container running locally by changing into scPhylogenomics mounted volume
```bash
docker exec -ti scphylogenomics bash
cd /home/rstudio/scPhylogenomics/
```

Module customization in Rstudio can performed by navigating to `localhost:8787` in the Linux or Mac OS X user browser. The username for login will be `rstudio`, and the password will be the same as the password set with the `docker run` command above (i.e., `pass`).

### Analysis in the Singuarity container

While Docker is recommended for local development and environments where root access is available, **Singularity is often the better choice in high-performance computing (HPC) or shared cluster environments where Docker is not supported or root privileges are restricted**. The Singularity image is directly built from the DockerHub image, ensuring the same software stack and reproducibility across platforms.

##### 1). Clone the scPhylogenomics repository from GitHub
Clone repository and change into the local `scPhylogenomics` directory:
```bash
git clone https://github.com/ewafula/scPhylogenomics.git
cd scPhylogenomics/
````

##### 2). Pull the latest scPhylogenomics Singularity image

Download the Singularity image (`.sif`) file from DockerHub. It is recommended to store this in the `images/` subdirectory of the repository:

```bash
mkdir -p images
cd images/
singularity pull images/scphylogenomics.sif docker://ewafula/scphylogenomics:latest
cd ../
```

##### 3). Run scPhylogenomics module analyses

Analyses are run using `singularity exec` while binding your local `scPhylogenomics` repository into the container. This ensures that all inputs, scripts, and outputs are accessible inside and outside the container.

Here are a few examples of module analyses:

###### data-preprocessing module

```bash
cd analyses/data-preprocessing/

singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/data-preprocessing \
  ../../images/scphylogenomics.sif \
  bash run-data-preprocessing.sh
```

###### cell-typing module

```bash
cd analyses/cell-typing/

singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/cell-typing \
  ../../images/scphylogenomics.sif \
  Rscript --vanilla 01-consolidate-samples.R \
    --project AML --integrate FALSE --normalize_method LogNormalize --components 20 --resolution 1.4

singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/cell-typing \
  ../../images/scphylogenomics.sif \
  Rscript --vanilla 02-annotate-cell-types.R \
    --project AML --assay RNA --annot_method mapping --ref_data inputs/AML-Le1-cell-annotation.tsv.gz
```

###### ploidy-inference module

```bash
cd analyses/ploidy-inference/

singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/ploidy-inference \
  ../../images/scphylogenomics.sif \
  Rscript --vanilla 01-classify-cell-ploidy.R \
    --project AML --annot_object AML-mapping-annotations.rds --nthreads 30
```

###### snv-calling module

```bash
cd analyses/snv-calling/

singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/snv-calling \
  ../../images/scphylogenomics.sif \
  Rscript --vanilla 01-get-cell-type-barcodes.R --project AML --cell_types all
  
singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/snv-calling \
  ../../images/scphylogenomics.sif \
  python 02-snv-calling.py --project AML --num_threads 30
  
 singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/snv-calling \
  ../../images/scphylogenomics.sif \
  Rscript --vanilla 03-summarize-sample-snv-calls.R --project AML --sample Le1 --cell_type all-cell-types 
```

###### phylogeny-inference module

```bash
cd analyses/phylogeny-inference/

singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/phylogeny-inference \
  ../../images/scphylogenomics.sif \
  python3 01-phylogeny-inference.py --project AML --gap_site_thresh 0.995 --gap_seq_thresh 0.985 --method iqtree --keep_temp
  
singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/phylogeny-inference \
  ../../images/scphylogenomics.sif \
  python3 02-wnmf-clonal-clustering.py --project AML --threads 30 --snv_threshold 0.05 --cell_threshold 0.05 --min_k 2 --min_k 10
  
singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/phylogeny-inference \
  ../../images/scphylogenomics.sif \
  python3 02-wnmf-clonal-clustering.py --project AML --threads 30 --snv_threshold 0.05 --cell_threshold 0.05 --force_k 3  
  
singularity exec \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  --pwd /home/rstudio/scPhylogenomics/analyses/phylogeny-inference \
  ../../images/scphylogenomics.sif \
  Rscript --vanilla 03-infer-clonal-phylogeny.R --project AML --rescale --percentile_outlier 0.99 --annot_type Clone
```

##### 4). Optional: Interactive use with Singularity shell

For development, debugging, or exploratory analysis, users may wish to open an interactive shell inside the Singularity container. This is equivalent to running an interactive session with Docker:

```bash
singularity shell \
  --bind /mnt/isilon/dbhi_bfx/wafulae/scPhylogenomics:/home/rstudio/scPhylogenomics \
  images/scphylogenomics.sif
```

From inside the container, navigate to the project directory:

```bash
cd /home/rstudio/scPhylogenomics
```

You can then run module scripts, R, Python, or bash commands interactively. Changes to files in the bound directory will be reflected on the host system.

## License
scPhylogenomics is distributed under the GNU GPL v3. For more information, see [license(LICENSE)
