FROM --platform=linux/amd64 rocker/tidyverse:4.4.0
LABEL maintainer="ewafula@gmail.edu"
WORKDIR /rocker-build/

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

### Install apt-getable packages to start
#########################################

# Installing all apt required packages at once
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    jq \
    libgmp3-dev \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev \
    libpoppler-cpp-dev \
    libglpk-dev \
    libncurses5 \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libsqlite3-dev \
    libgdbm-dev \
    libdb5.3-dev \
    libbz2-dev \
    libexpat1-dev \
    liblzma-dev \
    libffi-dev \
    libuuid1 \
    gsl-bin \
    libgsl-dev \
    wget \
    xorg \
    zlib1g-dev \
    sendmail \
    mailutils \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and install Python 3.11
RUN cd /usr/src && \
    wget https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tgz && \
    tar xzf Python-3.11.0.tgz && \
    cd Python-3.11.0 && \
    ./configure --enable-optimizations && \
    make altinstall && \
    rm -rf /usr/src/Python-3.11.0.tgz

# Setup the default python commands to use Python 3.11
RUN ln -s /usr/local/bin/python3.11 /usr/local/bin/python3 && \
    ln -s /usr/local/bin/python3.11 /usr/local/bin/python
RUN python3 -m pip install --upgrade pip

# Set working directory
WORKDIR /home/rstudio

# Install python packages
##########################

# Install python3 tools and ALL dependencies
RUN pip3 install \
    "matplotlib==3.10.1" \
    "numpy==2.2.4" \
    "pandas==2.2.3" \
    "scanpy==1.11.1" \
    "scikit-learn==1.5.2" \
    "scipy==1.15.2" \
    "seaborn==0.13.2" \
    "multiqc==1.28" \
    "tqdm==4.67.1" \
    "biopython==1.85" \
    "fastcluster==1.3.0" \
    "setuptools==78.1.0" \
    "umap-learn==0.5.7" \
    "pysam==0.23.3" \
    "utils==1.0.2" \
    "wheel==0.45.1" \
    && rm -rf /root/.cache/pip/wheels

# Standalone tools and libraries
################################

# Add Cell Ranger
RUN wget -O cellranger-10.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1768972784&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=ol8sP1WNEWi2iew-UjxRt4dE7k39Bmf2zLdnuM-5BN7awPtoFySz8WVLKR9MYBL6nQP4aHJRXWPD~FlJ070agwFpLMqSxQ0MJWsshaczUNzK6IkSsl9f4wtEvOiP4erJZM2AQMoyNHNSZ0ZksFcvZq-gA~1IrJVFVeaJLlEToimDQC91levhwNW2-xgukACW4P1NOLFgpYbw-YkgAd2T4FjcDfOUVO~D-zC1tAnf6zBB152ZHhw4Kwsg~~kflxXGdBJv9qdM8FuirQHeZfFrjNgp1rYt21orSazThdOnE8wEHw0pBUkMbQU3867W62ux4jxwQKz-Sb8fynmDvRgPCQ__" && \
    tar -zxvf cellranger-10.0.0.tar.gz && rm -f ccellranger-10.0.0.tar.gz && \
    mv cellranger-10.0.0 /usr/local/bin/ && \
    ln -s /usr/local/bin/cellranger-10.0.0/cellranger /usr/local/bin/cellranger


# Add FastTree
RUN wget -O v2.2.0.tar.gz https://github.com/morgannprice/fasttree/archive/refs/tags/v2.2.0.tar.gz && \
    tar -zxvf v2.2.0.tar.gz && rm -f v2.2.0.tar.gz && \
    mv fasttree-2.2.0/FastTree /usr/local/bin/FastTree && rm -rf fasttree-2.2.0

# Add iqtree3
RUN wget -O iqtree-3.0.1-Linux.tar.gz https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Linux.tar.gz && \
    tar -zxvf iqtree-3.0.1-Linux.tar.gz && rm -f iqtree-3.0.1-Linux.tar.gz && \
    mv iqtree-3.0.1-Linux/bin/iqtree3* /usr/local/bin/ && rm -rf iqtree-3.0.1-Linux

#### R packages
###############

# Define GITHUB_PAT argument to avoid rate limits
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.20')"

# Install matrixStats 1.4.1
RUN R -e "remotes::install_github('HenrikBengtsson/matrixStats', ref='develop')"

# Install R packages
RUN R -e 'BiocManager::install(c( \
  "multtest", \
  "DropletUtils", \
  "SingleCellExperiment", \
  "SummarizedExperiment", \
  "glmGamPoi", \
  "celldex", \
  "SingleR", \
  "scran", \
  "Biobase", \
  "GenomicRanges", \
  "GenomeInfoDb", \
  "GenomicAlignments", \
  "IRanges", \
  "S4Vectors", \
  "BiocGenerics", \
  "Matrix", \
  "MatrixGenerics", \
  "DirichletMultinomial", \
  "TFBSTools", \
  "biomaRt", \
  "BSgenome.Hsapiens.UCSC.hg38", \
  "EnsDb.Hsapiens.v86", \
  "harmony", \
  "SeuratObject", \
  "sctransform", \
  "future", \
  "devtools", \
  "tidyverse", \
  "SoupX", \
  "patchwork", \
  "Seurat", \
  "flexmix", \
  "optparse", \
  "glue", \
  "data.table", \
  "scCustomize", \
  "viridis", \
  "ComplexHeatmap", \
  "dittoSeq", \
  "Nebulosa", \
  "ggpubr", \
  "hdf5r", \
  "liger", \
  "ape", \
  "adephylo", \
  "phytools", \
  "symphony", \
  "RColorBrewer", \
  "AUCell", \
  "UCell", \
  "doMC", \
  "BiocNeighbors", \
  "uwot" \
  ))'

# Install DoubletFinder
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"

# install signac
RUN R -e "remotes::install_github('stuart-lab/signac', ref = 'develop')"

# Install Azimuth
RUN R -e "remotes::install_github('satijalab/seurat-data')"
RUN R -e "remotes::install_github('satijalab/azimuth', ref = 'master')"

# Intall ggtree
RUN R -e "remotes::install_github('YuLab-SMU/yulab.utils')"
RUN R -e "remotes::install_github('YuLab-SMU/ggfun')"
RUN R -e "BiocManager::install(c('treeio','ggtree'), update=TRUE, ask=FALSE)"

# Install SeuratWrappers
RUN R -e "options(timeout=9999999)"
RUN R -e "remotes::install_github('satijalab/seurat-wrappers')"

# Install CopyKAT
RUN R -e "remotes::install_github('navinlabcode/copykat')"

# Install phylotools
RUN R -e "remotes::install_github('helixcn/phylotools', build_vignettes = TRUE)"

# Install BoneMarrowMap
RUN R -e "remotes::install_github('andygxzeng/BoneMarrowMap')"

#########################################################
# Install Miniconda and Custom Environments
#########################################################
WORKDIR /opt

# Install Miniconda to /opt/conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.5.0-0-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Add conda to path temporarily for environment creation
ENV PATH="/opt/conda/bin:$PATH"

# 1. Copy the YAML files exported from your local machine
COPY cellsnp_env.yml /tmp/cellsnp_env.yml
COPY snpmanifold_env.yml /tmp/snpmanifold_env.yml

# 2. Create the environments from the YAML files
# This creates isolated environments in /opt/conda/envs/cellsnp_env and /opt/conda/envs/snpmanifold_env
RUN conda env create -f /tmp/cellsnp_env.yml && \
    conda env create -f /tmp/snpmanifold_env.yml && \
    conda clean -afy && \
    rm /tmp/cellsnp_env.yml /tmp/snpmanifold_env.yml

# 3. Ensure permissions are open (so non-root users can activate if needed)
RUN chmod -R 755 /opt/conda

##############################
# Reset shell/env for the rest
##############################
# Reset PATH so system Python 3.11 is default, not Conda's base python
ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# Reset the frontend variable
ENV DEBIAN_FRONTEND=

WORKDIR /rocker-build/