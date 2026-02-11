FROM --platform=linux/amd64 rocker/tidyverse:4.4.0
LABEL maintainer="ewafula@gmail.edu"

# Set initial working directory for build artifacts
WORKDIR /rocker-build/

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

### 1. System Dependencies
##########################
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
    cmake \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

### 2. Python 3.11 Setup
########################
RUN cd /usr/src && \
    wget https://www.python.org/ftp/python/3.11.0/Python-3.11.0.tgz && \
    tar xzf Python-3.11.0.tgz && \
    cd Python-3.11.0 && \
    ./configure --enable-optimizations && \
    make altinstall && \
    rm -rf /usr/src/Python-3.11.0.tgz

# Setup default python links
RUN ln -s /usr/local/bin/python3.11 /usr/local/bin/python3 && \
    ln -s /usr/local/bin/python3.11 /usr/local/bin/python
RUN python3 -m pip install --upgrade pip

# Install Python packages
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

### 3. Standalone Tools
#######################

# Cell Ranger
RUN wget -O cellranger-10.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1769654122&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=DYiC5v~8mxDGHwg75DE26edwR3YpbyUWGhL59bLq0LoQEW2AHVMciEyxJp7p2yyh6L3PQMKhBO6JzX5tpCafOk-p-km1Cut4xMphJutS4hd1V2hmoUAVGDEGPwafUiTzEkcDGz9mm3iWemWuJPH5v6tQDZ5VCvGvJJDZWM7D~hHN58yr5G~HF-4j0fkLfReLw9R6lbLDLgu3b~zlMBgBeA5dYK2oEyr2cF~zI-gxx4Enj6Q9D0shzIbNLINcYy1mpc~VB7Hx9jifrDxw4MoBTFG8Q0EMNMZ40B2JhX1v~L2wC1qeG9MUkVmRdL11lb-KLPuyfDnhUh0~lyBGZdCbZw__" && \
    tar -zxvf cellranger-10.0.0.tar.gz && rm -f cellranger-10.0.0.tar.gz && \
    mv cellranger-10.0.0 /usr/local/bin/ && \
    ln -s /usr/local/bin/cellranger-10.0.0/cellranger /usr/local/bin/cellranger

# FastTree
RUN wget -O v2.2.0.tar.gz https://github.com/morgannprice/fasttree/archive/refs/tags/v2.2.0.tar.gz && \
    tar -zxvf v2.2.0.tar.gz && rm -f v2.2.0.tar.gz && \
    mv fasttree-2.2.0/FastTree /usr/local/bin/FastTree && rm -rf fasttree-2.2.0

# iqtree3
RUN wget -O iqtree-3.0.1-Linux.tar.gz https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Linux.tar.gz && \
    tar -zxvf iqtree-3.0.1-Linux.tar.gz && rm -f iqtree-3.0.1-Linux.tar.gz && \
    mv iqtree-3.0.1-Linux/bin/iqtree3* /usr/local/bin/ && rm -rf iqtree-3.0.1-Linux

### 4. R Packages
#################

ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT

# Configure Bioconductor & Basic Setup
RUN R -e "options(repos = BiocManager::repositories()); install.packages('BiocManager'); BiocManager::install(version = '3.20')"

# Install matrixStats (dev version)
RUN R -e "remotes::install_github('HenrikBengtsson/matrixStats', ref='develop', upgrade='never')"

# Bulk Install Bioc/CRAN Packages
# Note: Added 'upgrade=FALSE' to prevent interactive prompts
RUN R -e "BiocManager::install(c( \
  'multtest', 'DropletUtils', 'SingleCellExperiment', 'SummarizedExperiment', \
  'glmGamPoi', 'celldex', 'SingleR', 'scran', 'Biobase', 'GenomicRanges', \
  'GenomeInfoDb', 'GenomicAlignments', 'IRanges', 'S4Vectors', 'BiocGenerics', \
  'Matrix', 'MatrixGenerics', 'DirichletMultinomial', 'TFBSTools', 'biomaRt', \
  'BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86', 'harmony', 'SeuratObject', \
  'sctransform', 'future', 'devtools', 'tidyverse', 'SoupX', 'patchwork', \
  'Seurat', 'flexmix', 'optparse', 'glue', 'data.table', 'scCustomize', \
  'viridis', 'ComplexHeatmap', 'dittoSeq', 'Nebulosa', 'ggpubr', 'hdf5r', \
  'liger', 'ape', 'adephylo', 'phytools', 'symphony', 'RColorBrewer', \
  'AUCell', 'UCell', 'doMC', 'BiocNeighbors', 'uwot' \
  ), update=FALSE, ask=FALSE)"

# GitHub Packages
# IMPORTANT: 'upgrade="never"' is added to all these to prevent build failures/hangs
RUN R -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade='never')"
RUN R -e "remotes::install_github('stuart-lab/signac', ref = 'develop', upgrade='never')"
RUN R -e "remotes::install_github('satijalab/seurat-data', upgrade='never')"
RUN R -e "remotes::install_github('satijalab/azimuth', ref = 'master', upgrade='never')"

# ggtree stack
RUN R -e "remotes::install_github('YuLab-SMU/yulab.utils', upgrade='never')"
RUN R -e "remotes::install_github('YuLab-SMU/ggfun', upgrade='never')"
RUN R -e "BiocManager::install(c('treeio','ggtree'), update=TRUE, ask=FALSE)"

# SeuratWrappers (FIXED)
# Combined timeout and install into one block so settings persist
RUN R -e "options(timeout=9999999); remotes::install_github('satijalab/seurat-wrappers', upgrade='never')"

# Other GitHub Tools
RUN R -e "remotes::install_github('navinlabcode/copykat', upgrade='never')"
RUN R -e "remotes::install_github('helixcn/phylotools', build_vignettes = TRUE, upgrade='never')"
RUN R -e "remotes::install_github('andygxzeng/BoneMarrowMap', upgrade='never')"

### 5. Miniconda & Environments
###############################
WORKDIR /opt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.5.0-0-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

COPY cellsnp_env.yml /tmp/cellsnp_env.yml
COPY snpmanifold_env.yml /tmp/snpmanifold_env.yml

# Set massive timeouts for Pip and Conda to prevent read errors
# PIP_DEFAULT_TIMEOUT=1000 (seconds) prevents time associated errors y
# CONDA_REMOTE_READ_TIMEOUT_SECS helps with the conda solver
ENV PIP_DEFAULT_TIMEOUT=1000 \
    CONDA_REMOTE_READ_TIMEOUT_SECS=1000

RUN conda env create -f /tmp/cellsnp_env.yml && \
    conda env create -f /tmp/snpmanifold_env.yml && \
    conda clean -afy && \
    rm /tmp/cellsnp_env.yml /tmp/snpmanifold_env.yml

RUN chmod -R 755 /opt/conda

# Reset Environment
ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV DEBIAN_FRONTEND=
WORKDIR /home/rstudio

### 6. FINAL VERIFICATION (The "Smoke Test")
############################################
# This block will cause the build to fail if critical tools are missing.
# We test the tools that are most likely to fail silently.

RUN R -e "\
    check_pkgs <- c('SeuratWrappers', 'Seurat', 'Signac', 'DoubletFinder', 'Azimuth', 'copykat'); \
    missing <- c(); \
    for (pkg in check_pkgs) { \
      if (!require(pkg, character.only = TRUE)) { \
        missing <- c(missing, pkg); \
      } \
    }; \
    if (length(missing) > 0) { \
      stop(paste('FATAL: The following packages failed to install:', paste(missing, collapse=', '))); \
    }"

# Check Command Line tools
RUN cellranger --version && \
    echo "Verification Complete: All critical tools are present."