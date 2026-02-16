# ============================================================
# rstudio-server_microglial-identity: RStudio Server +
# Bioconductor environment for microglial RNA-seq analysis
# Docker Hub: stodo3569/rstudio-server_microglial-identity:0.0
# ============================================================
# Base:  bioconductor/bioconductor_docker (rocker/rstudio +
#        Bioconductor system deps pre-installed)
# R/Bioc packages: tidyverse, DESeq2, edgeR, sva, DaMiRseq,
#        WGCNA, org.Hs.eg.db, googledrive, patchwork, gt, ...
#
# Build:
#   docker build -t stodo3569/rstudio-server_microglial-identity:0.0 \
#     -f Dockerfile_rstudio-server_microglial-identity .
#
# Run:
#   docker run -d -p 8787:8787 \
#     -e DISABLE_AUTH=true \
#     -v "$(pwd)/silken-forest-424410-r0-96bfeb3bd96c.json:/home/rstudio/key.json" \
#     stodo3569/rstudio-server_microglial-identity:0.0
#
# Push:
#   docker login
#   docker push stodo3569/rstudio-server_microglial-identity:0.0
#
# Access RStudio at http://localhost:8787
# ============================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL maintainer="stodo3569"
LABEL description="RStudio Server environment for microglial-identity RNA-seq analysis"

# ── System dependencies ──────────────────────────────────────
# Most are already in the base image; add any extras needed by
# WGCNA (GO.db compile), gt (chromote/webshot), googledrive, etc.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcairo2-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

# ── Fix xfun/knitr version mismatch in base image ───────────
# The base image ships an older xfun missing the 'attr' export,
# which breaks Hmisc, WGCNA, and DaMiRseq at install/load time.
RUN R -e "install.packages(c('xfun', 'knitr'), repos = 'https://cloud.r-project.org')"

# ── Bioconductor packages ───────────────────────────────────
# Install Bioconductor packages first (longer compile times)
RUN R -e "\
  BiocManager::install(c( \
    'DESeq2', \
    'edgeR', \
    'sva', \
    'DaMiRseq', \
    'org.Hs.eg.db', \
    'SummarizedExperiment', \
    'impute', \
    'preprocessCore', \
    'GO.db' \
  ), ask = FALSE, update = FALSE) \
"

# ── CRAN packages ────────────────────────────────────────────
RUN R -e "\
  install.packages(c( \
    'tidyverse', \
    'WGCNA', \
    'Hmisc', \
    'patchwork', \
    'gt', \
    'broom', \
    'googledrive', \
    'svglite' \
  ), repos = 'https://cloud.r-project.org', Ncpus = parallel::detectCores()) \
"

# ── Verify key packages load ────────────────────────────────
RUN R -e "\
  suppressPackageStartupMessages({ \
    library(tidyverse);   library(DESeq2);   library(edgeR); \
    library(sva);         library(DaMiRseq); library(WGCNA); \
    library(org.Hs.eg.db); library(patchwork); library(gt); \
    library(googledrive); library(Hmisc);    library(broom); \
    library(svglite) \
  }); \
  cat('All packages loaded successfully.\n') \
"

# ── Working directory & permissions ──────────────────────────
RUN mkdir -p /home/rstudio/all/Data /home/rstudio/all/figures \
    && chown -R rstudio:rstudio /home/rstudio

WORKDIR /home/rstudio

EXPOSE 8787
