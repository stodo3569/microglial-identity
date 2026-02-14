# ============================================================
# fastp-tools: FASTQ quality control and trimming environment
# Docker Hub: stodo3569/fastp_tools:0.0
# ============================================================
# Includes: fastp (>=0.23.0), FastQC (>=0.12.0), MultiQC (>=1.14),
#           parallel, pigz, curl
# Build:    docker build -t stodo3569/fastp_tools:0.0 .
#
# Usage:
#   docker run --rm \
#     -u "$(id -u):$(id -g)" \
#     -v /data:/data \
#     stodo3569/fastp_tools:0.0 \
#     bash /data/scripts/<your_script>.sh [OPTIONS]
#
# Note on /data:
#   The container expects a host directory mounted to /data.
#   Input FASTQ files (from geo-tools downloads) and scripts are
#   read from /data, and all QC reports and trimmed outputs are
#   written back to /data.
# ============================================================

FROM mambaorg/micromamba:1.5-jammy

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    parallel \
    pigz \
    curl \
    default-jre \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Install bioinformatics tools via conda
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    fastp>=0.23.0 \
    fastqc>=0.12.0 \
    multiqc>=1.14 \
    && micromamba clean --all --yes

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
