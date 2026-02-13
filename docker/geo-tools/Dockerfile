# ============================================================
# geo-tools: SRA/GEO dataset download environment
# Docker Hub: stodo3569/geo-tools:0.0
# ============================================================
# Includes: sra-tools (>=3.0.0), pysradb, parallel, pigz, curl
# Build:    docker build -t stodo3569/geo-tools:0.0 .
#
# Usage:
#   docker run --rm \
#     -u "$(id -u):$(id -g)" \
#     --shm-size=80g \
#     -v /data:/data \
#     stodo3569/geo-tools:0.0 \
#     bash /data/scripts/<your_script>.sh --input-file /data/<your_input>.txt
#
# Note on /data:
#   The container expects a host directory mounted to /data.
#   This serves as the shared working directory for input files,
#   scripts, and downloaded outputs. SRA toolkit is configured
#   to cache prefetch data to /data/ncbi.
# ============================================================

FROM mambaorg/micromamba:1.5-jammy

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    parallel \
    pigz \
    curl \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    sra-tools>=3.0.0 \
    pysradb \
    && micromamba clean --all --yes

# Set up SRA toolkit config — downloads go to /data/ncbi
RUN mkdir -p /home/$MAMBA_USER/.ncbi && \
    echo '/repository/user/main/public/root = "/data/ncbi"' > /home/$MAMBA_USER/.ncbi/user-settings.mkfg && \
    echo '/repository/user/cache-disabled = "true"' >> /home/$MAMBA_USER/.ncbi/user-settings.mkfg

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
