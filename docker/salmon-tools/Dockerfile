# ============================================================
# salmon-tools: Salmon transcript quantification environment
#               with pre-built decoy-aware transcriptome index
# Docker Hub: stodo3569/salmon-tools:0.1
# ============================================================
# Includes: salmon (>=1.10.0), parallel, procps
# Index:    GENCODE v49 / GRCh38 decoy-aware (built at image time)
# Build:    docker build -t stodo3569/salmon-tools:0.1 .
#
# Usage:
#   # Using the built-in GENCODE v49 index (default):
#   docker run --rm \
#     -u "$(id -u):$(id -g)" \
#     -v /data:/data \
#     stodo3569/salmon-tools:0.1 \
#     bash /data/geo_scripts/salmon_quant_docker.sh \
#       --index /opt/salmon_index \
#       --input-file /data/datasets.txt
#
#   # Using a custom index instead:
#   docker run --rm \
#     -u "$(id -u):$(id -g)" \
#     -v /data:/data \
#     stodo3569/salmon-tools:0.1 \
#     bash /data/geo_scripts/salmon_quant_docker.sh \
#       --index /data/my_custom_index \
#       --input-file /data/datasets.txt
#
# Pre-built index details:
#   Location inside container:  /opt/salmon_index
#   Reference:                  GENCODE Release 49 (GRCh38 primary assembly)
#   Transcriptome:              gencode.v49.transcripts.fa.gz
#   Decoy sequences:            GRCh38.primary_assembly.genome.fa.gz
#   Salmon flags used:          --gencode  (GENCODE-style header parsing)
#
#   The decoy-aware index concatenates the transcriptome and genome
#   into a single "gentrome" FASTA, with chromosome names listed in
#   decoys.txt.  This lets salmon map reads that originate from
#   genomic (non-transcriptomic) loci without incorrectly assigning
#   them to transcripts, reducing false-positive quantification.
#
# Note on /data:
#   The container expects a host directory mounted to /data.
#   Input trimmed FASTQ files (from fastp_tools processing) are
#   read from /data, and all quantification outputs (quant.sf,
#   aux_info/, logs/) are written back to
#   /data/<Output_Directory>/Aligned_data/.
#
# Note on image size:
#   This image is substantially larger than salmon-tools:0.0
#   because it ships the pre-built salmon index (~15-17 GB).
#   If image size is a concern, use salmon-tools:0.0 and supply
#   your own index via --index.
# ============================================================

FROM mambaorg/micromamba:1.5-jammy

USER root

# Install system dependencies
# - parallel:        GNU parallel for concurrent sample processing
# - procps:          /proc utilities used by RAM monitoring in the quant script
# - wget:            download GENCODE reference files during index build
# - ca-certificates: trusted root CAs required for HTTPS downloads (wget)
RUN apt-get update && apt-get install -y --no-install-recommends \
    parallel \
    procps \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Install bioinformatics tools via conda
RUN micromamba install -y -n base -c bioconda -c conda-forge \
    salmon>=1.10.0 \
    && micromamba clean --all --yes

# --------------------------------------------------------------------------
# Build the GENCODE v49 decoy-aware salmon index
# --------------------------------------------------------------------------
# This single RUN block downloads, builds, and cleans up in one layer
# so that the intermediate files (genome FASTA, gentrome, etc.) do not
# inflate the final image.  Only /opt/salmon_index is retained.
#
# Steps mirror the manual build procedure:
#   1. Download transcriptome and genome from GENCODE (Release 49 / GRCh38)
#   2. Extract chromosome names from the genome → decoys.txt
#   3. Concatenate transcriptome + genome → gentrome.fa.gz
#   4. Run salmon index with --gencode and decoy sequences
#   5. Remove all intermediate files
# --------------------------------------------------------------------------

# The conda base environment bin must be on PATH for root to find salmon
ENV PATH="/opt/conda/bin:$PATH"

USER root

RUN mkdir -p /opt/salmon_build && cd /opt/salmon_build \
    \
    # 1. Download GENCODE v49 reference files
    && echo ">>> Downloading GENCODE v49 transcriptome..." \
    && wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz \
    \
    && echo ">>> Downloading GRCh38 primary assembly genome..." \
    && wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz \
    \
    # 2. Create decoys.txt (chromosome names from the genome)
    && echo ">>> Extracting decoy sequence names..." \
    && gunzip -k GRCh38.primary_assembly.genome.fa.gz \
    && grep "^>" GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 | sed 's/>//g' > decoys.txt \
    && rm GRCh38.primary_assembly.genome.fa \
    \
    # 3. Create gentrome (transcriptome + genome concatenation)
    && echo ">>> Building gentrome FASTA..." \
    && cat gencode.v49.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz \
    \
    # 4. Build salmon index
    && echo ">>> Building salmon index (this will take a while)..." \
    && salmon index \
         -t gentrome.fa.gz \
         -d decoys.txt \
         -p $(nproc) \
         -i /opt/salmon_index \
         --gencode \
    \
    # 5. Clean up all intermediate files
    && echo ">>> Cleaning up build artifacts..." \
    && cd / && rm -rf /opt/salmon_build

USER $MAMBA_USER

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["bash"]
