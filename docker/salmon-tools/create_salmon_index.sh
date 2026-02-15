#!/bin/bash

################################################################################
# Salmon Transcriptome Index Builder
################################################################################
#
# Purpose:
#   Downloads GENCODE reference files and builds a Salmon index with decoy
#   sequences (the full genome) for accurate RNA-seq quantification. Uses the
#   "gentrome" approach, which concatenates the transcriptome with the primary
#   genome assembly so that reads originating from intergenic or intronic
#   regions are captured by decoy sequences rather than incorrectly assigned
#   to transcripts.
#
# Reference: GENCODE Release 49 / GRCh38 (human)
# Docker image: stodo3569/salmon-tools:0.0
#
# Usage:
#   Run from inside the directory that will hold the index, e.g.:
#
#   mkdir -p /data/salmon_transcriptome && cd /data/salmon_transcriptome
#   bash /data/<scripts_dir>/salmon_index_build.sh
#
# Pipeline stages:
#   1. Download    — fetch transcriptome FASTA and genome FASTA from GENCODE
#   2. Decoys      — extract chromosome names from the genome for decoy list
#   3. Gentrome    — concatenate transcriptome + genome into a single FASTA
#   4. Index       — run salmon index with decoy-aware mapping (via Docker)
#
# Outputs (written to the current working directory):
#   gencode.v49.transcripts.fa.gz           — GENCODE transcriptome (cDNA)
#   GRCh38.primary_assembly.genome.fa.gz    — GRCh38 primary assembly genome
#   GRCh38.primary_assembly.genome.fa       — decompressed genome (temporary)
#   decoys.txt                              — chromosome/scaffold names for Salmon
#   gentrome.fa.gz                          — concatenated transcriptome + genome
#   salmon_index/                           — final Salmon index directory
#
# Dependencies:
#   wget:          downloading reference files from GENCODE/EBI
#   gunzip:        decompressing genome FASTA for decoy extraction
#   Docker:        running the salmon-tools container
#   salmon (>=1.10): provided by the stodo3569/salmon-tools:0.0 image
#
# Resource notes:
#   - Downloads total ~1 GB compressed (~3.2 GB decompressed genome)
#   - Index build requires ~80 GB RAM and takes 30–60 min with 12 threads
#   - Final salmon_index/ directory is ~17 GB
#   - Adjust -p (threads) in Step 4 to match your available CPUs
#
# Note on /data:
#   The Docker container expects a host directory mounted at /data.
#   All reference files and the output index are read from and written
#   to /data/salmon_transcriptome/ within that mount point.
#
################################################################################

set -e
set -u
set -o pipefail

################################################################################
# Step 1 — Download reference files from GENCODE (Release 49 / GRCh38)
################################################################################
# The transcriptome FASTA contains all annotated transcript sequences (cDNA).
# The primary assembly genome excludes alt/patch scaffolds, which avoids
# mapping ambiguity from redundant sequences.

echo "Downloading GENCODE v49 transcriptome (cDNA)..."
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz

echo "Downloading GRCh38 primary assembly genome..."
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz

################################################################################
# Step 2 — Create the decoys.txt file
################################################################################
# Salmon's decoy-aware mode requires a list of sequence names that should be
# treated as decoys (i.e., the genome chromosomes/scaffolds). Reads mapping
# better to these decoys are discarded rather than being spuriously assigned
# to transcripts. The -k flag keeps the original .gz file intact.

echo "Extracting chromosome names for decoy list..."
gunzip -k GRCh38.primary_assembly.genome.fa.gz
grep "^>" GRCh38.primary_assembly.genome.fa \
  | cut -d " " -f 1 \
  | sed 's/>//g' \
  > decoys.txt

echo "  Found $(wc -l < decoys.txt) decoy sequences"

################################################################################
# Step 3 — Create the "gentrome" file
################################################################################
# The gentrome is a single FASTA that places the transcriptome BEFORE the
# genome. Order matters: Salmon indexes transcripts first, then genome
# sequences are loaded as decoys. Both inputs remain gzipped — Salmon's
# indexer handles compressed input natively.

echo "Creating gentrome (transcriptome + genome)..."
cat gencode.v49.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz

################################################################################
# Step 4 — Build the Salmon index
################################################################################
# Flags:
#   -t   Input FASTA (gentrome: transcriptome + genome)
#   -d   Decoy sequence names (chromosome list from Step 2)
#   -p   Number of threads (adjust to your machine; 12 is a good default)
#   -i   Output index directory
#   --gencode   Trims version suffixes from GENCODE transcript IDs
#               (e.g. ENST00000456328.2 → ENST00000456328) for cleaner
#               downstream ID matching

echo "Building Salmon index (this may take 30–60 minutes)..."
docker run --rm \
  -v /data:/data \
  stodo3569/salmon-tools:0.0 \
  salmon index \
    -t /data/salmon_transcriptome/gentrome.fa.gz \
    -d /data/salmon_transcriptome/decoys.txt \
    -p 12 \
    -i /data/salmon_transcriptome/salmon_index \
    --gencode

echo "Done. Salmon index is at: /data/salmon_transcriptome/salmon_index/"
