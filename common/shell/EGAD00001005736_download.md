# Downloading and processing data from the European Genome-Phenome Archive (EGA)

## Overview

This document describes the steps to download **selected samples** from the EGA dataset **EGAD00001005736** (Young et al. 2021 — *Profiling molecular heterogeneity in human primary microglia*) using the **pyEGA3** download client on an Ubuntu instance, and convert the downloaded BAM files to paired-end FASTQ format.

We download a subset of 21 selected samples (BAM files) from this dataset, identified by their EGAF accession IDs:

| EGAF ID | EGAF ID | EGAF ID |
|---------|---------|---------|
| EGAF00003091586 | EGAF00003091723 | EGAF00003091811 |
| EGAF00003091616 | EGAF00003091744 | EGAF00003091817 |
| EGAF00003091638 | EGAF00003091748 | EGAF00003091826 |
| EGAF00003091649 | EGAF00003091755 | EGAF00003091833 |
| EGAF00003091659 | EGAF00003091775 | EGAF00003091854 |
| EGAF00003091662 | EGAF00003091776 | EGAF00003091860 |
|  | EGAF00003091795 | EGAF00003091861 |
|  |  | EGAF00003091866 |

The sample selection is based on the metadata file `exvivo_metadata_final.csv`.

## Prerequisites

- An EGA account with **data access** granted to the dataset(s) of interest
- An Ubuntu instance (tested on Ubuntu 24.04, Oracle Cloud)
- Sufficient disk space (the selected 21 samples total ~57 GB as BAMs, plus ~57 GB for FASTQs, plus temporary space during conversion)

## 1. Install system dependencies

```bash
sudo apt update
sudo apt install python3-pip
sudo apt install python3.12-venv
sudo apt install samtools
```

## 2. Set up a Python virtual environment

A virtual environment is required because Ubuntu 24.04 prevents system-wide pip installs.

```bash
# Create the virtual environment
python3 -m venv ~/ega_env

# Activate it (required every time you open a new terminal session)
source ~/ega_env/bin/activate
```

## 3. Install pyEGA3 and fix tqdm compatibility

```bash
# Install the EGA download client
pip install pyega3

# Upgrade tqdm to fix a known progress bar bug in pyega3 5.2.0
pip install --upgrade tqdm
```

## 4. Create a credentials file

Create `credentials.json` with your EGA **download** credentials (not submission credentials):

```json
{
    "username": "your.email@example.com",
    "password": "YOUR_EGA_PASSWORD"
}
```

## 5. Verify access

```bash
# List datasets you have access to
pyega3 -cf credentials.json datasets

# List files within a specific dataset
pyega3 -cf credentials.json files EGAD00001005736
```

## 6. Download files

### Single file test

Test with a single file first to confirm everything works:

```bash
pyega3 -cf credentials.json fetch EGAF00003091616 --output-dir /data/Young_2021_EGAD00001005736/Raw_data/EGAF00003091616/
```

### Batch download script

The script below downloads the 21 selected EGAF files in parallel (14 at a time), each into its own subdirectory.

Save as `download_ega.sh`:

```bash
#!/bin/bash

# Activate the virtual environment
source ~/ega_env/bin/activate

# EGAF IDs to download (21 selected samples from Young et al. 2021, EGAD00001005736)
EGAF_LIST=(
EGAF00003091616 EGAF00003091659 EGAF00003091775 EGAF00003091662
EGAF00003091826 EGAF00003091795 EGAF00003091860 EGAF00003091748
EGAF00003091811 EGAF00003091586 EGAF00003091854 EGAF00003091833
EGAF00003091649 EGAF00003091866 EGAF00003091861 EGAF00003091723
EGAF00003091755 EGAF00003091744 EGAF00003091817 EGAF00003091638
EGAF00003091776
)

# Path to your credentials file (use absolute path)
CREDENTIALS="/data/geo_scripts/credentials.json"

# Base output directory
BASEDIR="/data/Young_2021_EGAD00001005736/Raw_data"

# Download function for a single EGAF file
download_one() {
    local EGAF=$1
    local OUTDIR="${BASEDIR}/${EGAF}"
    mkdir -p "$OUTDIR"
    echo "[$(date)] Starting $EGAF"
    pyega3 -cf "$CREDENTIALS" fetch "$EGAF" --output-dir "$OUTDIR"
    echo "[$(date)] Finished $EGAF"
}

export -f download_one
export CREDENTIALS BASEDIR

# Run downloads in parallel (-P 14 = 14 concurrent downloads, one per CPU)
# Reduce -P if you encounter 500 errors from EGA
printf '%s\n' "${EGAF_LIST[@]}" | xargs -P 14 -I {} bash -c 'download_one "$@"' _ {}
```

Make executable and run in the background:

```bash
chmod +x download_ega.sh
nohup ./download_ega.sh > download_log.txt 2>&1 &
```

Monitor progress:

```bash
tail -f download_log.txt
```

### Clean up nested directories

pyEGA3 creates an extra subdirectory inside each output folder. To flatten:

```bash
cd /data/Young_2021_EGAD00001005736/Raw_data/
for dir in EGAF*/; do
    mv "$dir"/EGAF*/* "$dir"/ 2>/dev/null
    rm -r "$dir"/EGAF*/
done
```

## 7. Verify downloads

```bash
# Count BAM files (expect 21 selected samples)
find /data/Young_2021_EGAD00001005736/Raw_data/ -name "*.bam" | wc -l

# Check for empty files
find /data/Young_2021_EGAD00001005736/Raw_data/ -name "*.bam" -empty
```

Note: pyEGA3 automatically verifies MD5 checksums during download. If a file shows "Download complete" in the log, it has passed integrity checks.

## 8. Convert BAM to FASTQ

The downloaded BAM files are paired-end. We convert them to paired FASTQ files using `samtools`.

Save as `bam_to_fastq.sh`:

```bash
#!/bin/bash

source ~/ega_env/bin/activate

BASEDIR="/data/Young_2021_EGAD00001005736/Raw_data"

convert_one() {
    local BAM=$1
    local DIR=$(dirname "$BAM")
    local NAME=$(basename "$BAM" .bam)

    echo "[$(date)] Processing $NAME..."

    # Sort by read name (required for proper paired-end extraction)
    samtools sort -n -@ 2 -o "${DIR}/${NAME}.namesorted.bam" "$BAM"

    # Convert to paired-end FASTQ
    samtools fastq -@ 2 \
        -1 "${DIR}/${NAME}_1.fastq.gz" \
        -2 "${DIR}/${NAME}_2.fastq.gz" \
        -0 /dev/null \
        -s /dev/null \
        "${DIR}/${NAME}.namesorted.bam"

    # Remove intermediate name-sorted BAM
    rm "${DIR}/${NAME}.namesorted.bam"

    echo "[$(date)] Finished $NAME"
}

export -f convert_one

# Run 7 conversions in parallel (2 threads each = 14 CPUs total)
find "$BASEDIR" -name "*.bam" ! -name "*.namesorted.bam" | \
    xargs -P 7 -I {} bash -c 'convert_one "$@"' _ {}
```

Make executable and run in the background:

```bash
chmod +x bam_to_fastq.sh
nohup ./bam_to_fastq.sh > fastq_conversion_log.txt 2>&1 &
```

Monitor progress:

```bash
tail -f fastq_conversion_log.txt
```

> **Note on disk space:** The name-sorting step creates temporary files roughly the same size as the original BAMs. With 7 parallel jobs, ensure you have at least ~40 GB of free space. If disk space is tight, reduce parallelism to `-P 3`.

## 9. Verify FASTQ conversion

```bash
# Count FASTQ files (expect 42: 21 pairs × 2)
find /data/Young_2021_EGAD00001005736/Raw_data/ -name "*.fastq.gz" | wc -l

# Check no files are empty
find /data/Young_2021_EGAD00001005736/Raw_data/ -name "*.fastq.gz" -empty

# Verify gzip integrity
find /data/Young_2021_EGAD00001005736/Raw_data/ -name "*.fastq.gz" | \
    xargs -P 14 -I {} bash -c 'gzip -t "$1" && echo "OK: $1" || echo "CORRUPT: $1"' _ {}

# Check paired files have equal read counts
for dir in /data/Young_2021_EGAD00001005736/Raw_data/EGAF*/; do
    R1=$(zcat "$dir"/*_1.fastq.gz | wc -l)
    R2=$(zcat "$dir"/*_2.fastq.gz | wc -l)
    NAME=$(basename "$dir")
    if [ "$R1" -eq "$R2" ]; then
        echo "OK: $NAME  ($((R1/4)) read pairs)"
    else
        echo "MISMATCH: $NAME  R1=$((R1/4)) R2=$((R2/4))"
    fi
done
```

## Final directory structure

```
/data/Young_2021_EGAD00001005736/
└── Raw_data/
    ├── EGAF00003091586/
    │   ├── MGB_073.bam
    │   ├── MGB_073.bam.md5
    │   ├── MGB_073_1.fastq.gz
    │   └── MGB_073_2.fastq.gz
    ├── EGAF00003091616/
    │   ├── MGB_017.bam
    │   ├── MGB_017.bam.md5
    │   ├── MGB_017_1.fastq.gz
    │   └── MGB_017_2.fastq.gz
    └── ...
```

## Useful commands

```bash
# Re-activate venv after a new SSH session
source ~/ega_env/bin/activate

# Check running downloads
ps aux | grep pyega3

# Kill all running downloads
pkill -f pyega3

# Kill all running BAM conversions
pkill -f samtools

# Check disk space
df -h /data/
```

## References

- [pyEGA3 GitHub repository](https://github.com/EGA-archive/ega-download-client)
- [EGA PyEGA3 documentation](https://ega-archive.org/access/download/files/pyega3/)
- [EGA video tutorial](https://embl-ebi.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=be79bb93-1737-4f95-b80f-ab4300aa6f5a)
