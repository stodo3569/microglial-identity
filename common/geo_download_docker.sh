#!/bin/bash

################################################################################
# GEO/SRA Data Download Pipeline
################################################################################
#
# Purpose:
#   Downloads RNA-seq (or other NGS) data from NCBI GEO/SRA, validates file
#   integrity, extracts FASTQ files, and compresses them. Supports both GSE
#   and SRP accessions with optional sample filtering (GSM/SRX).
#
# Docker image: stodo3569/geo-tools:0.0
#
# Usage:
#   docker run --rm \
#     -u "$(id -u):$(id -g)" \
#     --shm-size=80g \
#     -v /data:/data \
#     stodo3569/geo-tools:0.0 \
#     bash /data/<scripts_dir>/geo_download_docker.sh [OPTIONS]
#
# Modes:
#   # Single dataset — download all samples
#   bash geo_download_docker.sh "Study_Name" "GSE123456" "all"
#
#   # Single dataset — specific GSM samples
#   bash geo_download_docker.sh "Study_Name" "GSE123456" "GSM111,GSM222"
#
#   # Batch mode — multiple datasets from input file
#   bash geo_download_docker.sh --input-file datasets.txt
#
#   # Batch mode — parallel processing (2 datasets simultaneously)
#   bash geo_download_docker.sh --input-file datasets.txt --parallel 2
#
# Input file format (tab-separated):
#   Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC_file]
#   Lines with the same Output_Directory+Accession are automatically aggregated.
#
# Pipeline stages:
#   1. Prefetch   — downloads .sra/.sralite files from NCBI
#   2. Validate   — runs vdb-validate to check file integrity
#   3. Extract    — converts SRA to FASTQ via fasterq-dump + pigz compression
#   Failed samples are automatically retried sequentially, then with disk temp.
#
# Outputs (written to /data/<Output_Directory>/):
#   Raw_data/<GSM>/<SRR>.fastq.gz   — compressed FASTQ files
#   download_summary.txt             — summary of results
#   error_log.txt                    — detailed error log
#   prefetch_failed.txt              — samples that failed download
#   validation_failed.txt            — samples that failed integrity check
#   extraction_failed.txt            — samples that failed FASTQ extraction
#   successful_samples.txt           — samples that completed all stages
#
# Dependencies (provided by stodo3569/geo-tools Docker image):
#   sra-tools (>=3.0.0): prefetch, vdb-validate, fasterq-dump
#   pysradb:             automatic GSE→GSM→SRR accession resolution
#   GNU parallel:        parallel download and extraction
#   pigz:                parallel gzip compression
#   curl:                network utilities
#
# Resource handling:
#   CPU and memory are detected automatically and allocated adaptively
#   across pipeline stages. --shm-size controls the RAM disk used for
#   fast temporary extraction (recommended: ≥50% of system RAM).
#
# Note on /data:
#   The container expects a host directory mounted at /data. All inputs,
#   scripts, and outputs are read from and written to this mount point.
#   SRA toolkit caches prefetched data to /data/ncbi.
#
################################################################################

set -e
set -u
set -o pipefail

################################################################################
# System Resource Detection
################################################################################

detect_system_resources() {
    if [[ -f /proc/cpuinfo ]]; then
        TOTAL_CPUS=$(grep -c ^processor /proc/cpuinfo)
    else
        TOTAL_CPUS=$(nproc 2>/dev/null || echo 4)
    fi
    
    # NEW: Adjust for parallel datasets to prevent CPU oversubscription
    local available_cpus=$TOTAL_CPUS
    if [[ -n "${PARALLEL_DATASETS:-}" ]] && [[ $PARALLEL_DATASETS -gt 1 ]]; then
        available_cpus=$((TOTAL_CPUS / PARALLEL_DATASETS))
        log_message "CPU allocation adjusted for parallel processing:"
        log_message "  Total system CPUs: $TOTAL_CPUS"
        log_message "  Parallel datasets: $PARALLEL_DATASETS"
        log_message "  CPUs per dataset: ~${available_cpus}"
    fi
    
    # Smart CPU reservation based on allocated CPU count
    if [[ $available_cpus -gt 8 ]]; then
        CPUS=$((available_cpus - 2))  # Reserve 2 CPUs on large allocations
    elif [[ $available_cpus -gt 4 ]]; then
        CPUS=$((available_cpus - 1))  # Reserve 1 CPU on medium allocations
    else
        CPUS=$available_cpus  # Use all allocated CPUs on small allocations
    fi
    
    # Ensure minimum of 2 CPUs per dataset
    if [[ $CPUS -lt 2 ]]; then
        CPUS=2
        log_warning "Very low CPU allocation detected - setting minimum of 2 CPUs"
    fi
    
    # Detect total RAM
    if [[ -f /proc/meminfo ]]; then
        TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        TOTAL_MEM_GB=$((TOTAL_MEM_KB / 1024 / 1024))
        
        # Detect available RAM (more accurate than total)
        AVAILABLE_MEM_KB=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print $2}')
        if [[ -n "$AVAILABLE_MEM_KB" ]]; then
            AVAILABLE_MEM_GB=$((AVAILABLE_MEM_KB / 1024 / 1024))
        else
            # Fallback: estimate 70% of total is available
            AVAILABLE_MEM_GB=$((TOTAL_MEM_GB * 70 / 100))
        fi
        
        # Adjust available memory for parallel datasets
        if [[ -n "${PARALLEL_DATASETS:-}" ]] && [[ $PARALLEL_DATASETS -gt 1 ]]; then
            AVAILABLE_MEM_GB=$((AVAILABLE_MEM_GB / PARALLEL_DATASETS))
            log_message "  Memory per dataset: ~${AVAILABLE_MEM_GB}GB"
        fi
    else
        TOTAL_MEM_GB=8
        AVAILABLE_MEM_GB=6
    fi
    
    # Check for RAM disk availability
    RAM_DISK_AVAILABLE=false
    RAM_DISK_SIZE=0
    if [[ -d "/dev/shm" ]] && [[ -w "/dev/shm" ]]; then
        RAM_DISK_SIZE=$(df -BG /dev/shm 2>/dev/null | tail -1 | awk '{print $2}' | sed 's/G//' || echo 0)
        if [[ $RAM_DISK_SIZE -gt 2 ]]; then
            RAM_DISK_AVAILABLE=true
        fi
    fi
    
    log_message "System Resources Detected:"
    log_message "  Total CPUs: $TOTAL_CPUS"
    log_message "  CPUs for this dataset: $CPUS"
    log_message "  Total Memory: ${TOTAL_MEM_GB}GB"
    log_message "  Available Memory: ${AVAILABLE_MEM_GB}GB"
    if [[ "$RAM_DISK_AVAILABLE" == true ]]; then
        log_message "  RAM Disk: ${RAM_DISK_SIZE}GB available at /dev/shm ✓"
    else
        log_message "  RAM Disk: Not available"
    fi
}

################################################################################
# Adaptive Resource Calculation
################################################################################

calculate_optimal_parallelism() {
    local stage=$1
    local num_samples=${2:-10}
    
    case $stage in
        prefetch)
            # Network/disk I/O bound - scale conservatively
            if [[ $CPUS -ge 32 ]]; then
                echo 8
            elif [[ $CPUS -ge 16 ]]; then
                echo 6
            elif [[ $CPUS -ge 8 ]]; then
                echo 4
            elif [[ $CPUS -ge 4 ]]; then
                echo 3
            else
                echo 2
            fi
            ;;
            
        validation)
            # CPU-bound - use all available CPUs
            echo $CPUS
            ;;
            
        extraction)
            # SIMPLE STRATEGY: 1 file per 2 CPUs
            # 14 CPUs ÷ 2 = 7 parallel files
            local cpus_per_file=2
            local optimal_jobs=$((CPUS / cpus_per_file))
            
            # Ensure reasonable bounds
            [[ $optimal_jobs -lt 1 ]] && optimal_jobs=1
            
            # Don't spawn more jobs than samples
            if [[ $num_samples -lt $optimal_jobs ]]; then
                echo $num_samples
            else
                echo $optimal_jobs
            fi
            ;;
            
        *)
            echo 2
            ;;
    esac
}

calculate_extraction_threads() {
    local parallel_jobs=$1
    
    # 1 thread for extraction (1 of the 2 CPUs per file)
    echo 1
}

calculate_compression_threads() {
    local parallel_jobs=$1
    local extraction_threads=$2
    
    # 1 thread for compression (the other 1 of the 2 CPUs per file)
    echo 1
}

calculate_fasterq_mem_for_jobs() {
    local parallel_jobs=$1
    local available_mem=$2  # in GB
    
    # Calculate memory per job dynamically
    # Reserve 10% for system overhead
    local usable_mem=$((available_mem * 90 / 100))
    local mem_per_job=$((usable_mem / parallel_jobs))
    
    # Apply sensible limits
    if [[ $mem_per_job -gt 64 ]]; then
        mem_per_job=64  # Cap at 64GB per job (fasterq-dump doesn't need more)
    elif [[ $mem_per_job -lt 2 ]]; then
        mem_per_job=2   # Minimum 2GB per job
    fi
    
    echo "${mem_per_job}G"
}

calculate_fasterq_bufsize() {
    local mem_value=$1  # e.g., "16G"
    local mem_gb=${mem_value%G}  # Remove 'G' suffix
    
    # Bufsize should be ~10% of memory
    local bufsize_mb=$((mem_gb * 100))
    
    # Cap at 4GB
    if [[ $bufsize_mb -gt 4096 ]]; then
        bufsize_mb=4096
    elif [[ $bufsize_mb -lt 256 ]]; then
        bufsize_mb=256
    fi
    
    echo "${bufsize_mb}MB"
}

calculate_fasterq_cache() {
    local mem_value=$1  # e.g., "16G"
    local mem_gb=${mem_value%G}
    
    # Cache should be ~5% of memory
    local cache_mb=$((mem_gb * 50))
    
    # Cap at 2GB
    if [[ $cache_mb -gt 2048 ]]; then
        cache_mb=2048
    elif [[ $cache_mb -lt 128 ]]; then
        cache_mb=128
    fi
    
    echo "${cache_mb}MB"
}

optimize_fasterq_settings() {
    # DEPRECATED - now calculated dynamically per extraction phase
    # This is kept for backward compatibility but values will be overridden
    FASTERQ_MEM="16G"
    FASTERQ_BUFSIZE="2048MB"
    FASTERQ_CACHE="1024MB"
}

################################################################################
# Configuration
################################################################################

BASE_PATH="${GEO_BASE_PATH:-/data}"
MAX_SIZE="u"
PARALLEL_DATASETS=1

OUTPUT_DIR=""
ACCESSION=""
GSM_SAMPLES="all"
NGC_FILE=""
ACCESSION_TYPE=""
MANUAL_MODE=false
INPUT_FILE=""

CPUS=4
TOTAL_CPUS=4
TOTAL_MEM_GB=8
AVAILABLE_MEM_GB=6
RAM_DISK_AVAILABLE=false
RAM_DISK_SIZE=0

# Adaptive fasterq-dump settings (calculated based on available RAM)
FASTERQ_MEM="2G"
FASTERQ_BUFSIZE="512MB"
FASTERQ_CACHE="256MB"

# Error tracking files
ERROR_LOG=""
PREFETCH_FAILED=""
VALIDATION_FAILED=""
EXTRACTION_FAILED=""
SUCCESSFUL_SAMPLES=""

################################################################################
# Logging Functions
################################################################################

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log_error() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*"
    echo "$msg" >&2
    if [[ -n "$ERROR_LOG" ]]; then
        echo "$msg" >> "$ERROR_LOG"
    fi
}

log_warning() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $*"
    echo "$msg"
    if [[ -n "$ERROR_LOG" ]]; then
        echo "$msg" >> "$ERROR_LOG"
    fi
}

log_success() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✓ $*"
}

log_step_error() {
    local step="$1"
    local sample="$2"
    local run="$3"
    local error_msg="${4:-Unknown error}"
    
    local msg="[$step] Failed: $sample / $run - $error_msg"
    log_error "$msg"
}

################################################################################
# Error Tracking Functions
################################################################################

initialize_error_tracking() {
    ERROR_LOG="error_log.txt"
    PREFETCH_FAILED="prefetch_failed.txt"
    VALIDATION_FAILED="validation_failed.txt"
    EXTRACTION_FAILED="extraction_failed.txt"
    SUCCESSFUL_SAMPLES="successful_samples.txt"
    
    # Clear existing error files
    > "$ERROR_LOG"
    > "$PREFETCH_FAILED"
    > "$VALIDATION_FAILED"
    > "$EXTRACTION_FAILED"
    > "$SUCCESSFUL_SAMPLES"
    
    log_message "Error tracking initialized"
    log_message "  Error log: $ERROR_LOG"
    log_message "  Failed files will be tracked per stage"
}

mark_prefetch_failed() {
    local sample="$1"
    local run="$2"
    echo -e "${sample}\t${run}" >> "$PREFETCH_FAILED"
}

mark_validation_failed() {
    local sample="$1"
    local run="$2"
    echo -e "${sample}\t${run}" >> "$VALIDATION_FAILED"
}

mark_extraction_failed() {
    local sample="$1"
    local run="$2"
    echo -e "${sample}\t${run}" >> "$EXTRACTION_FAILED"
}

mark_successful() {
    local sample="$1"
    local run="$2"
    echo -e "${sample}\t${run}" >> "$SUCCESSFUL_SAMPLES"
}

is_sample_failed() {
    local sample="$1"
    local run="$2"
    
    # Check if sample failed in any previous stage
    if [[ -f "$PREFETCH_FAILED" ]] && grep -q "^${sample}[[:space:]]${run}$" "$PREFETCH_FAILED"; then
        return 0
    fi
    
    if [[ -f "$VALIDATION_FAILED" ]] && grep -q "^${sample}[[:space:]]${run}$" "$VALIDATION_FAILED"; then
        return 0
    fi
    
    return 1
}

################################################################################
# Argument Parsing
################################################################################

show_help() {
    cat << EOF
GEO/SRA Download Script - Docker/Oracle Version with Error Logging

USAGE:
  Single dataset:
    $0 OUTPUT_DIR ACCESSION GSM_SAMPLES [NGC_FILE]
    
  Multiple datasets from file:
    $0 --input-file FILE [--parallel N]

EXAMPLES:
  # Download all samples
  $0 "Friedman_2019" "GSE125050" "all"
  
  # Download specific GSM samples (comma-separated)
  $0 "Friedman_2019" "GSE125050" "GSM3559136,GSM3559137,GSM3559138"
  
  # Download GSMs from file
  $0 "Friedman_2019" "GSE125050" "gsm_list.txt"
  
  # With NGC file
  $0 "Gosselin_2017" "SRP115940" "all" "/data/prj_30034_D31481.ngc"
  
  # Batch processing (sequential)
  $0 --input-file datasets.txt
  
  # Batch processing (parallel - 2 datasets at once)
  $0 --input-file datasets.txt --parallel 2

PARALLEL PROCESSING:
  --parallel N will download N datasets simultaneously
  CPUs will be automatically divided among datasets to prevent oversubscription
  
  Examples for different system sizes:
    8 CPUs:  --parallel 1 (each dataset uses ~6 CPUs)
    16 CPUs: --parallel 2 (each dataset uses ~7 CPUs)
    32 CPUs: --parallel 4 (each dataset uses ~7 CPUs)

ERROR TRACKING:
  All errors are logged to error_log.txt
  Failed files are tracked in:
    - prefetch_failed.txt
    - validation_failed.txt
    - extraction_failed.txt
  
  Only files that pass validation will be extracted.
  Processing continues for other files even when some fail.

INPUT FILE FORMAT (One GSM per line):
  Output_Directory<TAB>GSE_accession<TAB>GSM_sample<TAB>[NGC_file]
  
  Multiple lines with the same Output_Directory will be combined.
  
  Examples:
    # Download all samples
    Friedman_2019	GSE125050	all
    
    # Download specific samples (one per line)
    Matlik_2023	GSE227729	GSM7106344
    Matlik_2023	GSE227729	GSM7106387
    Matlik_2023	GSE227729	GSM7106390
    
    # With NGC file
    Gosselin_2017	SRP115940	all	/data/prj_30034_D31481.ngc

ALTERNATIVE FORMAT (Comma-separated):
  You can also use comma-separated GSMs in one line:
    Matlik_2023	GSE227729	GSM7106344,GSM7106387,GSM7106390

GSM LIST FILE FORMAT (one GSM per line):
    GSM7106344
    GSM7106387
    GSM7106390

EOF
}

parse_arguments() {
    if [[ $# -eq 0 ]]; then
        show_help
        exit 1
    fi
    
    if [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
        show_help
        exit 0
    fi
    
    if [[ "$1" == "--input-file" ]]; then
        if [[ $# -lt 2 ]]; then
            log_error "Option --input-file requires a filename"
            exit 1
        fi
        
        INPUT_FILE="$2"
        shift 2
        
        while [[ $# -gt 0 ]]; do
            case $1 in
                --parallel)
                    if [[ $# -lt 2 ]]; then
                        log_error "Option --parallel requires a number"
                        exit 1
                    fi
                    PARALLEL_DATASETS="$2"
                    shift 2
                    ;;
                *)
                    log_error "Unknown option: $1"
                    exit 1
                    ;;
            esac
        done
        
        return 0
    fi
    
    # Single dataset mode
    if [[ $# -lt 3 ]]; then
        log_error "Insufficient arguments"
        log_error "Usage: $0 OUTPUT_DIR ACCESSION GSM_SAMPLES [NGC_FILE]"
        exit 1
    fi
    
    OUTPUT_DIR="$1"
    ACCESSION="$2"
    GSM_SAMPLES="$3"
    NGC_FILE="${4:-}"
    
    # Auto-detect accession type
    if [[ $ACCESSION =~ ^GSE ]]; then
        ACCESSION_TYPE="gse"
    elif [[ $ACCESSION =~ ^SRP ]]; then
        ACCESSION_TYPE="srp"
    else
        log_error "Unknown accession type: $ACCESSION (expected GSE* or SRP*)"
        exit 1
    fi
}

################################################################################
# Dependency Checking
################################################################################

check_dependencies() {
    log_message "Checking dependencies..."
    
    local deps=(prefetch vdb-validate fasterq-dump parallel pigz)
    local missing=()
    
    if ! command -v pysradb &> /dev/null; then
        log_warning "pysradb not found - will use manual mode if needed"
    fi
    
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            missing+=("$dep")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing required dependencies: ${missing[*]}"
        exit 1
    fi
    
    log_success "All required dependencies found"
}

################################################################################
# Pre-flight Checks
################################################################################

check_disk_space() {
    local num_samples=$1
    
    # Estimate: ~3GB per sample (1.5GB SRA + 1.5GB FASTQ compressed)
    local estimated_gb=$((num_samples * 3))
    
    # Add 20% buffer
    local required_gb=$((estimated_gb + estimated_gb / 5))
    
    local available_gb=$(df -BG "$BASE_PATH" | tail -1 | awk '{print $4}' | sed 's/G//')
    
    log_message "Disk space check:"
    log_message "  Location: $BASE_PATH"
    log_message "  Available: ${available_gb}GB"
    log_message "  Estimated needed: ~${required_gb}GB (${num_samples} samples)"
    
    if [[ $available_gb -lt $required_gb ]]; then
        log_warning "Low disk space detected!"
        log_warning "  Available: ${available_gb}GB"
        log_warning "  Estimated need: ${required_gb}GB"
        log_warning "Proceeding anyway, but monitor disk usage..."
    else
        log_success "Sufficient disk space available"
    fi
}

################################################################################
# Directory Setup
################################################################################

create_output_structure() {
    log_message "Creating output directory structure..."
    
    local full_path="$BASE_PATH/$OUTPUT_DIR"
    mkdir -p "$full_path"
    cd "$full_path" || exit 1
    
    log_message "Working directory: $(pwd)"
    
    # Initialize error tracking in the output directory
    initialize_error_tracking
}

################################################################################
# GSM Sample Filtering
################################################################################

parse_gsm_samples() {
    local gsm_input="$1"
    
    # If "all" or empty, return empty (means download all)
    if [[ -z "$gsm_input" ]] || [[ "$gsm_input" == "all" ]]; then
        echo ""
        return 0
    fi
    
    # Check if it's a file
    if [[ -f "$gsm_input" ]]; then
        # Read GSM list from file
        grep "^GSM" "$gsm_input" | tr '\n' ',' | sed 's/,$//'
        return 0
    fi
    
    # Otherwise, treat as comma-separated list
    echo "$gsm_input"
}

create_filtered_gsm_file() {
    local full_gsm_file="$1"
    local filtered_file="$2"
    local gsm_filter="$3"
    
    # If no filter, copy all
    if [[ -z "$gsm_filter" ]]; then
        cp "$full_gsm_file" "$filtered_file"
        local total=$(tail -n +2 "$filtered_file" | wc -l)
        log_message "Selected all $total GSM samples"
        return 0
    fi
    
    # Parse comma-separated GSMs
    IFS=',' read -ra GSM_ARRAY <<< "$gsm_filter"
    
    # Create filtered file with header
    head -1 "$full_gsm_file" > "$filtered_file"
    
    local kept=0
    local total=$(tail -n +2 "$full_gsm_file" | wc -l)
    
    for gsm in "${GSM_ARRAY[@]}"; do
        gsm=$(echo "$gsm" | tr -d '[:space:]')
        
        if grep -q "[[:space:]]${gsm}$\|[[:space:]]${gsm}[[:space:]]" "$full_gsm_file"; then
            grep "[[:space:]]${gsm}$\|[[:space:]]${gsm}[[:space:]]" "$full_gsm_file" >> "$filtered_file"
            ((kept++))
        else
            log_warning "GSM not found: $gsm"
        fi
    done
    
    if [[ $kept -eq 0 ]]; then
        log_error "No matching GSM samples found"
        return 1
    fi
    
    log_message "Selected $kept of $total GSM samples"
    return 0
}

filter_srr_mappings() {
    local gsm_file="$1"
    local full_srr_file="$2"
    local filtered_srr_file="$3"
    
    > "$filtered_srr_file"
    
    while read -r line; do
        gsm=$(echo "$line" | awk '{print $2}')
        grep "^${gsm}[[:space:]]" "$full_srr_file" >> "$filtered_srr_file" 2>/dev/null || true
    done < <(tail -n +2 "$gsm_file")
    
    local count=$(wc -l < "$filtered_srr_file")
    log_message "Kept $count SRR accessions for selected GSM samples"
    
    if [[ $count -eq 0 ]]; then
        log_error "No SRR accessions found for selected GSM samples"
        return 1
    fi
    
    return 0
}

################################################################################
# GSE Workflow Functions
################################################################################

get_sample_accessions_gse() {
    log_message "Getting GSM accessions for $ACCESSION..."
    
    # Parse GSM filter
    local gsm_filter=$(parse_gsm_samples "$GSM_SAMPLES")
    
    if [[ -n "$gsm_filter" ]]; then
        log_message "GSM filter active - will download only selected samples"
    else
        log_message "No GSM filter - will download all samples"
    fi
    
    # Check for manual files
    if [[ -f "tmp_gse_to_gsm.txt" ]] && [[ -f "tmp_gsm_to_srr.txt" ]]; then
        log_message "Found existing manual mapping files"
        
        if [[ -n "$gsm_filter" ]]; then
            mv tmp_gse_to_gsm.txt tmp_gse_to_gsm_full.txt
            create_filtered_gsm_file tmp_gse_to_gsm_full.txt tmp_gse_to_gsm.txt "$gsm_filter" || exit 1
            
            mv tmp_gsm_to_srr.txt tmp_gsm_to_srr_full.txt
            filter_srr_mappings tmp_gse_to_gsm.txt tmp_gsm_to_srr_full.txt tmp_gsm_to_srr.txt || exit 1
        fi
        
        MANUAL_MODE=true
        validate_manual_files
        return 0
    fi
    
    if ! command -v pysradb &> /dev/null; then
        log_error "pysradb not available and no manual files found"
        log_error "Please create tmp_gse_to_gsm.txt and tmp_gsm_to_srr.txt manually"
        exit 1
    fi
    
    log_message "Attempting automatic GSM fetch with pysradb..."
    
    if ! pysradb gse-to-gsm "$ACCESSION" > tmp_gse_to_gsm_full.txt 2>&1; then
        log_warning "pysradb gse-to-gsm failed for $ACCESSION"
        handle_pysradb_failure
        return 1
    fi
    
    local line_count=$(wc -l < tmp_gse_to_gsm_full.txt)
    if [[ $line_count -lt 2 ]]; then
        log_warning "pysradb returned insufficient data (only $line_count lines)"
        handle_pysradb_failure
        return 1
    fi
    
    # Apply filter
    if [[ -n "$gsm_filter" ]]; then
        create_filtered_gsm_file tmp_gse_to_gsm_full.txt tmp_gse_to_gsm.txt "$gsm_filter" || exit 1
    else
        mv tmp_gse_to_gsm_full.txt tmp_gse_to_gsm.txt
    fi
    
    # Create directories
    awk 'NR>1 {print $2}' tmp_gse_to_gsm.txt | \
        parallel --tmpdir . -j $CPUS -I% --max-args 1 mkdir -p "Raw_data/%"
    
    return 0
}

get_run_accessions_gse() {
    log_message "Getting SRR accessions..."
    
    if [[ -f "tmp_gsm_to_srr.txt" ]] && [[ -s "tmp_gsm_to_srr.txt" ]]; then
        log_message "Using existing tmp_gsm_to_srr.txt"
        return 0
    fi
    
    if ! command -v pysradb &> /dev/null; then
        log_error "pysradb not available and tmp_gsm_to_srr.txt not found"
        exit 1
    fi
    
    log_message "Fetching SRR accessions (this may take a while)..."
    
    > pre_tmp_gsm_to_srr.txt
    
    local gsm_count=0
    
    awk 'NR>1 {print $2}' tmp_gse_to_gsm.txt | while read -r gsm; do
        ((gsm_count++))
        log_message "  Processing $gsm ($gsm_count)..."
        
        if ! pysradb gsm-to-srr "$gsm" >> pre_tmp_gsm_to_srr.txt 2>&1; then
            log_warning "    Failed to get SRR for $gsm"
        fi
        
        sleep 2
    done
    
    if ! grep -E "^GSM" pre_tmp_gsm_to_srr.txt > tmp_gsm_to_srr.txt 2>/dev/null; then
        log_error "No valid GSM-to-SRR mappings found"
        handle_pysradb_failure
        return 1
    fi
    
    local mapping_count=$(wc -l < tmp_gsm_to_srr.txt)
    log_success "Successfully mapped $mapping_count samples to SRR accessions"
    
    if [[ $mapping_count -eq 0 ]]; then
        handle_pysradb_failure
        return 1
    fi
    
    return 0
}

handle_pysradb_failure() {
    log_warning "========================================="
    log_warning "AUTOMATIC MODE FAILED - MANUAL FILES NEEDED"
    log_warning "========================================="
    log_warning ""
    log_warning "pysradb cannot retrieve data for $ACCESSION"
    log_warning ""
    log_warning "NEXT STEPS:"
    log_warning "1. Create manual mapping files in $(pwd)/:"
    log_warning "   - tmp_gse_to_gsm.txt"
    log_warning "   - tmp_gsm_to_srr.txt"
    log_warning ""
    log_warning "2. Get data from:"
    log_warning "   GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$ACCESSION"
    log_warning "   SRA: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=$ACCESSION"
    log_warning ""
    log_warning "3. Rerun the same command - it will auto-detect manual files"
    log_warning "========================================="
    
    MANUAL_MODE=true
    exit 1
}

validate_manual_files() {
    log_message "Validating manual mapping files..."
    
    if [[ ! -f "tmp_gse_to_gsm.txt" ]] || [[ ! -f "tmp_gsm_to_srr.txt" ]]; then
        log_error "Manual mode requires both tmp_gse_to_gsm.txt and tmp_gsm_to_srr.txt"
        exit 1
    fi
    
    local gsm_count=$(tail -n +2 tmp_gse_to_gsm.txt | wc -l)
    local srr_count=$(wc -l < tmp_gsm_to_srr.txt)
    
    log_message "Found manual files:"
    log_message "  GSM accessions: $gsm_count"
    log_message "  SRR mappings: $srr_count"
    
    if [[ $gsm_count -eq 0 ]] || [[ $srr_count -eq 0 ]]; then
        log_error "Manual mapping files are empty"
        exit 1
    fi
    
    awk 'NR>1 {print $2}' tmp_gse_to_gsm.txt | \
        parallel --tmpdir . -j $CPUS -I% --max-args 1 mkdir -p "Raw_data/%"
}

################################################################################
# SRP/SRX Sample Filtering
################################################################################

parse_srx_samples() {
    local srx_input="$1"
    
    # If "all" or empty, return empty (means download all)
    if [[ -z "$srx_input" ]] || [[ "$srx_input" == "all" ]]; then
        echo ""
        return 0
    fi
    
    # Check if it's a file
    if [[ -f "$srx_input" ]]; then
        # Read SRX list from file
        grep "^SRX" "$srx_input" | tr '\n' ',' | sed 's/,$//'
        return 0
    fi
    
    # Otherwise, treat as comma-separated list
    echo "$srx_input"
}

detect_sample_type() {
    local sample_input="$1"
    
    # Check if input contains SRX accessions
    if [[ "$sample_input" == "all" ]] || [[ -z "$sample_input" ]]; then
        echo "unknown"
        return 0
    fi
    
    # Check first accession in comma-separated list
    local first_sample=$(echo "$sample_input" | cut -d',' -f1 | tr -d '[:space:]')
    
    if [[ "$first_sample" =~ ^SRX[0-9]+ ]]; then
        echo "srx"
    elif [[ "$first_sample" =~ ^GSM[0-9]+ ]]; then
        echo "gsm"
    else
        echo "unknown"
    fi
}

get_srx_to_srr_mappings() {
    local srp_accession="$1"
    
    log_message "Fetching SRX to SRR mappings for $srp_accession using pysradb..."
    
    if ! command -v pysradb &> /dev/null; then
        log_warning "pysradb not found - cannot auto-fetch SRX mappings"
        return 1
    fi
    
    # Try to get metadata from pysradb
    local temp_metadata="temp_srp_metadata.tsv"
    
    if pysradb metadata "$srp_accession" --detailed --saveto "$temp_metadata" 2>/dev/null; then
        log_success "Retrieved metadata from pysradb"
        
        # Extract SRX and SRR columns (columns may vary, so we search for them)
        # Format: experiment_accession (SRX) and run_accession (SRR)
        if awk -F'\t' 'NR==1 {
            for(i=1; i<=NF; i++) {
                if($i == "experiment_accession") srx_col=i;
                if($i == "run_accession") srr_col=i;
            }
        }
        NR>1 && srx_col && srr_col {
            print $srx_col "\t" $srr_col
        }' "$temp_metadata" > "full_srx_to_srr.tsv"; then
            
            local mapping_count=$(wc -l < "full_srx_to_srr.tsv")
            log_message "Extracted $mapping_count SRX to SRR mappings"
            rm -f "$temp_metadata"
            return 0
        else
            log_warning "Could not parse pysradb output"
            rm -f "$temp_metadata"
            return 1
        fi
    else
        log_warning "pysradb metadata retrieval failed for $srp_accession"
        rm -f "$temp_metadata"
        return 1
    fi
}

create_filtered_srx_file() {
    local full_srx_file="$1"
    local filtered_file="$2"
    local srx_filter="$3"
    
    # If no filter, copy all
    if [[ -z "$srx_filter" ]]; then
        cp "$full_srx_file" "$filtered_file"
        local total=$(wc -l < "$filtered_file")
        log_message "Selected all $total SRX samples"
        return 0
    fi
    
    # Parse comma-separated SRX accessions
    IFS=',' read -ra SRX_ARRAY <<< "$srx_filter"
    
    # Create filtered file
    > "$filtered_file"
    
    local kept=0
    local total=$(wc -l < "$full_srx_file")
    
    for srx in "${SRX_ARRAY[@]}"; do
        srx=$(echo "$srx" | tr -d '[:space:]')
        
        if grep -q "^${srx}[[:space:]]" "$full_srx_file"; then
            grep "^${srx}[[:space:]]" "$full_srx_file" >> "$filtered_file"
            ((kept++))
        else
            log_warning "SRX not found in mappings: $srx"
        fi
    done
    
    if [[ $kept -eq 0 ]]; then
        log_error "No matching SRX samples found"
        return 1
    fi
    
    log_message "Selected $kept of $total SRX samples"
    return 0
}

################################################################################
# SRP Workflow Functions
################################################################################

get_sample_accessions_srp() {
    log_message "Setting up SRP workflow for $ACCESSION..."
    
    # Detect if we're using SRX or traditional manual mode
    local sample_type=$(detect_sample_type "$GSM_SAMPLES")
    
    if [[ "$sample_type" == "srx" ]] || [[ "$sample_type" == "unknown" && "$GSM_SAMPLES" == "all" ]]; then
        # NEW MODE: Auto-fetch SRX mappings and filter
        log_message "SRX-based sample selection mode"
        
        # Parse SRX filter
        local srx_filter=$(parse_srx_samples "$GSM_SAMPLES")
        
        if [[ -n "$srx_filter" ]]; then
            log_message "SRX filter active - will download only selected samples: $srx_filter"
        else
            log_message "No SRX filter - will download all samples"
        fi
        
        # Try to get mappings automatically
        if get_srx_to_srr_mappings "$ACCESSION"; then
            # Filter the mappings if needed
            if create_filtered_srx_file "full_srx_to_srr.tsv" "srx_to_srr.tsv" "$srx_filter"; then
                log_success "Created filtered srx_to_srr.tsv mapping file"
            else
                log_error "Failed to filter SRX mappings"
                exit 1
            fi
        else
            # Fallback to manual mode
            log_error "Could not auto-fetch SRX mappings"
            log_error ""
            log_error "MANUAL MODE REQUIRED:"
            log_error "1. Create srx_to_srr.tsv in $(pwd)/"
            log_error "2. Format: SRX<tab>SRR (one per line)"
            log_error "3. Get data from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=$ACCESSION"
            log_error "4. Include only the SRX samples you want to download"
            log_error "5. Rerun the same command"
            exit 1
        fi
        
    else
        # LEGACY MODE: Check for existing manual file
        if [[ ! -f "srx_to_srr.tsv" ]]; then
            log_error "SRP workflow requires srx_to_srr.tsv file"
            log_error "Get data from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=$ACCESSION"
            exit 1
        fi
        
        log_message "Using existing srx_to_srr.tsv file (legacy mode)"
    fi
    
    local srx_count=$(wc -l < srx_to_srr.tsv)
    log_message "Processing $srx_count SRX samples"
    
    # Create directories for each SRX sample
    awk '{print $1}' srx_to_srr.tsv | \
        parallel --tmpdir . -j $CPUS -I% --max-args 1 mkdir -p "Raw_data/%"
    
    MANUAL_MODE=true
}

################################################################################
# Download Functions with Error Tracking
################################################################################

prefetch_data() {
    local sample_col=$1
    local run_col=$2
    local mapping_file=$3
    
    log_message "========================================="
    log_message "STEP 1: Prefetching SRA files"
    log_message "========================================="
    
    local ngc_arg=""
    if [[ -n "$NGC_FILE" ]]; then
        if [[ ! -f "$NGC_FILE" ]]; then
            log_error "NGC file specified but not found: $NGC_FILE"
            exit 1
        fi
        ngc_arg="--ngc $NGC_FILE"
        log_message "Using NGC authentication: $NGC_FILE"
    fi
    
    # Clean up any stale lock files from previous runs
    log_message "Cleaning up stale lock files..."
    find Raw_data -name "*.sra.lock" -mmin +60 -delete 2>/dev/null || true
    local locks_found=$(find Raw_data -name "*.sra.lock" 2>/dev/null | wc -l)
    if [[ $locks_found -gt 0 ]]; then
        log_warning "Found $locks_found recent lock files - may indicate concurrent downloads"
    fi
    
    # Create prefetch list
    local prefetch_list="samples_to_prefetch.txt"
    awk -v s="$sample_col" -v r="$run_col" '{print $s, $r}' "$mapping_file" > "$prefetch_list"
    
    local total_samples=$(wc -l < "$prefetch_list")
    
    # ADAPTIVE: Calculate optimal prefetch parallelism
    local prefetch_jobs=$(calculate_optimal_parallelism "prefetch" $total_samples)
    
    # Safety check for prefetch_jobs
    if [[ -z "$prefetch_jobs" ]] || ! [[ "$prefetch_jobs" =~ ^[0-9]+$ ]]; then
        log_warning "prefetch_jobs not set properly, defaulting to 4"
        prefetch_jobs=4
    fi
    
    log_message "Prefetching data for $total_samples samples..."
    log_message "  Parallel downloads: $prefetch_jobs (adaptive: $CPUS CPUs available)"
    log_message "  Errors will be logged, processing continues for other files"
    log_message ""
    
    # Export functions for parallel
    export -f log_message log_error log_step_error mark_prefetch_failed
    export ERROR_LOG PREFETCH_FAILED ngc_arg MAX_SIZE
    
    # Parallel prefetch
    cat "$prefetch_list" | \
        parallel --colsep ' ' -j $prefetch_jobs --line-buffer --tagstring '[{1}/{2}]' \
        '
        sample={1}
        run={2}
        
        echo "[$(date +%H:%M:%S)] Prefetching..."
        
        if prefetch {2} '"$ngc_arg"' --output-directory "Raw_data/${sample}" --max-size '"$MAX_SIZE"' 2>> "'"$ERROR_LOG"'"; then
            echo "[$(date +%H:%M:%S)] ✓ Complete"
        else
            echo "[$(date +%H:%M:%S)] ✗ Failed"
            log_step_error "PREFETCH" "${sample}" "{2}" "prefetch command failed"
            mark_prefetch_failed "${sample}" "{2}"
        fi
        '
    
    rm -f "$prefetch_list"
    
    # Count failures
    local failed_count=0
    if [[ -f "$PREFETCH_FAILED" ]]; then
        failed_count=$(wc -l < "$PREFETCH_FAILED")
    fi
    
    log_message ""
    if [[ $failed_count -gt 0 ]]; then
        log_warning "$failed_count sample(s) failed prefetch - see $PREFETCH_FAILED"
        log_message "Continuing with successful samples..."
    else
        log_success "All samples prefetched successfully"
    fi
}

validate_data() {
    local sample_col=$1
    local run_col=$2
    local mapping_file=$3
    
    log_message "========================================="
    log_message "STEP 2: Validating SRA files"
    log_message "========================================="
    log_message "Only samples that passed prefetch will be validated"
    
    # Build list of files to validate (excluding prefetch failures)
    local validate_list="samples_to_validate.txt"
    > "$validate_list"
    
    local skipped=0
    
    # Use process substitution instead of pipe to avoid subshell
    while read -r sample run; do
        if is_sample_failed "$sample" "$run"; then
            skipped=$((skipped + 1))
            log_message "  ⊗ Skipping validation (failed prefetch): $sample / $run"
        else
            # Support both .sra and .sralite (SRA Normalized Format)
            local sra_file="Raw_data/${sample}/${run}/${run}.sra"
            local sralite_file="Raw_data/${sample}/${run}/${run}.sralite"
            local file_to_validate=""
            
            if [[ -f "$sra_file" ]]; then
                file_to_validate="$sra_file"
            elif [[ -f "$sralite_file" ]]; then
                file_to_validate="$sralite_file"
            fi
            
            if [[ -n "$file_to_validate" ]]; then
                echo -e "${sample}\t${run}\t${file_to_validate}" >> "$validate_list"
            else
                log_step_error "VALIDATION" "$sample" "$run" "SRA file not found (checked .sra and .sralite)"
                mark_validation_failed "$sample" "$run"
            fi
        fi
    done < <(awk -v s="$sample_col" -v r="$run_col" '{print $s, $r}' "$mapping_file")
    
    local to_validate=$(wc -l < "$validate_list" 2>/dev/null || echo 0)
    
    if [[ $to_validate -eq 0 ]]; then
        log_warning "No files to validate"
        rm -f "$validate_list"
        return 0
    fi
    
    # ADAPTIVE: Use optimal parallelism for validation
    local validation_jobs
    validation_jobs=$(calculate_optimal_parallelism "validation" $to_validate) || {
        log_warning "Could not calculate optimal parallelism, using default"
        validation_jobs=$CPUS
    }
    
    log_message "Validating $to_validate files..."
    log_message "  Parallel validations: $validation_jobs (adaptive: using all $CPUS CPUs)"
    log_message ""
    
    # Export functions for parallel
    export -f log_message log_step_error mark_validation_failed
    export ERROR_LOG VALIDATION_FAILED
    
    cat "$validate_list" | \
        parallel --colsep '\t' -j $validation_jobs --line-buffer --tagstring '[{1}/{2}]' \
        '
        sample={1}
        run={2}
        sra_file={3}
        
        echo "[$(date +%H:%M:%S)] Validating..."
        
        if vdb-validate "$sra_file" 2>> "'"$ERROR_LOG"'"; then
            echo "[$(date +%H:%M:%S)] ✓ Valid"
        else
            echo "[$(date +%H:%M:%S)] ✗ Corrupt"
            log_step_error "VALIDATION" "${sample}" "{2}" "vdb-validate detected corruption"
            mark_validation_failed "${sample}" "{2}"
        fi
        '
    
    rm -f "$validate_list"
    
    # Count results
    local failed_count=0
    if [[ -f "$VALIDATION_FAILED" ]]; then
        failed_count=$(wc -l < "$VALIDATION_FAILED")
    fi
    
    local validated=$((to_validate - failed_count))
    
    log_message ""
    log_message "Validation Summary:"
    log_message "  Validated successfully: $validated"
    log_message "  Failed validation: $failed_count"
    log_message "  Skipped (prefetch failed): $skipped"
    
    if [[ $failed_count -gt 0 ]]; then
        log_warning "Files with validation errors will NOT be extracted"
        log_warning "See $VALIDATION_FAILED for list of failed files"
    fi
}

extract_and_compress() {
    local sample_col=$1
    local run_col=$2
    local mapping_file=$3
    
    log_message "========================================="
    log_message "STEP 3: Extracting and compressing FASTQ files"
    log_message "========================================="
    log_message "Only validated files will be extracted"
    
    local ngc_arg=""
    if [[ -n "$NGC_FILE" ]]; then
        ngc_arg="--ngc $NGC_FILE"
    fi
    
    local work_dir=$(pwd)
    
    # Create list of samples to process (excluding failed ones)
    local samples_to_process="$work_dir/samples_to_extract.txt"
    > "$samples_to_process"
    
    local skipped=0
    
    # Use process substitution to avoid subshell
    while read -r sample run; do
        if is_sample_failed "$sample" "$run"; then
            skipped=$((skipped + 1))
            log_message "  ⊗ Skipping (previous stage failed): $sample / $run"
        else
            # Support both .sra and .sralite formats
            local sra_file="$work_dir/Raw_data/${sample}/${run}/${run}.sra"
            local sralite_file="$work_dir/Raw_data/${sample}/${run}/${run}.sralite"
            
            if [[ -f "$sra_file" ]] || [[ -f "$sralite_file" ]]; then
                echo -e "${sample}\t${run}" >> "$samples_to_process"
            else
                log_step_error "EXTRACTION" "$sample" "$run" "SRA file not found (checked .sra and .sralite)"
                mark_extraction_failed "$sample" "$run"
            fi
        fi
    done < <(awk -v s="$sample_col" -v r="$run_col" '{print $s, $r}' "$mapping_file")
    
    local total_to_extract=$(wc -l < "$samples_to_process" 2>/dev/null || echo 0)
    
    if [[ $total_to_extract -eq 0 ]]; then
        log_warning "No files to extract"
        rm -rf "$temp_dir" 2>/dev/null
        return 0
    fi
    
    # ADAPTIVE: Calculate optimal extraction settings with error handling
    local parallel_jobs threads_per_job pigz_threads
    
    parallel_jobs=$(calculate_optimal_parallelism "extraction" $total_to_extract) || {
        log_warning "Could not calculate optimal parallelism, using default"
        parallel_jobs=2
    }
    
    threads_per_job=$(calculate_extraction_threads $parallel_jobs) || {
        log_warning "Could not calculate extraction threads, using default"
        threads_per_job=2
    }
    
    pigz_threads=$(calculate_compression_threads $parallel_jobs $threads_per_job) || {
        log_warning "Could not calculate compression threads, using default"
        pigz_threads=2
    }
    
    # ADAPTIVE: Optimize fasterq-dump memory settings
    optimize_fasterq_settings
    
    # Safety check: Ensure we don't oversubscribe CPUs
    local total_cpu_usage=$((parallel_jobs * (threads_per_job + pigz_threads)))
    if [[ $total_cpu_usage -gt $CPUS ]]; then
        log_warning "CPU oversubscription detected! Reducing parallel jobs..."
        log_warning "  Would use: $total_cpu_usage cores"
        log_warning "  Available: $CPUS cores"
        
        # Recalculate with lower parallelism
        parallel_jobs=$((CPUS / (threads_per_job + pigz_threads)))
        [[ $parallel_jobs -lt 1 ]] && parallel_jobs=1
        
        total_cpu_usage=$((parallel_jobs * (threads_per_job + pigz_threads)))
        log_message "  Adjusted to: $parallel_jobs parallel jobs (${total_cpu_usage} cores)"
    fi
    
    # CRITICAL: Cap parallel jobs based on RAM disk size if using /dev/shm
    # Each large sample needs ~10GB of temp space
    if [[ "$RAM_DISK_AVAILABLE" == true ]]; then
        local temp_space_per_job=10  # GB estimate for large samples
        local max_jobs_for_ram=$((RAM_DISK_SIZE / temp_space_per_job))
        
        if [[ $max_jobs_for_ram -lt 1 ]]; then
            max_jobs_for_ram=1
        fi
        
        if [[ $parallel_jobs -gt $max_jobs_for_ram ]]; then
            log_warning "Reducing parallel jobs from $parallel_jobs to $max_jobs_for_ram to fit in RAM disk"
            log_warning "  RAM disk: ${RAM_DISK_SIZE}GB, need ~${temp_space_per_job}GB per job"
            log_warning "  Estimated need: $((parallel_jobs * temp_space_per_job))GB, available: ${RAM_DISK_SIZE}GB"
            parallel_jobs=$max_jobs_for_ram
        fi
    fi
    
    # ADAPTIVE: Use RAM disk if available (major speed boost)
    local temp_dir
    if [[ "$RAM_DISK_AVAILABLE" == true ]]; then
        temp_dir="/dev/shm/sra_extraction_$$"
        log_message "Using RAM disk for temp: $temp_dir (${RAM_DISK_SIZE}GB available)"
    else
        temp_dir="$work_dir/tmp_extraction"
        log_message "Using disk for temp: $temp_dir"
    fi
    mkdir -p "$temp_dir"
    
    # DYNAMIC: Calculate optimal memory settings based on actual parallel jobs
    # This gives each job maximum memory instead of using pre-calculated static values
    FASTERQ_MEM=$(calculate_fasterq_mem_for_jobs $parallel_jobs $AVAILABLE_MEM_GB)
    FASTERQ_BUFSIZE=$(calculate_fasterq_bufsize "$FASTERQ_MEM")
    FASTERQ_CACHE=$(calculate_fasterq_cache "$FASTERQ_MEM")
    
    log_message "Adaptive extraction configuration:"
    log_message "  Samples to extract: $total_to_extract"
    log_message "  Parallel jobs: $parallel_jobs (based on $CPUS CPUs, ${AVAILABLE_MEM_GB}GB RAM)"
    log_message "  Extraction threads per job: $threads_per_job"
    log_message "  Compression threads per job: $pigz_threads"
    log_message "  Memory per job: $FASTERQ_MEM (dynamically calculated for $parallel_jobs jobs)"
    log_message "  Buffer size: $FASTERQ_BUFSIZE, Cache: $FASTERQ_CACHE"
    log_message "  Total CPU usage: ~$((parallel_jobs * (threads_per_job + pigz_threads)))/$CPUS cores"
    log_message "  Fast compression enabled (pigz -1)"
    log_message ""
    
    # Safety check: ensure parallel_jobs is set and is a number
    if [[ -z "$parallel_jobs" ]] || ! [[ "$parallel_jobs" =~ ^[0-9]+$ ]]; then
        log_warning "parallel_jobs not set properly, defaulting to 2"
        parallel_jobs=2
    fi
    
    log_message "Starting parallel extraction with -j $parallel_jobs (progress shown below)..."
    log_message ""
    
    # Export functions and variables for parallel
    export -f log_message log_error log_step_error mark_extraction_failed mark_successful is_sample_failed
    export ERROR_LOG EXTRACTION_FAILED SUCCESSFUL_SAMPLES PREFETCH_FAILED VALIDATION_FAILED
    export work_dir ngc_arg temp_dir threads_per_job pigz_threads
    export FASTERQ_MEM FASTERQ_BUFSIZE FASTERQ_CACHE
    
    # Process in parallel using GNU parallel - THIS IS THE KEY SPEEDUP
    cat "$samples_to_process" | \
        parallel --colsep '\t' -j $parallel_jobs --line-buffer --tagstring '[{1}/{2}]' \
        '
        sample={1}
        run={2}
        srr_dir="'"$work_dir"'/Raw_data/${sample}/${run}"
        
        echo "[$(date +%H:%M:%S)] Starting extraction..."
        
        # Determine which SRA format exists (.sra or .sralite)
        if [[ -f "$srr_dir/${run}.sra" ]]; then
            sra_input="${run}.sra"
        elif [[ -f "$srr_dir/${run}.sralite" ]]; then
            sra_input="${run}.sralite"
        else
            echo "[$(date +%H:%M:%S)] ✗ Failed - SRA file not found"
            log_step_error "EXTRACTION" "${sample}" "${run}" "No .sra or .sralite file found"
            mark_extraction_failed "${sample}" "${run}"
            exit 1
        fi
        
        # CRITICAL: Each job needs its own temp directory to avoid conflicts!
        job_temp_dir="'"$temp_dir"'/${sample}_${run}"
        mkdir -p "$job_temp_dir"
        
        # NOTE: We do NOT use --disk-limit parameter!
        # Diagnostic testing proved that ANY value for --disk-limit causes immediate failure
        # The OS already sets proper limits (995GB) which work perfectly
        # See diagnostic output: without --disk-limit = SUCCESS, with --disk-limit = FAIL
        if (cd "$srr_dir" && \
            fasterq-dump "$sra_input" '"$ngc_arg"' \
                --split-3 \
                --threads '"$threads_per_job"' \
                --temp "$job_temp_dir" \
                --mem '"$FASTERQ_MEM"' \
                --bufsize '"$FASTERQ_BUFSIZE"' \
                --curcache '"$FASTERQ_CACHE"' \
                --progress 2>> "'"$ERROR_LOG"'" && \
            echo "[$(date +%H:%M:%S)] Compressing FASTQ files..." && \
            pigz -1 -p '"$pigz_threads"' *.fastq 2>> "'"$ERROR_LOG"'" && \
            mv *.fastq.gz .. 2>> "'"$ERROR_LOG"'" && \
            cd .. && \
            rm -rf "${run}" "$job_temp_dir" 2>> "'"$ERROR_LOG"'"); then
            
            echo "[$(date +%H:%M:%S)] ✓ Complete"
            mark_successful "${sample}" "${run}"
        else
            echo "[$(date +%H:%M:%S)] ✗ Failed"
            log_step_error "EXTRACTION" "${sample}" "${run}" "Extraction or compression failed"
            mark_extraction_failed "${sample}" "${run}"
            # Clean up temp directory on failure
            rm -rf "$job_temp_dir" 2>/dev/null
        fi
        '
    
    rm -rf "$temp_dir"
    rm -f "$samples_to_process"
    
    # Count results
    local successful=0
    local failed_count=0
    
    [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
    [[ -f "$EXTRACTION_FAILED" ]] && failed_count=$(wc -l < "$EXTRACTION_FAILED")
    
    # SMART RETRY: If any samples failed, retry them one at a time
    # This gives speed for samples that work, reliability for problematic ones
    if [[ $failed_count -gt 0 ]]; then
        log_message ""
        log_warning "Initial parallel extraction had $failed_count failures"
        log_message "Retrying failed samples sequentially (one at a time) for reliability..."
        log_message ""
        
        # Create list of failed samples for retry
        local retry_list="$work_dir/samples_to_retry.txt"
        > "$retry_list"
        
        # Extract GSM and SRR from failed list (remove from EXTRACTION_FAILED)
        while read -r line; do
            local sample=$(echo "$line" | awk '{print $1}')
            local run=$(echo "$line" | awk '{print $2}')
            
            # Check if file still exists
            local sra_file="$work_dir/Raw_data/${sample}/${run}/${run}.sra"
            local sralite_file="$work_dir/Raw_data/${sample}/${run}/${run}.sralite"
            
            if [[ -f "$sra_file" ]] || [[ -f "$sralite_file" ]]; then
                echo -e "${sample}\t${run}" >> "$retry_list"
            fi
        done < "$EXTRACTION_FAILED"
        
        local retry_count=$(wc -l < "$retry_list")
        
        if [[ $retry_count -gt 0 ]]; then
            log_message "Retrying $retry_count samples with sequential processing..."
            
            # Clear the extraction failed list - we're retrying these
            > "$EXTRACTION_FAILED"
            
            # DYNAMIC: Recalculate memory for sequential mode (1 job = maximum memory!)
            local RETRY_MEM=$(calculate_fasterq_mem_for_jobs 1 $AVAILABLE_MEM_GB)
            local RETRY_BUFSIZE=$(calculate_fasterq_bufsize "$RETRY_MEM")
            local RETRY_CACHE=$(calculate_fasterq_cache "$RETRY_MEM")
            
            log_message "Sequential retry configuration:"
            log_message "  Jobs: 1 (sequential)"
            log_message "  Memory per job: $RETRY_MEM (maximum available for single job)"
            log_message "  Buffer: $RETRY_BUFSIZE, Cache: $RETRY_CACHE"
            log_message ""
            
            # Sequential extraction (parallel_jobs = 1)
            local retry_temp_dir
            if [[ "$RAM_DISK_AVAILABLE" == true ]]; then
                retry_temp_dir="/dev/shm/sra_retry_$$"
            else
                retry_temp_dir="$work_dir/tmp_retry"
            fi
            mkdir -p "$retry_temp_dir"
            
            # Process one at a time
            local retry_success=0
            local retry_fail=0
            
            while IFS=$'\t' read -r sample run; do
                local srr_dir="$work_dir/Raw_data/${sample}/${run}"
                
                log_message "[RETRY] Processing ${sample}/${run} (sequential mode)..."
                
                # Determine file format
                local sra_input=""
                if [[ -f "$srr_dir/${run}.sra" ]]; then
                    sra_input="${run}.sra"
                elif [[ -f "$srr_dir/${run}.sralite" ]]; then
                    sra_input="${run}.sralite"
                else
                    log_error "[RETRY] ${sample}/${run} - SRA file not found"
                    mark_extraction_failed "${sample}" "${run}"
                    ((retry_fail++))
                    continue
                fi
                
                # Create unique temp for this retry
                local retry_job_temp="${retry_temp_dir}/${sample}_${run}"
                mkdir -p "$retry_job_temp"
                
                # Extract with maximum available memory (sequential = 1 job)
                if (cd "$srr_dir" && \
                    fasterq-dump "$sra_input" $ngc_arg \
                        --split-3 \
                        --threads 2 \
                        --temp "$retry_job_temp" \
                        --mem "$RETRY_MEM" \
                        --bufsize "$RETRY_BUFSIZE" \
                        --curcache "$RETRY_CACHE" \
                        --progress 2>> "$ERROR_LOG" && \
                    echo "[$(date +%H:%M:%S)] Compressing FASTQ files..." && \
                    pigz -1 -p 2 *.fastq 2>> "$ERROR_LOG" && \
                    mv *.fastq.gz .. 2>> "$ERROR_LOG" && \
                    cd .. && \
                    rm -rf "${run}" "$retry_job_temp" 2>> "$ERROR_LOG"); then
                    
                    log_message "[RETRY] ✓ ${sample}/${run} succeeded"
                    mark_successful "${sample}" "${run}"
                    ((retry_success++))
                else
                    log_error "[RETRY] ✗ ${sample}/${run} failed again"
                    mark_extraction_failed "${sample}" "${run}"
                    rm -rf "$retry_job_temp" 2>/dev/null
                    ((retry_fail++))
                fi
            done < "$retry_list"
            
            # Clean up retry temp
            rm -rf "$retry_temp_dir"
            rm -f "$retry_list"
            
            log_message ""
            log_message "Sequential Retry Results:"
            log_message "  Succeeded on retry: $retry_success"
            log_message "  Still failed: $retry_fail"
            
            # PHASE 3: Ultimate fallback - retry persistent failures with DISK temp (995GB!)
            # If samples fail even with full RAM disk, they need more space than RAM provides
            if [[ $retry_fail -gt 0 ]]; then
                log_message ""
                log_warning "Phase 3: ${retry_fail} samples still failing - switching to DISK temp (995GB available)"
                log_message "This is slower but handles even the largest samples..."
                log_message ""
                
                # Create list of samples that failed sequential retry
                local disk_retry_list="$work_dir/samples_disk_retry.txt"
                > "$disk_retry_list"
                
                while read -r line; do
                    local sample=$(echo "$line" | awk '{print $1}')
                    local run=$(echo "$line" | awk '{print $2}')
                    
                    # Check if file still exists
                    local sra_file="$work_dir/Raw_data/${sample}/${run}/${run}.sra"
                    local sralite_file="$work_dir/Raw_data/${sample}/${run}/${run}.sralite"
                    
                    if [[ -f "$sra_file" ]] || [[ -f "$sralite_file" ]]; then
                        echo -e "${sample}\t${run}" >> "$disk_retry_list"
                    fi
                done < "$EXTRACTION_FAILED"
                
                local disk_retry_count=$(wc -l < "$disk_retry_list")
                
                if [[ $disk_retry_count -gt 0 ]]; then
                    log_message "Retrying $disk_retry_count samples with DISK temp (/data)..."
                    
                    # Clear extraction failed - we're retrying these one last time
                    > "$EXTRACTION_FAILED"
                    
                    # DYNAMIC: Recalculate memory for disk retry (1 job = maximum memory!)
                    local DISK_RETRY_MEM=$(calculate_fasterq_mem_for_jobs 1 $AVAILABLE_MEM_GB)
                    local DISK_RETRY_BUFSIZE=$(calculate_fasterq_bufsize "$DISK_RETRY_MEM")
                    local DISK_RETRY_CACHE=$(calculate_fasterq_cache "$DISK_RETRY_MEM")
                    
                    log_message "Disk retry configuration:"
                    log_message "  Jobs: 1 (sequential)"
                    log_message "  Memory per job: $DISK_RETRY_MEM (maximum available)"
                    log_message "  Temp: /data disk (995GB available)"
                    log_message ""
                    
                    # Use DISK temp directory (not RAM)
                    local disk_temp_dir="$work_dir/tmp_disk_retry"
                    mkdir -p "$disk_temp_dir"
                    
                    local disk_success=0
                    local disk_fail=0
                    
                    while IFS=$'\t' read -r sample run; do
                        local srr_dir="$work_dir/Raw_data/${sample}/${run}"
                        
                        log_message "[DISK-RETRY] Processing ${sample}/${run} with disk temp..."
                        
                        # Determine file format
                        local sra_input=""
                        if [[ -f "$srr_dir/${run}.sra" ]]; then
                            sra_input="${run}.sra"
                        elif [[ -f "$srr_dir/${run}.sralite" ]]; then
                            sra_input="${run}.sralite"
                        else
                            log_error "[DISK-RETRY] ${sample}/${run} - SRA file not found"
                            mark_extraction_failed "${sample}" "${run}"
                            ((disk_fail++))
                            continue
                        fi
                        
                        # Create unique temp on DISK (not RAM)
                        local disk_job_temp="${disk_temp_dir}/${sample}_${run}"
                        mkdir -p "$disk_job_temp"
                        
                        # Extract using DISK temp (995GB available!)
                        # With maximum memory per job
                        if (cd "$srr_dir" && \
                            fasterq-dump "$sra_input" $ngc_arg \
                                --split-3 \
                                --threads 2 \
                                --temp "$disk_job_temp" \
                                --mem "$DISK_RETRY_MEM" \
                                --bufsize "$DISK_RETRY_BUFSIZE" \
                                --curcache "$DISK_RETRY_CACHE" \
                                --progress 2>> "$ERROR_LOG" && \
                            echo "[$(date +%H:%M:%S)] Compressing FASTQ files..." && \
                            pigz -1 -p 2 *.fastq 2>> "$ERROR_LOG" && \
                            mv *.fastq.gz .. 2>> "$ERROR_LOG" && \
                            cd .. && \
                            rm -rf "${run}" "$disk_job_temp" 2>> "$ERROR_LOG"); then
                            
                            log_message "[DISK-RETRY] ✓ ${sample}/${run} succeeded with disk temp!"
                            mark_successful "${sample}" "${run}"
                            ((disk_success++))
                        else
                            log_error "[DISK-RETRY] ✗ ${sample}/${run} failed even with disk temp"
                            log_error "This sample may be corrupted or have other issues"
                            mark_extraction_failed "${sample}" "${run}"
                            rm -rf "$disk_job_temp" 2>/dev/null
                            ((disk_fail++))
                        fi
                    done < "$disk_retry_list"
                    
                    # Clean up disk retry temp
                    rm -rf "$disk_temp_dir"
                    rm -f "$disk_retry_list"
                    
                    log_message ""
                    log_message "Disk Retry Results (Phase 3):"
                    log_message "  Succeeded with disk temp: $disk_success"
                    log_message "  Still failed (likely corrupted): $disk_fail"
                fi
            fi
            
            # Update final counts
            successful=0
            failed_count=0
            [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
            [[ -f "$EXTRACTION_FAILED" ]] && failed_count=$(wc -l < "$EXTRACTION_FAILED")
        fi
    fi
    
    log_message ""
    log_message "Extraction Summary:"
    log_message "  Extracted successfully: $successful"
    log_message "  Failed extraction: $failed_count"
    log_message "  Skipped (previous errors): $skipped"
    
    if [[ $failed_count -gt 0 ]]; then
        log_warning "Some files failed extraction - see $EXTRACTION_FAILED"
    fi
}

################################################################################
# Cleanup Handler
################################################################################

cleanup_on_exit() {
    local exit_code=$?
    
    # Clean up any temporary files
    rm -f samples_to_prefetch.txt 2>/dev/null
    rm -f samples_to_validate.txt 2>/dev/null
    rm -f samples_to_extract.txt 2>/dev/null
    
    # Clean up RAM disk temp directory if it exists
    if [[ -d "/dev/shm/sra_extraction_$$" ]]; then
        rm -rf "/dev/shm/sra_extraction_$$" 2>/dev/null
    fi
    
    # Only log if we're in an actual run (not during argument parsing)
    if [[ -n "${OUTPUT_DIR:-}" ]] && [[ -d "${BASE_PATH}/${OUTPUT_DIR}" ]]; then
        if [[ $exit_code -ne 0 ]]; then
            log_error "Script exited with errors (code: $exit_code)"
            if [[ -n "${ERROR_LOG:-}" ]]; then
                log_error "Check $ERROR_LOG for details"
            fi
        fi
    fi
}

trap cleanup_on_exit EXIT INT TERM

################################################################################
# Cleanup and Summary
################################################################################

cleanup_temp_files() {
    log_message "Cleaning up temporary files..."
    
    rm -f pre_tmp_gsm_to_srr.txt
    rm -f tmp_gse_to_gsm_full.txt
    rm -f tmp_gsm_to_srr_full.txt
    
    log_success "Cleanup completed"
}

generate_summary() {
    log_message "Generating download summary..."
    
    local sample_count=$(find Raw_data -maxdepth 1 -type d 2>/dev/null | tail -n +2 | wc -l)
    local fastq_count=$(find Raw_data -name "*.fastq.gz" 2>/dev/null | wc -l)
    local total_size=$(du -sh Raw_data 2>/dev/null | cut -f1)
    
    # Count errors by stage
    local prefetch_errors=0
    local validation_errors=0
    local extraction_errors=0
    local successful=0
    
    [[ -f "$PREFETCH_FAILED" ]] && prefetch_errors=$(wc -l < "$PREFETCH_FAILED")
    [[ -f "$VALIDATION_FAILED" ]] && validation_errors=$(wc -l < "$VALIDATION_FAILED")
    [[ -f "$EXTRACTION_FAILED" ]] && extraction_errors=$(wc -l < "$EXTRACTION_FAILED")
    [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
    
    cat > download_summary.txt << EOF
Download Summary
================================================================================
Accession: $ACCESSION ($ACCESSION_TYPE)
Output Directory: $OUTPUT_DIR
Date: $(date)
GSM Filter: ${GSM_SAMPLES}
Mode: $(if [[ "$MANUAL_MODE" == true ]]; then echo "Manual"; else echo "Automatic"; fi)
NGC File: ${NGC_FILE:-None}
System CPUs: $TOTAL_CPUS
CPUs Used: $CPUS

Results:
--------
Samples Downloaded: $sample_count
FASTQ Files: $fastq_count
Total Size: ${total_size:-Unknown}

Processing Status:
------------------
Successful (complete pipeline): $successful
Failed at Prefetch: $prefetch_errors
Failed at Validation: $validation_errors
Failed at Extraction: $extraction_errors

Error Logs:
-----------
Main error log: $ERROR_LOG
Prefetch failures: $PREFETCH_FAILED
Validation failures: $VALIDATION_FAILED
Extraction failures: $EXTRACTION_FAILED

Successfully Processed Files:
-----------------------------
EOF
    
    if [[ -f "$SUCCESSFUL_SAMPLES" ]]; then
        cat "$SUCCESSFUL_SAMPLES" >> download_summary.txt
    fi
    
    echo "" >> download_summary.txt
    echo "FASTQ Files:" >> download_summary.txt
    echo "------------" >> download_summary.txt
    find Raw_data -name "*.fastq.gz" | sort >> download_summary.txt
    
    log_success "Summary saved to download_summary.txt"
    
    echo ""
    log_message "========================================="
    log_message "DOWNLOAD COMPLETED"
    log_message "========================================="
    log_message "Accession: $ACCESSION"
    log_message "Samples: $sample_count"
    log_message "FASTQ files: $fastq_count"
    log_message "Total size: ${total_size:-Unknown}"
    log_message ""
    log_message "Processing Summary:"
    log_message "  Successful: $successful"
    log_message "  Failed at prefetch: $prefetch_errors"
    log_message "  Failed at validation: $validation_errors"
    log_message "  Failed at extraction: $extraction_errors"
    log_message ""
    log_message "Location: $(pwd)"
    log_message "Error log: $ERROR_LOG"
    log_message "========================================="
    
    if [[ $((prefetch_errors + validation_errors + extraction_errors)) -gt 0 ]]; then
        log_warning "Some files failed processing - check error logs for details"
    fi
}

################################################################################
# Single Dataset Processing
################################################################################

process_single_dataset() {
    log_message "========================================="
    log_message "Processing Single Dataset"
    log_message "========================================="
    log_message "Output Directory: $OUTPUT_DIR"
    log_message "Accession: $ACCESSION ($ACCESSION_TYPE)"
    log_message "GSM Samples: $GSM_SAMPLES"
    log_message "NGC File: ${NGC_FILE:-None}"
    log_message "========================================="
    
    detect_system_resources
    check_dependencies
    create_output_structure
    
    if [[ "$ACCESSION_TYPE" == "gse" ]]; then
        log_message "Processing GSE accession..."
        
        get_sample_accessions_gse || exit 1
        get_run_accessions_gse || exit 1
        
        # Check disk space before starting
        local num_samples=$(wc -l < tmp_gsm_to_srr.txt)
        check_disk_space $num_samples
        
        prefetch_data 1 2 tmp_gsm_to_srr.txt
        validate_data 1 2 tmp_gsm_to_srr.txt
        extract_and_compress 1 2 tmp_gsm_to_srr.txt
        
    elif [[ "$ACCESSION_TYPE" == "srp" ]]; then
        log_message "Processing SRP accession..."
        
        get_sample_accessions_srp
        
        # Check disk space before starting
        local num_samples=$(wc -l < srx_to_srr.tsv)
        check_disk_space $num_samples
        
        prefetch_data 1 2 srx_to_srr.tsv
        validate_data 1 2 srx_to_srr.tsv
        extract_and_compress 1 2 srx_to_srr.tsv
    fi
    
    cleanup_temp_files
    generate_summary
    
    log_success "Dataset processing completed!"
}

################################################################################
# Batch Processing from Input File
################################################################################

process_input_file() {
    log_message "========================================="
    log_message "Processing Input File: $INPUT_FILE"
    log_message "========================================="
    
    if [[ ! -f "$INPUT_FILE" ]]; then
        log_error "Input file not found: $INPUT_FILE"
        exit 1
    fi
    
    detect_system_resources
    
    local total_lines=$(grep -v "^#" "$INPUT_FILE" | grep -v "^$" | wc -l)
    log_message "Found $total_lines lines in input file"
    log_message "Processing $PARALLEL_DATASETS dataset(s) at a time"
    log_message "========================================="
    
    if [[ $PARALLEL_DATASETS -eq 1 ]]; then
        process_sequential
    else
        process_parallel
    fi
}

process_sequential() {
    local dataset_num=0
    local success_count=0
    local fail_count=0
    
    log_message "Aggregating GSM samples by dataset..."
    
    # Use awk to aggregate - much more reliable than bash associative arrays
    local aggregated=$(mktemp)
    
    awk -F'\t' '
        # Skip comments and empty lines
        /^#/ { next }
        /^[[:space:]]*$/ { next }
        
        {
            # Clean up whitespace
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $3)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $4)
            
            # Skip if key fields are empty
            if ($1 == "" || $2 == "" || $3 == "") next
            
            # Create key using | as delimiter (simpler, less likely to cause regex issues)
            key = $1 "|" $2
            
            # Aggregate GSMs
            if (key in gsms) {
                gsms[key] = gsms[key] "," $3
            } else {
                gsms[key] = $3
                order[++n] = key
            }
            
            # Store NGC file
            if ($4 != "") {
                ngc[key] = $4
            }
        }
        
        END {
            for (i = 1; i <= n; i++) {
                key = order[i]
                # Use sub() to split at first |
                output_dir = key
                sub(/\|.*/, "", output_dir)
                accession = key
                sub(/[^|]*\|/, "", accession)
                
                gsm_list = gsms[key]
                ngc_file = (key in ngc) ? ngc[key] : ""
                
                printf "%s\t%s\t%s\t%s\n", output_dir, accession, gsm_list, ngc_file
            }
        }
    ' "$INPUT_FILE" > "$aggregated"
    
    local total_datasets=$(wc -l < "$aggregated")
    log_message "Aggregated into $total_datasets unique dataset(s)"
    log_message ""
    
    # Temporarily disable set -e to prevent issues with read returning 1 at EOF
    set +e
    
    # Process each aggregated dataset
    local line_num=0
    
    while IFS=$'\t' read -r output_dir accession gsm_list ngc_file || [[ -n "$output_dir" ]]; do
        ((line_num++))
        
        # Skip empty lines
        if [[ -z "$output_dir" ]]; then
            continue
        fi
        if [[ -z "$accession" ]]; then
            continue
        fi
        
        ((dataset_num++))
        
        log_message ""
        log_message "========================================="
        log_message "Dataset $dataset_num of $total_datasets: $output_dir"
        log_message "========================================="
        log_message "  Accession: $accession"
        log_message "  GSMs: ${gsm_list:0:100}..."
        log_message "  NGC: ${ngc_file:-none}"
        log_message ""
        
        # Instead of recursive call, directly process the dataset
        OUTPUT_DIR="$output_dir"
        ACCESSION="$accession"
        GSM_SAMPLES="$gsm_list"
        NGC_FILE="$ngc_file"
        
        # Auto-detect accession type
        if [[ $ACCESSION =~ ^GSE ]]; then
            ACCESSION_TYPE="gse"
        elif [[ $ACCESSION =~ ^SRP ]]; then
            ACCESSION_TYPE="srp"
        else
            log_error "Unknown accession type: $ACCESSION (expected GSE* or SRP*)"
            ((fail_count++))
            continue
        fi
        
        # Process this dataset directly
        if process_single_dataset; then
            ((success_count++))
            log_success "Dataset $dataset_num completed: $output_dir"
        else
            ((fail_count++))
            log_error "Dataset $dataset_num failed: $output_dir"
        fi
        
    done < "$aggregated"
    
    # Re-enable set -e
    set -e
    
    rm -f "$aggregated"
    
    log_message ""
    log_message "========================================="
    log_message "BATCH PROCESSING COMPLETE"
    log_message "========================================="
    log_message "Total datasets: $dataset_num"
    log_message "Successful: $success_count"
    log_message "Failed: $fail_count"
    log_message "========================================="
}

process_parallel() {
    log_message "Starting parallel processing with $PARALLEL_DATASETS jobs..."
    
    # First, aggregate GSMs by Output_Directory+Accession
    local temp_aggregated=$(mktemp)
    
    declare -A gsm_lists
    declare -A ngc_files
    
    while IFS=$'\t' read -r output_dir accession gsm_sample ngc_file || [[ -n "$output_dir" ]]; do
        [[ "$output_dir" =~ ^#.*$ ]] && continue
        [[ -z "$output_dir" ]] && continue
        [[ -z "$accession" ]] && continue
        [[ -z "$gsm_sample" ]] && continue
        
        local key="${output_dir}^^^${accession}"
        
        # Aggregate GSM samples
        if [[ -z "${gsm_lists[$key]:-}" ]]; then
            gsm_lists[$key]="${gsm_sample}"
        else
            gsm_lists[$key]="${gsm_lists[$key]},${gsm_sample}"
        fi
        
        # Store NGC file if provided
        if [[ -n "${ngc_file:-}" ]]; then
            ngc_files[$key]="${ngc_file}"
        fi
        
    done < <(grep -v "^#" "$INPUT_FILE" | grep -v "^$")
    
    # Write aggregated data to temp file
    for key in "${!gsm_lists[@]}"; do
        local output_dir="${key%^^^*}"
        local accession="${key#*^^^}"
        local gsm_list="${gsm_lists[$key]}"
        local ngc_file="${ngc_files[$key]:-}"
        
        printf "%s\t%s\t%s\t%s\n" "$output_dir" "$accession" "$gsm_list" "$ngc_file" >> "$temp_aggregated"
    done
    
    local total=$(wc -l < "$temp_aggregated")
    log_message "Aggregated into $total unique dataset(s)"
    
    # Export functions for parallel
    export -f log_message log_error log_warning log_success
    export BASE_PATH PARALLEL_DATASETS
    
    # Process using GNU parallel
    cat "$temp_aggregated" | \
        parallel --colsep '\t' -j $PARALLEL_DATASETS --line-buffer \
        bash "$0" {1} {2} {3} {4}
    
    rm -f "$temp_aggregated"
    
    log_message ""
    log_message "========================================="
    log_success "PARALLEL PROCESSING COMPLETE"
    log_message "========================================="
}

################################################################################
# Main Entry Point
################################################################################

main() {
    parse_arguments "$@"
    
    if [[ -n "$INPUT_FILE" ]]; then
        process_input_file
    else
        process_single_dataset
    fi
}

main "$@"
