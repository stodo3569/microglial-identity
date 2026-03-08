#!/bin/bash

################################################################################
# Salmon Quantification Script - Docker Version
#
# Mirrors the architecture of geo_download_docker.sh and fastp_trim_docker.sh
# for robust, parallel salmon quantification of trimmed FASTQ files.
#
# Features:
#   - RAM-aware adaptive parallelism (memory per job calculation)
#   - 3-phase retry strategy:
#       Phase 1: Parallel with calculated resources
#       Phase 2: Sequential with maximum resources (all CPUs + all RAM)
#       Phase 3: Minimal-memory fallback (reduced threads, no bias models)
#   - Comprehensive error logging and per-stage tracking
#   - Automatic paired-end / single-end detection from fastp output
#   - Correct handling of technical replicates (multiple SRR per sample)
#   - Batch processing of multiple studies from input file
#   - Skip already-quantified samples (unless --force)
#
# Usage:
#   # Quantify all samples in a single study
#   bash salmon_quant_docker.sh --index /data/salmon_index "Study_Name"
#
#   # Quantify specific samples
#   bash salmon_quant_docker.sh --index /data/salmon_index "Study_Name" "GSM123,GSM456"
#
#   # Batch from input file (same format as geo_download_docker.sh)
#   bash salmon_quant_docker.sh --index /data/salmon_index --input-file datasets.txt
#
# Docker example:
#   docker run --rm \
#     -u "$(id -u):$(id -g)" \
#     --shm-size=80g \
#     -e SALMON_SHOW_STDERR=true \
#     -e SALMON_HEARTBEAT_SEC=30 \
#     -v /data:/data \
#     stodo3569/salmon-tools:0.0 \
#     bash /data/geo_scripts/salmon_quant_docker.sh \
#       --index /data/salmon_index \
#       --input-file /data/datasets.txt
#
# Input file format (same as geo_download_docker.sh / fastp_trim_docker.sh):
#   Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC_file]
#
# Expected directory structure (created by fastp_trim_docker.sh):
#   /data/Study_Name/Trimmed_data/GSM123456/fastp_*.fastq.gz
#
# Technical replicate handling:
#   fastp produces files like:
#     Paired:  fastp_SRR111_1.fastq.gz  fastp_SRR111_2.fastq.gz
#              fastp_SRR222_1.fastq.gz  fastp_SRR222_2.fastq.gz
#     Single:  fastp_SRR111.fastq.gz    fastp_SRR222.fastq.gz
#
#   Salmon receives ALL R1 files via -1 and ALL R2 files via -2 (space-separated),
#   or ALL single-end files via -r. This correctly merges technical replicates.
#
# Output structure:
#   /data/Study_Name/Aligned_data/GSM123456/quant.sf
#   /data/Study_Name/Aligned_data/GSM123456/quant.genes.sf
#   /data/Study_Name/Aligned_data/GSM123456/aux_info/
#   /data/Study_Name/Aligned_data/GSM123456/cmd_info.json
#   /data/Study_Name/Aligned_data/GSM123456/lib_format_counts.json
#   /data/Study_Name/Aligned_data/GSM123456/logs/
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

    # Adjust for parallel datasets to prevent CPU oversubscription
    local available_cpus=$TOTAL_CPUS
    if [[ -n "${PARALLEL_DATASETS:-}" ]] && [[ $PARALLEL_DATASETS -gt 1 ]]; then
        available_cpus=$((TOTAL_CPUS / PARALLEL_DATASETS))
        log_message "CPU allocation adjusted for parallel processing:"
        log_message "  Total system CPUs: $TOTAL_CPUS"
        log_message "  Parallel datasets: $PARALLEL_DATASETS"
        log_message "  CPUs per dataset: ~${available_cpus}"
    fi

    # Smart CPU reservation
    if [[ $available_cpus -gt 8 ]]; then
        CPUS=$((available_cpus - 2))
    elif [[ $available_cpus -gt 4 ]]; then
        CPUS=$((available_cpus - 1))
    else
        CPUS=$available_cpus
    fi

    if [[ $CPUS -lt 2 ]]; then
        CPUS=2
        log_warning "Very low CPU allocation detected - setting minimum of 2 CPUs"
    fi

    # Detect total RAM
    if [[ -f /proc/meminfo ]]; then
        TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        TOTAL_MEM_GB=$((TOTAL_MEM_KB / 1024 / 1024))

        AVAILABLE_MEM_KB=$(grep MemAvailable /proc/meminfo 2>/dev/null | awk '{print $2}')
        if [[ -n "$AVAILABLE_MEM_KB" ]]; then
            AVAILABLE_MEM_GB=$((AVAILABLE_MEM_KB / 1024 / 1024))
        else
            # Fallback: estimate 70% of total is available
            AVAILABLE_MEM_GB=$((TOTAL_MEM_GB * 70 / 100))
        fi

        if [[ -n "${PARALLEL_DATASETS:-}" ]] && [[ $PARALLEL_DATASETS -gt 1 ]]; then
            AVAILABLE_MEM_GB=$((AVAILABLE_MEM_GB / PARALLEL_DATASETS))
            log_message "  Memory per dataset: ~${AVAILABLE_MEM_GB}GB"
        fi
    else
        TOTAL_MEM_GB=8
        AVAILABLE_MEM_GB=6
    fi

    # Estimate salmon index size for memory planning
    estimate_index_memory

    log_message "System Resources Detected:"
    log_message "  Total CPUs: $TOTAL_CPUS"
    log_message "  CPUs for this dataset: $CPUS"
    log_message "  Total Memory: ${TOTAL_MEM_GB}GB"
    log_message "  Available Memory: ${AVAILABLE_MEM_GB}GB"
    log_message "  Estimated index memory: ~${INDEX_MEM_GB}GB per salmon job"
}

################################################################################
# Index Memory Estimation
################################################################################

estimate_index_memory() {
    # Salmon loads the index into memory for each process.
    # A typical human transcriptome index is ~4-8GB in RAM.
    # We estimate from the on-disk size (in-memory is typically 1.5-3x disk size).

    INDEX_MEM_GB=6  # Conservative default

    if [[ -n "$SALMON_INDEX" ]] && [[ -d "$SALMON_INDEX" ]]; then
        local index_disk_mb=$(du -sm "$SALMON_INDEX" 2>/dev/null | cut -f1)
        if [[ -n "$index_disk_mb" ]] && [[ "$index_disk_mb" =~ ^[0-9]+$ ]]; then
            # In-memory footprint is roughly 2x the on-disk size
            local estimated_mem_mb=$((index_disk_mb * 2))
            INDEX_MEM_GB=$(( (estimated_mem_mb / 1024) + 1 ))  # Round up

            # Apply sensible bounds
            [[ $INDEX_MEM_GB -lt 2 ]] && INDEX_MEM_GB=2
            [[ $INDEX_MEM_GB -gt 32 ]] && INDEX_MEM_GB=32

            log_message "  Index on-disk: ~${index_disk_mb}MB -> estimated ~${INDEX_MEM_GB}GB in RAM"
        fi
    fi
}

################################################################################
# Adaptive Resource Calculation
################################################################################

calculate_salmon_threads() {
    # Salmon thread allocation per job
    # Salmon scales well up to ~8-12 threads; diminishing returns beyond that
    if [[ $CPUS -ge 32 ]]; then
        echo 8
    elif [[ $CPUS -ge 16 ]]; then
        echo 6
    elif [[ $CPUS -ge 8 ]]; then
        echo 4
    elif [[ $CPUS -ge 4 ]]; then
        echo 4
    else
        echo 2
    fi
}

calculate_optimal_parallelism() {
    local num_samples=${1:-10}

    # Salmon is BOTH CPU- and memory-bound.
    # Each salmon job loads the index into memory (~INDEX_MEM_GB) plus
    # working memory that scales with threads (~0.5-1GB per thread).
    #
    # Strategy: min(cpu_constraint, memory_constraint)

    # CPU constraint: total CPUs / threads per job
    local cpu_jobs=$((CPUS / SALMON_THREADS))
    [[ $cpu_jobs -lt 1 ]] && cpu_jobs=1

    # Memory constraint: each job needs index + working memory
    local mem_per_job=$(calculate_required_mem_per_job "$SALMON_THREADS")

    # Reserve 10% for system overhead (matching geo_download_docker.sh pattern)
    local usable_mem=$((AVAILABLE_MEM_GB * 90 / 100))
    local mem_jobs=1
    if [[ $usable_mem -gt 0 ]] && [[ $mem_per_job -gt 0 ]]; then
        mem_jobs=$((usable_mem / mem_per_job))
        [[ $mem_jobs -lt 1 ]] && mem_jobs=1
    fi

    local optimal_jobs=$cpu_jobs
    local limiting_factor="CPU"

    if [[ $mem_jobs -lt $optimal_jobs ]]; then
        optimal_jobs=$mem_jobs
        limiting_factor="RAM"
    fi

    # Keep function return value clean for command substitution by logging to stderr
    log_message "  Parallelism calculation:" >&2
    log_message "    CPU constraint: $cpu_jobs jobs ($CPUS CPUs / $SALMON_THREADS threads)" >&2
    log_message "    RAM constraint: $mem_jobs jobs (${usable_mem}GB usable / ~${mem_per_job}GB per job)" >&2
    log_message "    -> $optimal_jobs parallel jobs (limited by $limiting_factor)" >&2

    # Don't spawn more jobs than samples
    if [[ $num_samples -lt $optimal_jobs ]]; then
        echo $num_samples
    else
        echo $optimal_jobs
    fi
}

calculate_required_mem_per_job() {
    local threads=${1:-$SALMON_THREADS}

    # Working memory ~0.5GB per thread + 1GB overhead
    local working_mem_per_job=$(( (threads / 2) + 1 ))
    echo $((INDEX_MEM_GB + working_mem_per_job))
}

calculate_mem_per_job() {
    local parallel_jobs=$1
    local available_mem=$2  # in GB

    # Calculate memory per job dynamically (matching geo_download_docker.sh pattern)
    # Reserve 10% for system overhead
    local usable_mem=$((available_mem * 90 / 100))
    local mem_per_job=$((usable_mem / parallel_jobs))

    # Apply sensible limits for salmon
    if [[ $mem_per_job -gt 64 ]]; then
        mem_per_job=64  # Cap (salmon doesn't benefit from more)
    elif [[ $mem_per_job -lt 4 ]]; then
        mem_per_job=4   # Minimum: index + minimal working set
    fi

    echo "$mem_per_job"
}

################################################################################
# Configuration
################################################################################

BASE_PATH="${GEO_BASE_PATH:-/data}"
PARALLEL_DATASETS=1

STUDY_DIR=""
SAMPLE_FILTER=""
INPUT_FILE=""
FORCE_RERUN=false

CPUS=4
TOTAL_CPUS=4
TOTAL_MEM_GB=8
AVAILABLE_MEM_GB=6
INDEX_MEM_GB=6

# Salmon parameters (from salmon_quant_10.sh)
SALMON_INDEX=""
SALMON_THREADS=4
SALMON_LIBTYPE_PAIRED="A"
SALMON_LIBTYPE_SINGLE="A"
SALMON_EXTRA_ARGS="--validateMappings --seqBias --gcBias --posBias --dumpEq"

# Error tracking files (set per-study in initialize_error_tracking)
ERROR_LOG=""
QUANT_FAILED=""
SUCCESSFUL_SAMPLES=""
SKIPPED_SAMPLES=""
RAM_USAGE_LOG=""

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
    local error_msg="${3:-Unknown error}"

    local msg="[$step] Failed: $sample - $error_msg"
    log_error "$msg"
}

################################################################################
# Error Tracking Functions
################################################################################

initialize_error_tracking() {
    local study_path="$1"

    ERROR_LOG="${study_path}/salmon_error_log.txt"
    QUANT_FAILED="${study_path}/salmon_quant_failed.txt"
    SUCCESSFUL_SAMPLES="${study_path}/salmon_successful_samples.txt"
    SKIPPED_SAMPLES="${study_path}/salmon_skipped_samples.txt"
    RAM_USAGE_LOG="${study_path}/salmon_ram_usage.tsv"

    > "$ERROR_LOG"
    > "$QUANT_FAILED"
    > "$SUCCESSFUL_SAMPLES"
    > "$SKIPPED_SAMPLES"
    > "$RAM_USAGE_LOG"
    echo -e "sample\tphase\tthreads\tstatus\tpeak_kb\tpeak_gb" >> "$RAM_USAGE_LOG"

    log_message "Error tracking initialized"
    log_message "  Error log: $ERROR_LOG"
    log_message "  RAM usage log: $RAM_USAGE_LOG"
    log_message "  Failed samples will be tracked"
}

mark_quant_failed() {
    local sample="$1"
    echo "${sample}" >> "$QUANT_FAILED"
}

mark_successful() {
    local sample="$1"
    echo "${sample}" >> "$SUCCESSFUL_SAMPLES"
}

mark_skipped() {
    local sample="$1"
    echo "${sample}" >> "$SKIPPED_SAMPLES"
}

is_sample_failed() {
    local sample="$1"

    if [[ -f "$QUANT_FAILED" ]] && grep -q "^${sample}$" "$QUANT_FAILED"; then
        return 0
    fi

    return 1
}

################################################################################
# Argument Parsing
################################################################################

show_help() {
    cat << 'EOF'
Salmon Quantification Script - Docker Version

USAGE:
  Single study:
    $0 --index INDEX_PATH STUDY_DIR [SAMPLE_FILTER]

  Multiple studies from file:
    $0 --index INDEX_PATH --input-file FILE [--parallel N]

ARGUMENTS:
  STUDY_DIR       Name of the study directory under BASE_PATH (e.g., "Friedman_2019")
  SAMPLE_FILTER   Optional comma-separated GSM/SRX IDs to process (default: all)

REQUIRED OPTIONS:
  --index PATH        Path to the salmon index directory (REQUIRED)

OPTIONS:
  --input-file FILE   Process multiple studies from a tab-delimited file
                      (same format as geo_download_docker.sh)
  --parallel N        Process N studies simultaneously (default: 1)
  --force             Re-quantify samples even if output already exists
  --threads N         Override salmon threads per sample (default: auto-detected)
  --libtype-pe TYPE   Override paired-end library type (default: A = auto-detect)
  --libtype-se TYPE   Override single-end library type (default: A = auto-detect)
  --extra-args "..."  Override extra salmon arguments
                      (default: "--validateMappings --seqBias --gcBias --posBias --dumpEq")
  --base-path PATH    Override base data path (default: /data or $GEO_BASE_PATH)

EXAMPLES:
  # Quantify all samples in a study
  bash salmon_quant_docker.sh --index /data/salmon_index "Friedman_2019"

  # Quantify specific samples
  bash salmon_quant_docker.sh --index /data/salmon_index "Friedman_2019" "GSM3559136,GSM3559137"

  # Batch from input file
  bash salmon_quant_docker.sh --index /data/salmon_index --input-file /data/datasets.txt

  # Parallel studies
  bash salmon_quant_docker.sh --index /data/salmon_index --input-file /data/datasets.txt --parallel 2

DOCKER:
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    --shm-size=80g \
    -e SALMON_SHOW_STDERR=true \
    -e SALMON_HEARTBEAT_SEC=30 \
    -v /data:/data \
    stodo3569/salmon-tools:0.0 \
    bash /data/geo_scripts/salmon_quant_docker.sh \
      --index /data/salmon_index \
      --input-file /data/datasets.txt

RUNTIME ENV VARS:
  SALMON_SHOW_STDERR=true   Stream salmon stderr to terminal while also logging to error file
  SALMON_HEARTBEAT_SEC=30   Print "salmon still running..." heartbeat every N seconds

SALMON PARAMETERS (from salmon_quant_10.sh):
  --validateMappings     Enable selective alignment validation
  --seqBias              Correct for sequence-specific bias
  --gcBias               Correct for fragment GC bias
  --posBias              Correct for positional bias
  --dumpEq               Dump equivalence classes

RETRY STRATEGY (3-phase, from geo_download_docker.sh):
  Phase 1: Parallel quantification with adaptive resources
  Phase 2: Sequential retry with maximum resources (all CPUs + all RAM for 1 job)
  Phase 3: Minimal-memory fallback (2 threads, bias models disabled)

TECHNICAL REPLICATES:
  If a sample (e.g., GSM123456) has multiple runs (SRR111, SRR222), fastp
  produces separate trimmed files per run. This script collects ALL R1 files
  and ALL R2 files within a sample directory and passes them to salmon as
  space-separated lists, correctly merging technical replicates.

OUTPUT STRUCTURE:
  /data/Study_Name/
    Trimmed_data/GSM123456/fastp_*.fastq.gz  (input - unchanged)
    Aligned_data/GSM123456/quant.sf           (transcript quantification)
    Aligned_data/GSM123456/aux_info/          (auxiliary info)
    Aligned_data/GSM123456/cmd_info.json      (command used)
    Aligned_data/GSM123456/logs/              (salmon logs)

ERROR TRACKING:
  All errors are logged to salmon_error_log.txt per study
  Failed samples tracked in salmon_quant_failed.txt
  Processing continues for other samples when some fail
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

    local positional_args=()
    local found_input_file=false

    local args=("$@")
    local i=0
    while [[ $i -lt ${#args[@]} ]]; do
        case "${args[$i]}" in
            --input-file)
                found_input_file=true
                INPUT_FILE="${args[$((i+1))]}"
                i=$((i + 2)) ;;
            --index)
                SALMON_INDEX="${args[$((i+1))]}"
                i=$((i + 2)) ;;
            --parallel)
                PARALLEL_DATASETS="${args[$((i+1))]}"
                i=$((i + 2)) ;;
            --force)
                FORCE_RERUN=true
                i=$((i + 1)) ;;
            --threads)
                SALMON_THREADS="${args[$((i+1))]}"
                i=$((i + 2)) ;;
            --libtype-pe)
                SALMON_LIBTYPE_PAIRED="${args[$((i+1))]}"
                i=$((i + 2)) ;;
            --libtype-se)
                SALMON_LIBTYPE_SINGLE="${args[$((i+1))]}"
                i=$((i + 2)) ;;
            --extra-args)
                SALMON_EXTRA_ARGS="${args[$((i+1))]}"
                i=$((i + 2)) ;;
            --base-path)
                BASE_PATH="${args[$((i+1))]}"
                i=$((i + 2)) ;;
            --*)
                log_error "Unknown option: ${args[$i]}"
                exit 1 ;;
            *)
                positional_args+=("${args[$i]}")
                i=$((i + 1)) ;;
        esac
    done

    # Validate required --index
    if [[ -z "$SALMON_INDEX" ]]; then
        log_error "Missing required option: --index PATH"
        log_error "A salmon index is required for quantification."
        log_error ""
        log_error "To build an index:"
        log_error "  1. Create transcriptome FASTA with gffread:"
        log_error "     gffread -w transcriptome.fa -g genome.fa annotations.gtf"
        log_error "  2. Build salmon index:"
        log_error "     salmon index -t transcriptome.fa -i salmon_index -k 31"
        exit 1
    fi

    if [[ "$found_input_file" == true ]]; then
        return 0
    fi

    # Single study mode
    if [[ ${#positional_args[@]} -lt 1 ]]; then
        log_error "Missing required argument: STUDY_DIR"
        log_error "Usage: $0 --index INDEX_PATH STUDY_DIR [SAMPLE_FILTER] [OPTIONS]"
        exit 1
    fi

    STUDY_DIR="${positional_args[0]}"
    SAMPLE_FILTER="${positional_args[1]:-}"
}

################################################################################
# Dependency Checking
################################################################################

check_dependencies() {
    log_message "Checking dependencies..."

    local deps=(salmon parallel)
    local missing=()

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

    # Validate salmon index
    if [[ ! -d "$SALMON_INDEX" ]]; then
        log_error "Salmon index directory not found: $SALMON_INDEX"
        exit 1
    fi

    if [[ ! -f "$SALMON_INDEX/versionInfo.json" ]] && \
       [[ ! -f "$SALMON_INDEX/info.json" ]] && \
       [[ ! -f "$SALMON_INDEX/duplicate_clusters.tsv" ]]; then
        log_warning "Salmon index may be incomplete - expected index files not found in $SALMON_INDEX"
        log_warning "Proceeding anyway - salmon will report errors if the index is invalid"
    else
        log_success "Salmon index validated: $SALMON_INDEX"
    fi

    # Log salmon version for reproducibility
    local salmon_version=$(salmon --version 2>&1 | head -1)
    log_message "  Salmon version: $salmon_version"
}

################################################################################
# Pre-flight Checks
################################################################################

check_disk_space() {
    local study_path="$1"
    local num_samples=$2

    # Salmon output: ~50-200MB per sample, with --dumpEq ~200-500MB
    local estimated_mb=$((num_samples * 500))
    local required_gb=$(( (estimated_mb / 1024) + 1 ))

    local available_gb=$(df -BG "$study_path" | tail -1 | awk '{print $4}' | sed 's/G//')

    log_message "Disk space check:"
    log_message "  Location: $study_path"
    log_message "  Available: ${available_gb}GB"
    log_message "  Estimated needed: ~${required_gb}GB (${num_samples} samples)"

    if [[ $available_gb -lt $required_gb ]]; then
        log_warning "Low disk space detected!"
        log_warning "  Available: ${available_gb}GB, Estimated need: ${required_gb}GB"
        log_warning "Proceeding anyway, but monitor disk usage..."
    else
        log_success "Sufficient disk space available"
    fi
}

################################################################################
# Sample Discovery (from Trimmed_data)
################################################################################

discover_samples() {
    local study_path="$1"
    local sample_filter="$2"
    local output_file="$3"

    local trimmed_data_dir="$study_path/Trimmed_data"

    if [[ ! -d "$trimmed_data_dir" ]]; then
        log_error "Trimmed_data directory not found: $trimmed_data_dir"
        log_error "Has fastp_trim_docker.sh been run for this study?"
        return 1
    fi

    > "$output_file"

    local all_samples=()
    local empty_dirs=()
    while IFS= read -r sample_dir; do
        local sample_name=$(basename "$sample_dir")

        # Skip FastQC and MultiQC directories
        [[ "$sample_name" == "FastQC" ]] && continue
        [[ "$sample_name" == "MultiQC" ]] && continue

        local fq_count=$(find "$sample_dir" -maxdepth 1 -name "fastp_*.fastq.gz" 2>/dev/null | wc -l)
        if [[ $fq_count -gt 0 ]]; then
            all_samples+=("$sample_name")
        else
            empty_dirs+=("$sample_name")
        fi
    done < <(find "$trimmed_data_dir" -mindepth 1 -maxdepth 1 -type d | sort)

    if [[ ${#empty_dirs[@]} -gt 0 ]]; then
        log_warning "${#empty_dirs[@]} sample director(ies) contain NO fastp_*.fastq.gz files:"
        for empty_dir in "${empty_dirs[@]}"; do
            log_warning "  - $empty_dir (empty or missing trimmed FASTQ files)"
        done
        log_warning "These may indicate failed fastp trimming - check fastp_trim logs"
    fi

    if [[ ${#all_samples[@]} -eq 0 ]]; then
        log_error "No samples with trimmed FASTQ files found in $trimmed_data_dir"
        return 1
    fi

    log_message "Found ${#all_samples[@]} sample(s) with trimmed FASTQ files in Trimmed_data/"

    # Apply filter if provided
    if [[ -n "$sample_filter" ]] && [[ "$sample_filter" != "all" ]]; then
        IFS=',' read -ra FILTER_ARRAY <<< "$sample_filter"

        local kept=0
        for filter_sample in "${FILTER_ARRAY[@]}"; do
            filter_sample=$(echo "$filter_sample" | tr -d '[:space:]')

            local found=false
            for sample in "${all_samples[@]}"; do
                if [[ "$sample" == "$filter_sample" ]]; then
                    echo "$sample" >> "$output_file"
                    kept=$((kept + 1))
                    found=true
                    break
                fi
            done

            if [[ "$found" == false ]]; then
                log_warning "Filtered sample not found in Trimmed_data: $filter_sample"
            fi
        done

        log_message "Selected $kept of ${#all_samples[@]} samples after filtering"

        if [[ $kept -eq 0 ]]; then
            log_error "No matching samples found after filtering"
            return 1
        fi
    else
        for sample in "${all_samples[@]}"; do
            echo "$sample" >> "$output_file"
        done
        log_message "Selected all ${#all_samples[@]} samples"
    fi

    return 0
}

################################################################################
# Skip Already-Quantified Samples
################################################################################

filter_already_quantified() {
    local study_path="$1"
    local samples_file="$2"

    if [[ "$FORCE_RERUN" == true ]]; then
        log_message "Force mode enabled - all samples will be (re)quantified"
        return 0
    fi

    local aligned_dir="$study_path/Aligned_data"
    local filtered_file="${samples_file}.filtered"
    > "$filtered_file"

    local total=0
    local skipped=0
    local to_process=0

    while IFS= read -r sample; do
        [[ -z "$sample" ]] && continue
        total=$((total + 1))

        local sample_out_dir="$aligned_dir/$sample"

        if [[ -f "$sample_out_dir/quant.sf" ]]; then
            skipped=$((skipped + 1))
            log_message "  - Skipping $sample (already quantified: quant.sf found)"
            mark_skipped "$sample"
        else
            echo "$sample" >> "$filtered_file"
            to_process=$((to_process + 1))
        fi
    done < "$samples_file"

    mv "$filtered_file" "$samples_file"

    if [[ $skipped -gt 0 ]]; then
        log_message ""
        log_message "Skip Summary:"
        log_message "  Total samples found: $total"
        log_message "  Already quantified (skipped): $skipped"
        log_message "  To process: $to_process"
        log_message ""
        log_message "  Use --force to re-quantify all samples"
    fi

    if [[ $to_process -eq 0 ]]; then
        log_success "All $total samples already quantified - nothing to do"
        log_message "Use --force to re-quantify"
    fi

    return 0
}

################################################################################
# Paired-End / Single-End Detection (from fastp output)
################################################################################

detect_read_layout() {
    local sample_dir="$1"

    # fastp produces:
    #   Paired:  fastp_*_1.fastq.gz  and  fastp_*_2.fastq.gz
    #   Single:  fastp_*.fastq.gz  (without _1 or _2 suffix)

    local r1_count=$(find "$sample_dir" -maxdepth 1 -name "fastp_*_1.fastq.gz" 2>/dev/null | wc -l)
    local r2_count=$(find "$sample_dir" -maxdepth 1 -name "fastp_*_2.fastq.gz" 2>/dev/null | wc -l)

    if [[ $r1_count -gt 0 ]] && [[ $r2_count -gt 0 ]]; then
        echo "paired"
    elif [[ $r1_count -gt 0 ]] && [[ $r2_count -eq 0 ]]; then
        echo "single"
    else
        local total_count=$(find "$sample_dir" -maxdepth 1 -name "fastp_*.fastq.gz" 2>/dev/null | wc -l)
        if [[ $total_count -gt 0 ]]; then
            echo "single"
        else
            echo "none"
        fi
    fi
}

record_peak_ram_usage() {
    local sample_name="$1"
    local phase_label="$2"
    local salmon_threads="$3"
    local run_status="$4"
    local time_log="$5"

    local peak_kb="NA"
    local peak_gb="NA"

    if [[ -f "$time_log" ]]; then
        peak_kb=$(grep -E 'Maximum resident set size \(kbytes\):' "$time_log" | awk '{print $6}' | tail -1)
        if [[ -n "$peak_kb" ]] && [[ "$peak_kb" =~ ^[0-9]+$ ]]; then
            peak_gb=$(awk -v kb="$peak_kb" 'BEGIN {printf "%.2f", kb/1024/1024}')
        else
            peak_kb="NA"
        fi
    fi

    if [[ -n "${RAM_USAGE_LOG:-}" ]]; then
        echo -e "${sample_name}\t${phase_label}\t${salmon_threads}\t${run_status}\t${peak_kb}\t${peak_gb}" >> "$RAM_USAGE_LOG"
    fi
}

run_salmon_with_monitor() {
    local sample_name="$1"
    local error_log="$2"
    local time_log="$3"
    shift 3

    local status=0
    local heartbeat_secs="${SALMON_HEARTBEAT_SEC:-60}"
    local show_salmon_stderr="${SALMON_SHOW_STDERR:-false}"
    [[ "$heartbeat_secs" =~ ^[0-9]+$ ]] || heartbeat_secs=60
    [[ "$heartbeat_secs" -lt 10 ]] && heartbeat_secs=10

    if command -v /usr/bin/time >/dev/null 2>&1; then
        if [[ "$show_salmon_stderr" == "true" ]]; then
            /usr/bin/time -v -o "$time_log" "$@" 2> >(tee -a "$error_log" >&2) &
        else
            /usr/bin/time -v -o "$time_log" "$@" 2>> "$error_log" &
        fi
    else
        if [[ "$show_salmon_stderr" == "true" ]]; then
            "$@" 2> >(tee -a "$error_log" >&2) &
        else
            "$@" 2>> "$error_log" &
        fi
        echo "WARNING: /usr/bin/time not available; peak RAM metrics unavailable" >> "$error_log"
    fi
    local salmon_pid=$!

    (
        while kill -0 "$salmon_pid" 2>/dev/null; do
            sleep "$heartbeat_secs"
            if kill -0 "$salmon_pid" 2>/dev/null; then
                echo "[$(date +%H:%M:%S)] salmon still running..."
            fi
        done
    ) &
    local heartbeat_pid=$!

    wait "$salmon_pid" || status=$?
    wait "$heartbeat_pid" 2>/dev/null || true

    return $status
}

################################################################################
# Salmon Quantification - Core Processing
################################################################################

quant_single_sample() {
    local sample_name="$1"
    local trimmed_dir="$2"
    local aligned_dir="$3"
    local salmon_index="$4"
    local salmon_threads="$5"
    local libtype_paired="$6"
    local libtype_single="$7"
    local extra_args="$8"
    local error_log="$9"
    local phase_label="${10:-phase1}"

    local sample_trim_dir="$trimmed_dir/$sample_name"
    local sample_out_dir="$aligned_dir/$sample_name"
    local sample_time_log="$aligned_dir/.salmon_time_${sample_name}_${phase_label}.log"
    local extra_args_arr=()
    if [[ -n "$extra_args" ]]; then
        # shellcheck disable=SC2206
        extra_args_arr=($extra_args)
    fi

    # If re-running, remove previous output so salmon doesn't complain
    if [[ -d "$sample_out_dir" ]]; then
        rm -rf "$sample_out_dir"
    fi
    mkdir -p "$aligned_dir"

    local layout=$(detect_read_layout "$sample_trim_dir")

    if [[ "$layout" == "paired" ]]; then
        # Collect ALL R1 and R2 files (handles technical replicates)
        # Sort to ensure matching order between R1 and R2
        local r1_files=()
        local r2_files=()

        while IFS= read -r r1_file; do
            r1_files+=("$r1_file")

            local r1_basename=$(basename "$r1_file")
            local r2_basename="${r1_basename/_1.fastq.gz/_2.fastq.gz}"
            local r2_file="$sample_trim_dir/$r2_basename"

            if [[ -f "$r2_file" ]]; then
                r2_files+=("$r2_file")
            else
                echo "WARNING: Missing R2 for $r1_basename in $sample_name" >> "$error_log"
                unset 'r1_files[-1]'
            fi
        done < <(find "$sample_trim_dir" -maxdepth 1 -name "fastp_*_1.fastq.gz" | sort)

        if [[ ${#r1_files[@]} -eq 0 ]]; then
            echo "ERROR: No valid paired-end files found for $sample_name" >> "$error_log"
            return 1
        fi

        echo "Paired-end quantification: ${#r1_files[@]} run(s)" >> "$error_log"
        echo "  R1: ${r1_files[*]}" >> "$error_log"
        echo "  R2: ${r2_files[*]}" >> "$error_log"

        if ! run_salmon_with_monitor \
            "$sample_name" \
            "$error_log" \
            "$sample_time_log" \
            salmon quant \
            -l "$libtype_paired" \
            -i "$salmon_index" \
            -1 "${r1_files[@]}" \
            -2 "${r2_files[@]}" \
            "${extra_args_arr[@]}" \
            --threads "$salmon_threads" \
            -o "$sample_out_dir"; then
            record_peak_ram_usage "$sample_name" "$phase_label" "$salmon_threads" "failed" "$sample_time_log"
            return 1
        fi

    elif [[ "$layout" == "single" ]]; then
        local se_files=()

        while IFS= read -r fq_file; do
            se_files+=("$fq_file")
        done < <(find "$sample_trim_dir" -maxdepth 1 -name "fastp_*.fastq.gz" \
                     ! -name "fastp_*_2.fastq.gz" | sort)

        if [[ ${#se_files[@]} -eq 0 ]]; then
            echo "ERROR: No single-end files found for $sample_name" >> "$error_log"
            return 1
        fi

        echo "Single-end quantification: ${#se_files[@]} run(s)" >> "$error_log"
        echo "  Files: ${se_files[*]}" >> "$error_log"

        if ! run_salmon_with_monitor \
            "$sample_name" \
            "$error_log" \
            "$sample_time_log" \
            salmon quant \
            -l "$libtype_single" \
            -i "$salmon_index" \
            -r "${se_files[@]}" \
            "${extra_args_arr[@]}" \
            --threads "$salmon_threads" \
            -o "$sample_out_dir"; then
            record_peak_ram_usage "$sample_name" "$phase_label" "$salmon_threads" "failed" "$sample_time_log"
            return 1
        fi

    else
        echo "No trimmed FASTQ files found in $sample_trim_dir" >> "$error_log"
        record_peak_ram_usage "$sample_name" "$phase_label" "$salmon_threads" "failed" "$sample_time_log"
        return 1
    fi

    # Verify output
    if [[ ! -f "$sample_out_dir/quant.sf" ]]; then
        echo "ERROR: quant.sf not generated for $sample_name" >> "$error_log"
        record_peak_ram_usage "$sample_name" "$phase_label" "$salmon_threads" "failed" "$sample_time_log"
        return 1
    fi

    record_peak_ram_usage "$sample_name" "$phase_label" "$salmon_threads" "success" "$sample_time_log"

    return 0
}

################################################################################
# Parallel Quantification Orchestrator (with 3-Phase Retry)
################################################################################

run_quantification() {
    local study_path="$1"
    local samples_file="$2"

    log_message "========================================="
    log_message "STEP 1: Salmon Quantification"
    log_message "========================================="

    local trimmed_dir="$study_path/Trimmed_data"
    local aligned_dir="$study_path/Aligned_data"
    mkdir -p "$aligned_dir"

    local total_samples=$(wc -l < "$samples_file")

    # Calculate optimal parallelism (RAM-aware)
    SALMON_THREADS=$(calculate_salmon_threads)
    local parallel_jobs_raw
    parallel_jobs_raw="$(calculate_optimal_parallelism "$total_samples")"
    local parallel_jobs

    # Defensive parsing: keep only numeric value in case command substitution
    # accidentally captures log lines from helper functions.
    parallel_jobs="$(echo "$parallel_jobs_raw" | grep -Eo '[0-9]+' | tail -1)"
    if [[ -z "$parallel_jobs" ]]; then
        log_warning "Could not parse parallel_jobs from calculation output, defaulting to 1"
        parallel_jobs=1
    fi

    # Safety check: ensure we don't oversubscribe CPUs
    local total_cpu_usage=$((parallel_jobs * SALMON_THREADS))
    if [[ $total_cpu_usage -gt $CPUS ]]; then
        log_warning "CPU oversubscription detected! Reducing parallel jobs..."
        log_warning "  Would use: $total_cpu_usage cores, Available: $CPUS cores"
        parallel_jobs=$((CPUS / SALMON_THREADS))
        [[ $parallel_jobs -lt 1 ]] && parallel_jobs=1
        total_cpu_usage=$((parallel_jobs * SALMON_THREADS))
        log_message "  Adjusted to: $parallel_jobs parallel jobs (${total_cpu_usage} cores)"
    fi

    # Safety check for parallel_jobs
    if [[ -z "$parallel_jobs" ]] || ! [[ "$parallel_jobs" =~ ^[0-9]+$ ]]; then
        log_warning "parallel_jobs not set properly, defaulting to 1"
        parallel_jobs=1
        total_cpu_usage=$((parallel_jobs * SALMON_THREADS))
    fi

    # If RAM limits us to one job, prioritize speed by using more threads in that one job.
    if [[ $parallel_jobs -eq 1 ]]; then
        local boosted_threads=$CPUS
        [[ $boosted_threads -gt 12 ]] && boosted_threads=12
        [[ $boosted_threads -lt 2 ]] && boosted_threads=2

        if [[ $boosted_threads -gt $SALMON_THREADS ]]; then
            log_message "  Single-job mode detected; increasing threads per job: $SALMON_THREADS -> $boosted_threads"
            SALMON_THREADS=$boosted_threads
            total_cpu_usage=$((parallel_jobs * SALMON_THREADS))
        fi
    fi

    # DYNAMIC: distinguish estimated required RAM from assigned per-job RAM budget
    local estimated_required_mem_per_job=$(calculate_required_mem_per_job "$SALMON_THREADS")
    local assigned_mem_per_job=$(calculate_mem_per_job "$parallel_jobs" "$AVAILABLE_MEM_GB")

    log_message "Adaptive quantification configuration:"
    log_message "  Salmon index: $SALMON_INDEX"
    log_message "  Samples to process: $total_samples"
    log_message "  Parallel jobs: $parallel_jobs (based on $CPUS CPUs, ${AVAILABLE_MEM_GB}GB RAM)"
    log_message "  Salmon threads per job: $SALMON_THREADS"
    log_message "  Estimated required RAM/job: ~${estimated_required_mem_per_job}GB (index + working memory)"
    log_message "  Assigned RAM budget/job: ~${assigned_mem_per_job}GB (dynamic share for $parallel_jobs job(s))"
    log_message "  Total CPU usage: ~${total_cpu_usage}/${CPUS} cores"
    log_message "  Library type (PE): $SALMON_LIBTYPE_PAIRED"
    log_message "  Library type (SE): $SALMON_LIBTYPE_SINGLE"
    log_message "  Extra args: $SALMON_EXTRA_ARGS"
    log_message ""
    log_message "Errors will be logged, processing continues for other samples"
    log_message ""

    ############################################################################
    # PHASE 1: Parallel quantification with adaptive resources
    ############################################################################

    log_message "--- Phase 1: Parallel quantification ($parallel_jobs jobs) ---"
    log_message ""

    # Export functions and variables for GNU parallel
    export -f detect_read_layout record_peak_ram_usage run_salmon_with_monitor quant_single_sample
    export -f log_message log_error log_step_error mark_quant_failed mark_successful
    export ERROR_LOG QUANT_FAILED SUCCESSFUL_SAMPLES RAM_USAGE_LOG
    export SALMON_THREADS SALMON_LIBTYPE_PAIRED SALMON_LIBTYPE_SINGLE SALMON_EXTRA_ARGS SALMON_INDEX
    export SALMON_SHOW_STDERR SALMON_HEARTBEAT_SEC

    cat "$samples_file" | \
        parallel -j $parallel_jobs --line-buffer --tagstring '[{1}]' \
        '
        sample={1}
        trimmed_dir="'"$trimmed_dir"'"
        aligned_dir="'"$aligned_dir"'"

        echo "[$(date +%H:%M:%S)] Starting salmon quantification..."

        # Detect layout for logging
        layout=$(detect_read_layout "$trimmed_dir/$sample")
        echo "[$(date +%H:%M:%S)] Layout: $layout"

        # Count input files for logging
        if [[ "$layout" == "paired" ]]; then
            r1_count=$(find "$trimmed_dir/$sample" -maxdepth 1 -name "fastp_*_1.fastq.gz" 2>/dev/null | wc -l)
            echo "[$(date +%H:%M:%S)] Technical replicates: $r1_count run(s)"
        elif [[ "$layout" == "single" ]]; then
            se_count=$(find "$trimmed_dir/$sample" -maxdepth 1 -name "fastp_*.fastq.gz" ! -name "fastp_*_2.fastq.gz" 2>/dev/null | wc -l)
            echo "[$(date +%H:%M:%S)] Technical replicates: $se_count run(s)"
        fi

        if quant_single_sample \
            "$sample" \
            "$trimmed_dir" \
            "$aligned_dir" \
            "'"$SALMON_INDEX"'" \
            "'"$SALMON_THREADS"'" \
            "'"$SALMON_LIBTYPE_PAIRED"'" \
            "'"$SALMON_LIBTYPE_SINGLE"'" \
            "'"$SALMON_EXTRA_ARGS"'" \
            "'"$ERROR_LOG"'" \
            "phase1"; then

            echo "[$(date +%H:%M:%S)] ✓ Complete"
            mark_successful "$sample"
        else
            echo "[$(date +%H:%M:%S)] ✗ Failed"
            log_step_error "QUANT" "$sample" "salmon quantification failed (Phase 1: parallel)"
            mark_quant_failed "$sample"
        fi
        '

    # Count Phase 1 results
    local successful=0
    local failed_count=0

    [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
    [[ -f "$QUANT_FAILED" ]] && failed_count=$(wc -l < "$QUANT_FAILED")

    log_message ""
    log_message "Phase 1 Summary:"
    log_message "  Successful: $successful"
    log_message "  Failed: $failed_count"

    ############################################################################
    # PHASE 2: Sequential retry with maximum resources
    # Gives each failed sample ALL CPUs + ALL RAM (like geo_download_docker.sh
    # sequential retry that gives maximum memory to a single job)
    ############################################################################

    if [[ $failed_count -gt 0 ]]; then
        log_message ""
        log_warning "Initial parallel quantification had $failed_count failures"
        log_message "Retrying failed samples sequentially (one at a time) for reliability..."
        log_message ""

        # Create retry list from failed samples
        local retry_list=$(mktemp)
        cp "$QUANT_FAILED" "$retry_list"
        > "$QUANT_FAILED"  # Clear for retry tracking

        local retry_count=$(wc -l < "$retry_list")
        local retry_success=0
        local retry_fail=0

        # Sequential retry: give ALL CPUs to one sample (like geo_download Phase 2)
        local retry_threads=$CPUS
        [[ $retry_threads -gt 16 ]] && retry_threads=16  # Salmon diminishing returns past ~16

        # DYNAMIC: Recalculate memory for sequential mode (1 job = maximum memory!)
        local retry_mem=$(calculate_mem_per_job 1 $AVAILABLE_MEM_GB)

        log_message "Phase 2 configuration:"
        log_message "  Jobs: 1 (sequential)"
        log_message "  Threads per job: $retry_threads (maximum available)"
        log_message "  Memory per job: ~${retry_mem}GB (maximum available for single job)"
        log_message "  Samples to retry: $retry_count"
        log_message ""

        while IFS= read -r sample; do
            [[ -z "$sample" ]] && continue

            log_message "[RETRY] Processing $sample (sequential, threads: $retry_threads, mem: ~${retry_mem}GB)..."

            if quant_single_sample \
                "$sample" \
                "$trimmed_dir" \
                "$aligned_dir" \
                "$SALMON_INDEX" \
                "$retry_threads" \
                "$SALMON_LIBTYPE_PAIRED" \
                "$SALMON_LIBTYPE_SINGLE" \
                "$SALMON_EXTRA_ARGS" \
                "$ERROR_LOG" \
                "phase2"; then

                log_message "[RETRY] ✓ $sample succeeded"
                mark_successful "$sample"
                retry_success=$((retry_success + 1))
            else
                log_error "[RETRY] ✗ $sample failed again"
                mark_quant_failed "$sample"
                retry_fail=$((retry_fail + 1))
            fi
        done < "$retry_list"

        rm -f "$retry_list"

        log_message ""
        log_message "Phase 2 (Sequential Retry) Results:"
        log_message "  Succeeded on retry: $retry_success"
        log_message "  Still failed: $retry_fail"

        ########################################################################
        # PHASE 3: Minimal-memory fallback
        # For samples that fail even with full resources, the issue is likely
        # memory pressure from bias models. Salmon's --seqBias, --gcBias, and
        # --posBias each add significant memory overhead. Disabling them and
        # using minimal threads reduces memory footprint drastically.
        #
        # This is analogous to geo_download_docker.sh's Phase 3 disk-temp
        # fallback, which switches from RAM disk to disk to handle samples
        # that need more space than RAM provides.
        ########################################################################

        if [[ $retry_fail -gt 0 ]]; then
            log_message ""
            log_warning "Phase 3: ${retry_fail} samples still failing - switching to minimal-memory mode"
            log_message "This disables bias models (--seqBias, --gcBias, --posBias) to reduce RAM usage."
            log_message "Results will be usable but less accurate than full bias-corrected quantification."
            log_message ""

            # Create list of samples that failed Phase 2
            local phase3_list=$(mktemp)
            cp "$QUANT_FAILED" "$phase3_list"
            > "$QUANT_FAILED"  # Clear for final tracking

            local phase3_count=$(wc -l < "$phase3_list")
            local phase3_success=0
            local phase3_fail=0

            # Minimal configuration: 2 threads, no bias models
            local phase3_threads=2

            # Strip bias-related flags from extra args to reduce memory
            local phase3_extra_args=$(echo "$SALMON_EXTRA_ARGS" | \
                sed 's/--seqBias//g; s/--gcBias//g; s/--posBias//g' | \
                tr -s ' ')

            # DYNAMIC: Recalculate memory for minimal mode
            local phase3_mem=$(calculate_mem_per_job 1 $AVAILABLE_MEM_GB)

            log_message "Phase 3 configuration:"
            log_message "  Jobs: 1 (sequential)"
            log_message "  Threads: $phase3_threads (minimal)"
            log_message "  Memory per job: ~${phase3_mem}GB (maximum available)"
            log_message "  Bias models: DISABLED (reduced memory footprint)"
            log_message "  Extra args: $phase3_extra_args"
            log_message "  Samples to retry: $phase3_count"
            log_message ""

            while IFS= read -r sample; do
                [[ -z "$sample" ]] && continue

                log_message "[MINIMAL] Processing $sample (threads: $phase3_threads, no bias models)..."

                if quant_single_sample \
                    "$sample" \
                    "$trimmed_dir" \
                    "$aligned_dir" \
                    "$SALMON_INDEX" \
                    "$phase3_threads" \
                    "$SALMON_LIBTYPE_PAIRED" \
                    "$SALMON_LIBTYPE_SINGLE" \
                    "$phase3_extra_args" \
                    "$ERROR_LOG" \
                    "phase3"; then

                    log_message "[MINIMAL] ✓ $sample succeeded (note: quantified WITHOUT bias correction)"
                    echo "${sample}" >> "$SUCCESSFUL_SAMPLES"

                    # Write a notice file so downstream analysis knows this sample lacks bias correction
                    local notice_file="$aligned_dir/$sample/NO_BIAS_CORRECTION.txt"
                    cat > "$notice_file" << NOTICE_EOF
WARNING: This sample was quantified WITHOUT bias correction models.

The standard quantification (with --seqBias --gcBias --posBias) failed,
likely due to memory constraints. This sample was re-run in minimal-memory
mode (Phase 3) with bias models disabled.

The quantification is usable but may be less accurate than samples processed
with full bias correction. Consider re-running with more available RAM:
  salmon quant -l ${SALMON_LIBTYPE_PAIRED} -i ${SALMON_INDEX} ... $SALMON_EXTRA_ARGS

Phase 1 args: $SALMON_EXTRA_ARGS
Phase 3 args: $phase3_extra_args
Date: $(date)
NOTICE_EOF
                    phase3_success=$((phase3_success + 1))
                else
                    log_error "[MINIMAL] ✗ $sample failed even in minimal mode"
                    log_error "This sample may have corrupted input files or other issues"
                    mark_quant_failed "$sample"
                    phase3_fail=$((phase3_fail + 1))
                fi
            done < "$phase3_list"

            rm -f "$phase3_list"

            log_message ""
            log_message "Phase 3 Results (minimal-memory fallback):"
            log_message "  Succeeded (without bias correction): $phase3_success"
            log_message "  Still failed (likely corrupted): $phase3_fail"

            if [[ $phase3_success -gt 0 ]]; then
                log_warning "$phase3_success sample(s) quantified WITHOUT bias correction"
                log_warning "Check Aligned_data/*/NO_BIAS_CORRECTION.txt for affected samples"
            fi
        fi

        # Update final counts
        successful=0
        failed_count=0
        [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
        [[ -f "$QUANT_FAILED" ]] && failed_count=$(wc -l < "$QUANT_FAILED")
    fi

    log_message ""
    log_message "========================================="
    log_message "Quantification Summary (all phases):"
    log_message "  Quantified successfully: $successful"
    log_message "  Failed quantification: $failed_count"
    log_message "========================================="

    if [[ $failed_count -gt 0 ]]; then
        log_warning "Some samples failed quantification - see $QUANT_FAILED"
    else
        log_success "All samples quantified successfully"
    fi
}

################################################################################
# Summary Generation
################################################################################

generate_summary() {
    local study_path="$1"

    log_message "Generating quantification summary..."

    local aligned_dir="$study_path/Aligned_data"

    local sample_count=$(find "$aligned_dir" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
    local quant_count=$(find "$aligned_dir" -name "quant.sf" 2>/dev/null | wc -l)
    local total_size=$(du -sh "$aligned_dir" 2>/dev/null | cut -f1)
    local no_bias_count=$(find "$aligned_dir" -name "NO_BIAS_CORRECTION.txt" 2>/dev/null | wc -l)

    local quant_errors=0
    local successful=0
    local skipped=0
    local peak_ram_kb=""
    local peak_ram_gb="N/A"
    local peak_ram_sample="N/A"
    local peak_ram_phase="N/A"

    [[ -f "$QUANT_FAILED" ]] && quant_errors=$(wc -l < "$QUANT_FAILED")
    [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
    [[ -f "$SKIPPED_SAMPLES" ]] && skipped=$(wc -l < "$SKIPPED_SAMPLES")
    if [[ -s "$RAM_USAGE_LOG" ]]; then
        read -r peak_ram_kb peak_ram_sample peak_ram_phase < <(
            awk -F'\t' '
                $5 ~ /^[0-9]+$/ {
                    if ($5 > max) {
                        max = $5
                        sample = $1
                        phase = $2
                    }
                }
                END {
                    if (max > 0) print max, sample, phase
                }
            ' "$RAM_USAGE_LOG"
        )
        if [[ -n "$peak_ram_kb" ]] && [[ "$peak_ram_kb" =~ ^[0-9]+$ ]]; then
            peak_ram_gb=$(awk -v kb="$peak_ram_kb" 'BEGIN {printf "%.2f", kb/1024/1024}')
        fi
    fi

    # Collect library type info from salmon outputs
    local lib_summary=""
    if [[ -d "$aligned_dir" ]]; then
        lib_summary=$(
            for lf in "$aligned_dir"/*/lib_format_counts.json; do
                [[ -f "$lf" ]] || continue
                local sname=$(basename "$(dirname "$lf")")
                local expected=$(grep -o '"expected_format"[[:space:]]*:[[:space:]]*"[^"]*"' "$lf" 2>/dev/null | head -1 | cut -d'"' -f4)
                local compatible=$(grep -o '"compatible_fragment_ratio"[[:space:]]*:[[:space:]]*[0-9.]*' "$lf" 2>/dev/null | head -1 | awk -F: '{print $2}' | tr -d ' ')
                if [[ -n "$expected" ]]; then
                    echo "  $sname: format=$expected compatible_ratio=${compatible:-N/A}"
                fi
            done
        )
    fi

    cat > "$study_path/salmon_quant_summary.txt" << EOF
Salmon Quantification Summary
================================================================================
Study: $(basename "$study_path")
Date: $(date)
System CPUs: $TOTAL_CPUS
CPUs Used: $CPUS
Salmon Threads per Job: $SALMON_THREADS
Available Memory: ${AVAILABLE_MEM_GB}GB
Estimated Index Memory: ~${INDEX_MEM_GB}GB
Observed Peak RAM (RSS): ${peak_ram_gb}GB (sample: ${peak_ram_sample}, phase: ${peak_ram_phase})

Parameters:
-----------
Salmon Index: $SALMON_INDEX
Library Type (PE): $SALMON_LIBTYPE_PAIRED
Library Type (SE): $SALMON_LIBTYPE_SINGLE
Extra Args: $SALMON_EXTRA_ARGS
Force Re-run: $(if [[ "$FORCE_RERUN" == true ]]; then echo "Yes"; else echo "No"; fi)

Results:
--------
Sample Directories: $sample_count
quant.sf Files: $quant_count
Total Size: ${total_size:-Unknown}

Processing Status:
------------------
Successful (quantified): $successful
Skipped (already quantified): $skipped
Failed at Quantification: $quant_errors
Quantified without bias correction (Phase 3): $no_bias_count

Retry Strategy:
---------------
Phase 1: Parallel with adaptive resources (CPU+RAM aware)
Phase 2: Sequential retry with maximum resources (all CPUs + all RAM for 1 job)
Phase 3: Minimal-memory fallback (2 threads, bias models disabled)

Library Format Detection:
-------------------------
${lib_summary:-No library format data available}

Error Logs:
-----------
Main error log: $ERROR_LOG
Quant failures: $QUANT_FAILED
RAM usage log: $RAM_USAGE_LOG

Quantification Output Files:
-----------------------------
EOF

    find "$aligned_dir" -name "quant.sf" | sort >> "$study_path/salmon_quant_summary.txt"

    if [[ $no_bias_count -gt 0 ]]; then
        echo "" >> "$study_path/salmon_quant_summary.txt"
        echo "Samples Quantified WITHOUT Bias Correction:" >> "$study_path/salmon_quant_summary.txt"
        echo "---------------------------------------------" >> "$study_path/salmon_quant_summary.txt"
        find "$aligned_dir" -name "NO_BIAS_CORRECTION.txt" -exec dirname {} \; | \
            xargs -I{} basename {} >> "$study_path/salmon_quant_summary.txt"
    fi

    log_success "Summary saved to $study_path/salmon_quant_summary.txt"

    echo ""
    log_message "========================================="
    log_message "QUANTIFICATION COMPLETED"
    log_message "========================================="
    log_message "Study: $(basename "$study_path")"
    log_message "Samples: $sample_count"
    log_message "quant.sf files: $quant_count"
    log_message "Total size: ${total_size:-Unknown}"
    log_message ""
    log_message "Processing Summary:"
    log_message "  Successful: $successful"
    log_message "  Skipped (already quantified): $skipped"
    log_message "  Failed at quantification: $quant_errors"
    log_message "  Peak observed RAM (RSS): ${peak_ram_gb}GB (sample: ${peak_ram_sample}, phase: ${peak_ram_phase})"
    if [[ $no_bias_count -gt 0 ]]; then
        log_warning "  Without bias correction (Phase 3): $no_bias_count"
    fi
    log_message ""
    log_message "Location: $study_path/Aligned_data/"
    log_message "Error log: $ERROR_LOG"
    log_message "========================================="

    if [[ $quant_errors -gt 0 ]]; then
        log_warning "Some samples failed processing - check error logs for details"
    fi
}

################################################################################
# Cleanup Handler
################################################################################

cleanup_on_exit() {
    local exit_code=$?

    rm -f /tmp/salmon_samples_*.txt 2>/dev/null

    if [[ -n "${STUDY_DIR:-}" ]] && [[ $exit_code -ne 0 ]]; then
        if [[ -n "${ERROR_LOG:-}" ]]; then
            log_error "Script exited with errors (code: $exit_code)"
            log_error "Check $ERROR_LOG for details"
        fi
    fi
}

trap cleanup_on_exit EXIT INT TERM

################################################################################
# Single Study Processing
################################################################################

process_single_study() {
    local study_dir="$1"
    local sample_filter="${2:-}"

    local study_path="$BASE_PATH/$study_dir"

    log_message "========================================="
    log_message "Processing Study: $study_dir"
    log_message "========================================="
    log_message "Study Path: $study_path"
    log_message "Salmon Index: $SALMON_INDEX"
    log_message "Sample Filter: ${sample_filter:-all}"
    log_message "========================================="

    if [[ ! -d "$study_path" ]]; then
        log_error "Study directory not found: $study_path"
        return 1
    fi

    if [[ ! -d "$study_path/Trimmed_data" ]]; then
        log_error "Trimmed_data directory not found in $study_path"
        log_error "Has fastp_trim_docker.sh been run for this study?"
        return 1
    fi

    detect_system_resources
    check_dependencies

    initialize_error_tracking "$study_path"

    local samples_file="$study_path/salmon_samples_to_process.txt"
    if ! discover_samples "$study_path" "$sample_filter" "$samples_file"; then
        log_error "Failed to discover samples"
        return 1
    fi

    filter_already_quantified "$study_path" "$samples_file"

    local num_samples=$(wc -l < "$samples_file")

    if [[ $num_samples -eq 0 ]]; then
        rm -f "$samples_file"
        generate_summary "$study_path"
        return 0
    fi

    check_disk_space "$study_path" "$num_samples"

    # Run salmon quantification (includes 3-phase retry)
    run_quantification "$study_path" "$samples_file"

    rm -f "$samples_file"
    generate_summary "$study_path"

    log_success "Study processing completed: $study_dir"
    return 0
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

    local aggregated=$(mktemp)

    awk -F'\t' '
        /^#/ { next }
        /^[[:space:]]*$/ { next }
        {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $3)

            if ($1 == "" || $3 == "") next

            study = $1

            if ($3 == "all") {
                studies[study] = "all"
            } else if (study in studies && studies[study] == "all") {
                # Already set to all, skip
            } else if (study in studies) {
                studies[study] = studies[study] "," $3
            } else {
                studies[study] = $3
                order[++n] = study
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                printf "%s\t%s\n", order[i], studies[order[i]]
            }
        }
    ' "$INPUT_FILE" > "$aggregated"

    local total_studies=$(wc -l < "$aggregated")
    log_message "Found $total_studies unique study/studies to process"
    log_message "Processing $PARALLEL_DATASETS study/studies at a time"
    log_message "========================================="

    if [[ $PARALLEL_DATASETS -eq 1 ]]; then
        local study_num=0
        local success_count=0
        local fail_count=0

        set +e

        while IFS=$'\t' read -r study_dir gsm_list || [[ -n "$study_dir" ]]; do
            [[ -z "$study_dir" ]] && continue
            ((study_num++))

            log_message ""
            log_message "========================================="
            log_message "Study $study_num of $total_studies: $study_dir"
            log_message "  Samples: ${gsm_list:0:100}..."
            log_message "========================================="

            local filter=""
            [[ "$gsm_list" != "all" ]] && filter="$gsm_list"

            if process_single_study "$study_dir" "$filter"; then
                ((success_count++))
                log_success "Study $study_num completed: $study_dir"
            else
                ((fail_count++))
                log_error "Study $study_num failed: $study_dir"
            fi

        done < "$aggregated"

        set -e

        rm -f "$aggregated"

        log_message ""
        log_message "========================================="
        log_message "BATCH PROCESSING COMPLETE"
        log_message "========================================="
        log_message "Total studies: $study_num"
        log_message "Successful: $success_count"
        log_message "Failed: $fail_count"
        log_message "========================================="
    else
        log_message "Starting parallel processing with $PARALLEL_DATASETS jobs..."

        export -f log_message log_error log_warning log_success log_step_error
        export -f detect_system_resources check_dependencies estimate_index_memory
        export -f initialize_error_tracking mark_quant_failed mark_successful mark_skipped is_sample_failed
        export -f discover_samples filter_already_quantified detect_read_layout
        export -f record_peak_ram_usage run_salmon_with_monitor quant_single_sample
        export -f calculate_optimal_parallelism calculate_salmon_threads calculate_required_mem_per_job calculate_mem_per_job
        export -f check_disk_space run_quantification generate_summary
        export -f process_single_study
        export BASE_PATH PARALLEL_DATASETS FORCE_RERUN
        export SALMON_INDEX SALMON_THREADS SALMON_LIBTYPE_PAIRED SALMON_LIBTYPE_SINGLE SALMON_EXTRA_ARGS

        cat "$aggregated" | \
            parallel --colsep '\t' -j $PARALLEL_DATASETS --line-buffer \
            '
            study_dir={1}
            gsm_list={2}

            filter=""
            [[ "$gsm_list" != "all" ]] && filter="$gsm_list"

            bash "'"$0"'" \
                --index "'"$SALMON_INDEX"'" \
                "$study_dir" "$filter" \
                $(if [[ "'"$FORCE_RERUN"'" == true ]]; then echo "--force"; fi) \
                --threads "'"$SALMON_THREADS"'" \
                --libtype-pe "'"$SALMON_LIBTYPE_PAIRED"'" \
                --libtype-se "'"$SALMON_LIBTYPE_SINGLE"'" \
                --extra-args "'"$SALMON_EXTRA_ARGS"'" \
                --base-path "'"$BASE_PATH"'"
            '

        rm -f "$aggregated"

        log_message ""
        log_message "========================================="
        log_success "PARALLEL PROCESSING COMPLETE"
        log_message "========================================="
    fi
}

################################################################################
# Main Entry Point
################################################################################

main() {
    parse_arguments "$@"

    if [[ -n "$INPUT_FILE" ]]; then
        process_input_file
    else
        process_single_study "$STUDY_DIR" "$SAMPLE_FILTER"
    fi
}

main "$@"
