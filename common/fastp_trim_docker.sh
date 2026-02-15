#!/bin/bash

################################################################################
# FASTQ Quality Control & Trimming Pipeline
################################################################################
#
# Purpose:
#   Performs adapter trimming, quality filtering, and optional QC reporting
#   on FASTQ files produced by the GEO/SRA download pipeline
#   (geo_download_docker.sh). Automatically detects paired-end vs single-end
#   layout, runs fastp with sensible defaults, and optionally generates
#   per-file FastQC reports aggregated into a single MultiQC summary.
#
# Docker image: stodo3569/fastp_tools:0.0
#
# Usage:
#   docker run --rm \
#     -u "$(id -u):$(id -g)" \
#     -v /data:/data \
#     stodo3569/fastp_tools:0.0 \
#     bash /data/<scripts_dir>/fastp_trim_docker.sh [OPTIONS]
#
# Modes:
#   # Single study — trim all samples
#   bash fastp_trim_docker.sh "Study_Name"
#
#   # Single study — specific GSM samples
#   bash fastp_trim_docker.sh "Study_Name" "GSM111,GSM222"
#
#   # Batch mode — multiple studies from input file
#   bash fastp_trim_docker.sh --input-file datasets.txt
#
#   # Batch mode — parallel processing (2 studies simultaneously)
#   bash fastp_trim_docker.sh --input-file datasets.txt --parallel 2
#
#   # With FastQC + MultiQC reporting
#   bash fastp_trim_docker.sh --input-file datasets.txt --fastqc
#
# Input file format (tab-separated, same as geo_download_docker.sh):
#   Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC_file]
#   Lines with the same Output_Directory are automatically aggregated.
#
# Pipeline stages:
#   1. Trim      — adapter removal, quality filtering, polyX trimming (fastp)
#   2. FastQC    — per-file quality reports (optional, --fastqc)
#   3. MultiQC   — aggregate HTML report across all samples (optional, --fastqc)
#   Failed samples are automatically retried sequentially with full resources.
#
# Default fastp parameters:
#   --trim_poly_x              Trim polyX tails
#   --correction               Base correction for overlapping paired reads
#   --detect_adapter_for_pe    Auto-detect adapters for paired-end data
#   --length_required 36       Discard reads shorter than 36 bp
#
# Outputs (written to /data/<Output_Directory>/):
#   Trimmed_data/<GSM>/fastp_<SRR>_1.fastq.gz  — trimmed R1 (paired-end)
#   Trimmed_data/<GSM>/fastp_<SRR>_2.fastq.gz  — trimmed R2 (paired-end)
#   Trimmed_data/<GSM>/<SRR>.html               — fastp HTML report per run
#   Trimmed_data/<GSM>/<SRR>.json               — fastp JSON report per run
#   Trimmed_data/FastQC/                         — FastQC reports (if --fastqc)
#   Trimmed_data/MultiQC/multiqc_report.html     — MultiQC summary (if --fastqc)
#   fastp_trim_summary.txt                       — summary of results
#   fastp_error_log.txt                          — detailed error log
#   fastp_trim_failed.txt                        — samples that failed trimming
#   fastp_fastqc_failed.txt                      — samples that failed FastQC
#   fastp_successful_samples.txt                 — samples that completed all stages
#   fastp_skipped_samples.txt                    — samples skipped (already trimmed)
#
# Dependencies (provided by stodo3569/fastp_tools Docker image):
#   fastp (>=0.23.0):   adapter trimming and quality filtering
#   FastQC (>=0.12.0):  per-file quality reports (optional)
#   MultiQC (>=1.14):   aggregate QC summary (optional)
#   GNU parallel:       parallel sample processing
#   pigz:               parallel gzip compression
#
# Resource handling:
#   CPU and memory are detected automatically and allocated adaptively.
#   fastp threads per job scale with available CPUs (2–6 threads),
#   and the number of parallel jobs is calculated to avoid oversubscription.
#
# Expected input structure (created by geo_download_docker.sh):
#
#   /data/
#   └── <Output_Directory>/              e.g. Friedman_2019_GSE125050
#       └── Raw_data/
#           ├── <GSM>/                   e.g. GSM3559136
#           │   ├── <SRR>_1.fastq.gz    e.g. SRR8571953_1.fastq.gz
#           │   └── <SRR>_2.fastq.gz    e.g. SRR8571953_2.fastq.gz
#           └── ...
#
# Note on /data:
#   Same mount point as the download script. All study directories,
#   raw FASTQ inputs, and trimmed outputs live under /data/<Output_Directory>/.
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

    log_message "System Resources Detected:"
    log_message "  Total CPUs: $TOTAL_CPUS"
    log_message "  CPUs for this dataset: $CPUS"
    log_message "  Total Memory: ${TOTAL_MEM_GB}GB"
    log_message "  Available Memory: ${AVAILABLE_MEM_GB}GB"
}

################################################################################
# Adaptive Resource Calculation
################################################################################

calculate_optimal_parallelism() {
    local num_samples=${1:-10}

    # fastp is moderately CPU-bound. Each job gets FASTP_THREADS CPUs.
    # Strategy: total CPUs / threads_per_job = parallel jobs
    local optimal_jobs=$((CPUS / FASTP_THREADS))

    # Ensure reasonable bounds
    [[ $optimal_jobs -lt 1 ]] && optimal_jobs=1

    # Don't spawn more jobs than samples
    if [[ $num_samples -lt $optimal_jobs ]]; then
        echo $num_samples
    else
        echo $optimal_jobs
    fi
}

calculate_fastp_threads() {
    # fastp thread allocation per job
    # fastp is efficient with 2-8 threads; diminishing returns beyond that
    if [[ $CPUS -ge 32 ]]; then
        echo 6
    elif [[ $CPUS -ge 16 ]]; then
        echo 4
    elif [[ $CPUS -ge 8 ]]; then
        echo 4
    elif [[ $CPUS -ge 4 ]]; then
        echo 2
    else
        echo 2
    fi
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
RUN_FASTQC=false
RUN_MULTIQC=false

CPUS=4
TOTAL_CPUS=4
TOTAL_MEM_GB=8
AVAILABLE_MEM_GB=6

# fastp parameters (from fastp_trim_1.sh)
FASTP_THREADS=4
FASTP_LENGTH_REQUIRED=36
FASTP_EXTRA_ARGS="--trim_poly_x --correction --detect_adapter_for_pe"

# Error tracking files (set per-study in initialize_error_tracking)
ERROR_LOG=""
TRIM_FAILED=""
FASTQC_FAILED=""
SUCCESSFUL_SAMPLES=""
SKIPPED_SAMPLES=""

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

    ERROR_LOG="${study_path}/fastp_error_log.txt"
    TRIM_FAILED="${study_path}/fastp_trim_failed.txt"
    FASTQC_FAILED="${study_path}/fastp_fastqc_failed.txt"
    SUCCESSFUL_SAMPLES="${study_path}/fastp_successful_samples.txt"
    SKIPPED_SAMPLES="${study_path}/fastp_skipped_samples.txt"

    > "$ERROR_LOG"
    > "$TRIM_FAILED"
    > "$FASTQC_FAILED"
    > "$SUCCESSFUL_SAMPLES"
    > "$SKIPPED_SAMPLES"

    log_message "Error tracking initialized"
    log_message "  Error log: $ERROR_LOG"
    log_message "  Failed files will be tracked per stage"
}

mark_trim_failed() {
    local sample="$1"
    echo "${sample}" >> "$TRIM_FAILED"
}

mark_fastqc_failed() {
    local sample="$1"
    echo "${sample}" >> "$FASTQC_FAILED"
}

mark_successful() {
    local sample="$1"
    echo "${sample}" >> "$SUCCESSFUL_SAMPLES"
}

mark_skipped() {
    local sample="$1"
    echo "${sample}" >> "$SKIPPED_SAMPLES"
}

################################################################################
# Argument Parsing
################################################################################

show_help() {
    cat << 'EOF'
FASTP Quality Control & Trimming Script - Docker Version

USAGE:
  Single study:
    $0 STUDY_DIR [SAMPLE_FILTER]

  Multiple studies from file:
    $0 --input-file FILE [--parallel N] [--fastqc]

ARGUMENTS:
  STUDY_DIR       Name of the study directory under BASE_PATH (e.g., "Friedman_2019")
  SAMPLE_FILTER   Optional comma-separated GSM/SRX IDs to process (default: all)

OPTIONS:
  --input-file FILE   Process multiple studies from a tab-delimited file
                      (same format as geo_download_docker.sh)
  --parallel N        Process N studies simultaneously (default: 1)
  --force             Re-trim samples even if trimmed output already exists
  --fastqc            Run FastQC on trimmed files and generate MultiQC report
  --threads N         Override fastp threads per sample (default: auto-detected)
  --length N          Override minimum read length (default: 36)
  --base-path PATH    Override base data path (default: /data or $GEO_BASE_PATH)

EXAMPLES:
  # Trim all samples in a study
  bash fastp_trim_docker.sh "Friedman_2019"

  # Trim specific samples
  bash fastp_trim_docker.sh "Friedman_2019" "GSM3559136,GSM3559137"

  # Batch from input file
  bash fastp_trim_docker.sh --input-file /data/datasets.txt

  # With FastQC + MultiQC
  bash fastp_trim_docker.sh --input-file /data/datasets.txt --fastqc

  # Parallel studies
  bash fastp_trim_docker.sh --input-file /data/datasets.txt --parallel 2

DOCKER:
  docker run --rm \
    -u "$(id -u):$(id -g)" \
    -v /data:/data \
    stodo3569/fastp_tools:0.0 \
    bash /data/geo_scripts/fastp_trim_docker.sh --input-file /data/datasets.txt

FASTP PARAMETERS (from fastp_trim_1.sh):
  --trim_poly_x              Trim polyX tails
  --correction               Enable base correction for overlapping paired reads
  --detect_adapter_for_pe    Auto-detect adapters for paired-end data
  --length_required 36       Discard reads shorter than 36bp

OUTPUT STRUCTURE:
  /data/Study_Name/
    Raw_data/GSM123456/*.fastq.gz          (input - unchanged)
    Trimmed_data/GSM123456/fastp_*.fastq.gz (trimmed output)
    Trimmed_data/GSM123456/*.html           (fastp HTML reports)
    Trimmed_data/GSM123456/*.json           (fastp JSON reports)
    Trimmed_data/FastQC/                    (if --fastqc)
    Trimmed_data/MultiQC/multiqc_report.html (if --fastqc)

ERROR TRACKING:
  All errors are logged to fastp_error_log.txt per study
  Failed samples tracked in fastp_trim_failed.txt
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
                    PARALLEL_DATASETS="$2"; shift 2 ;;
                --force)
                    FORCE_RERUN=true; shift ;;
                --fastqc)
                    RUN_FASTQC=true; RUN_MULTIQC=true; shift ;;
                --threads)
                    FASTP_THREADS="$2"; shift 2 ;;
                --length)
                    FASTP_LENGTH_REQUIRED="$2"; shift 2 ;;
                --base-path)
                    BASE_PATH="$2"; shift 2 ;;
                *)
                    log_error "Unknown option: $1"; exit 1 ;;
            esac
        done

        return 0
    fi

    # Single study mode - parse all args, flags can appear anywhere
    local positional_args=()

    while [[ $# -gt 0 ]]; do
        case $1 in
            --force)
                FORCE_RERUN=true; shift ;;
            --fastqc)
                RUN_FASTQC=true; RUN_MULTIQC=true; shift ;;
            --threads)
                FASTP_THREADS="$2"; shift 2 ;;
            --length)
                FASTP_LENGTH_REQUIRED="$2"; shift 2 ;;
            --base-path)
                BASE_PATH="$2"; shift 2 ;;
            --*)
                log_error "Unknown option: $1"; exit 1 ;;
            *)
                positional_args+=("$1"); shift ;;
        esac
    done

    if [[ ${#positional_args[@]} -lt 1 ]]; then
        log_error "Missing required argument: STUDY_DIR"
        log_error "Usage: $0 STUDY_DIR [SAMPLE_FILTER] [OPTIONS]"
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

    local deps=(fastp parallel pigz)
    local missing=()

    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            missing+=("$dep")
        fi
    done

    if [[ "$RUN_FASTQC" == true ]]; then
        if ! command -v fastqc &> /dev/null; then
            log_warning "fastqc not found - skipping FastQC step"
            RUN_FASTQC=false
        fi
    fi

    if [[ "$RUN_MULTIQC" == true ]]; then
        if ! command -v multiqc &> /dev/null; then
            log_warning "multiqc not found - skipping MultiQC step"
            RUN_MULTIQC=false
        fi
    fi

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
    local study_path="$1"
    local num_samples=$2

    # Estimate: trimmed files ~80-90% of original size, plus reports
    # Very rough: ~2GB per sample for trimmed output
    local estimated_gb=$((num_samples * 2))
    local required_gb=$((estimated_gb + estimated_gb / 5))

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
# Sample Discovery
################################################################################

discover_samples() {
    local study_path="$1"
    local sample_filter="$2"
    local output_file="$3"

    local raw_data_dir="$study_path/Raw_data"

    if [[ ! -d "$raw_data_dir" ]]; then
        log_error "Raw_data directory not found: $raw_data_dir"
        log_error "Has geo_download_docker.sh been run for this study?"
        return 1
    fi

    > "$output_file"

    # Find all sample directories that contain .fastq.gz files
    local all_samples=()
    local empty_dirs=()
    while IFS= read -r sample_dir; do
        local sample_name=$(basename "$sample_dir")
        # Check that directory actually contains fastq.gz files
        local fq_count=$(find "$sample_dir" -maxdepth 1 -name "*.fastq.gz" 2>/dev/null | wc -l)
        if [[ $fq_count -gt 0 ]]; then
            all_samples+=("$sample_name")
        else
            empty_dirs+=("$sample_name")
        fi
    done < <(find "$raw_data_dir" -mindepth 1 -maxdepth 1 -type d | sort)

    # Warn about empty directories
    if [[ ${#empty_dirs[@]} -gt 0 ]]; then
        log_warning "${#empty_dirs[@]} sample director(ies) contain NO .fastq.gz files:"
        for empty_dir in "${empty_dirs[@]}"; do
            log_warning "  ⊘ $empty_dir (empty or missing FASTQ files)"
        done
        log_warning "These may indicate failed downloads - check geo_download logs"
    fi

    if [[ ${#all_samples[@]} -eq 0 ]]; then
        log_error "No samples with .fastq.gz files found in $raw_data_dir"
        return 1
    fi

    log_message "Found ${#all_samples[@]} sample(s) with FASTQ files in Raw_data/"

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
                log_warning "Filtered sample not found in Raw_data: $filter_sample"
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
# Skip Already-Trimmed Samples
################################################################################

filter_already_trimmed() {
    local study_path="$1"
    local samples_file="$2"

    # If --force is set, skip this check entirely
    if [[ "$FORCE_RERUN" == true ]]; then
        log_message "Force mode enabled - all samples will be (re)trimmed"
        return 0
    fi

    local trimmed_dir="$study_path/Trimmed_data"
    local filtered_file="${samples_file}.filtered"
    > "$filtered_file"

    local total=0
    local skipped=0
    local to_process=0

    while IFS= read -r sample; do
        [[ -z "$sample" ]] && continue
        total=$((total + 1))

        local sample_trim_dir="$trimmed_dir/$sample"

        # Check if trimmed output already exists (fastp_*.fastq.gz files)
        local trimmed_count=$(find "$sample_trim_dir" -maxdepth 1 -name "fastp_*.fastq.gz" 2>/dev/null | wc -l)

        if [[ $trimmed_count -gt 0 ]]; then
            skipped=$((skipped + 1))
            log_message "  ⊘ Skipping $sample (already trimmed: $trimmed_count output files found)"
            mark_skipped "$sample"
        else
            echo "$sample" >> "$filtered_file"
            to_process=$((to_process + 1))
        fi
    done < "$samples_file"

    # Replace the samples file with the filtered version
    mv "$filtered_file" "$samples_file"

    if [[ $skipped -gt 0 ]]; then
        log_message ""
        log_message "Skip Summary:"
        log_message "  Total samples found: $total"
        log_message "  Already trimmed (skipped): $skipped"
        log_message "  To process: $to_process"
        log_message ""
        log_message "  Use --force to re-trim all samples"
    fi

    if [[ $to_process -eq 0 ]]; then
        log_success "All $total samples already trimmed - nothing to do"
        log_message "Use --force to re-trim"
    fi

    return 0
}

################################################################################
# Paired-End / Single-End Detection
################################################################################

detect_read_layout() {
    local sample_dir="$1"

    # Check for paired-end indicators: *_1.fastq.gz and *_2.fastq.gz
    local r1_count=$(find "$sample_dir" -maxdepth 1 -name "*_1.fastq.gz" 2>/dev/null | wc -l)
    local r2_count=$(find "$sample_dir" -maxdepth 1 -name "*_2.fastq.gz" 2>/dev/null | wc -l)

    if [[ $r1_count -gt 0 ]] && [[ $r2_count -gt 0 ]]; then
        echo "paired"
    else
        # Check for any .fastq.gz files (single-end or unpaired)
        local total_count=$(find "$sample_dir" -maxdepth 1 -name "*.fastq.gz" 2>/dev/null | wc -l)
        if [[ $total_count -gt 0 ]]; then
            echo "single"
        else
            echo "none"
        fi
    fi
}

################################################################################
# FASTP Trimming - Core Processing
################################################################################

trim_single_sample() {
    local sample_name="$1"
    local raw_dir="$2"
    local trimmed_dir="$3"
    local fastp_threads="$4"
    local length_required="$5"
    local extra_args="$6"
    local error_log="$7"

    local sample_raw_dir="$raw_dir/$sample_name"
    local sample_trim_dir="$trimmed_dir/$sample_name"

    mkdir -p "$sample_trim_dir"

    # Detect read layout
    local layout=$(detect_read_layout "$sample_raw_dir")

    if [[ "$layout" == "paired" ]]; then
        # Process each paired-end read set
        local success=true

        for r1_file in "$sample_raw_dir"/*_1.fastq.gz; do
            [[ ! -f "$r1_file" ]] && continue

            local r1_basename=$(basename "$r1_file")
            local run_prefix=${r1_basename%_1.fastq.gz}

            local r2_file="$sample_raw_dir/${run_prefix}_2.fastq.gz"

            if [[ ! -f "$r2_file" ]]; then
                echo "WARNING: Missing R2 for $run_prefix, treating R1 as single-end" >> "$error_log"
                # Process as single-end
                if ! fastp \
                    -i "$r1_file" \
                    -o "$sample_trim_dir/fastp_${run_prefix}_1.fastq.gz" \
                    -h "$sample_trim_dir/${run_prefix}.html" \
                    -j "$sample_trim_dir/${run_prefix}.json" \
                    $extra_args \
                    --length_required "$length_required" \
                    --thread "$fastp_threads" 2>> "$error_log"; then
                    success=false
                fi
                continue
            fi

            # Paired-end fastp
            if ! fastp \
                -i "$r1_file" \
                -I "$r2_file" \
                -o "$sample_trim_dir/fastp_${run_prefix}_1.fastq.gz" \
                -O "$sample_trim_dir/fastp_${run_prefix}_2.fastq.gz" \
                -h "$sample_trim_dir/${run_prefix}.html" \
                -j "$sample_trim_dir/${run_prefix}.json" \
                $extra_args \
                --length_required "$length_required" \
                --thread "$fastp_threads" 2>> "$error_log"; then
                success=false
            fi
        done

        if [[ "$success" == true ]]; then
            return 0
        else
            return 1
        fi

    elif [[ "$layout" == "single" ]]; then
        # Process each single-end file
        local success=true

        for fq_file in "$sample_raw_dir"/*.fastq.gz; do
            [[ ! -f "$fq_file" ]] && continue

            # Skip R2 files if they happen to exist alongside non-paired files
            local fq_basename=$(basename "$fq_file")

            # Derive a clean prefix (strip .fastq.gz)
            local run_prefix=${fq_basename%.fastq.gz}

            if ! fastp \
                -i "$fq_file" \
                -o "$sample_trim_dir/fastp_${run_prefix}.fastq.gz" \
                -h "$sample_trim_dir/${run_prefix}.html" \
                -j "$sample_trim_dir/${run_prefix}.json" \
                $extra_args \
                --length_required "$length_required" \
                --thread "$fastp_threads" 2>> "$error_log"; then
                success=false
            fi
        done

        if [[ "$success" == true ]]; then
            return 0
        else
            return 1
        fi
    else
        echo "No FASTQ files found in $sample_raw_dir" >> "$error_log"
        return 1
    fi
}

################################################################################
# Parallel Trimming Orchestrator
################################################################################

run_trimming() {
    local study_path="$1"
    local samples_file="$2"

    log_message "========================================="
    log_message "STEP 1: FASTP Quality Control & Trimming"
    log_message "========================================="

    local raw_dir="$study_path/Raw_data"
    local trimmed_dir="$study_path/Trimmed_data"
    mkdir -p "$trimmed_dir"

    local total_samples=$(wc -l < "$samples_file")

    # Calculate optimal parallelism
    FASTP_THREADS=$(calculate_fastp_threads)
    local parallel_jobs=$(calculate_optimal_parallelism $total_samples)

    # Safety check: ensure we don't oversubscribe CPUs
    local total_cpu_usage=$((parallel_jobs * FASTP_THREADS))
    if [[ $total_cpu_usage -gt $CPUS ]]; then
        parallel_jobs=$((CPUS / FASTP_THREADS))
        [[ $parallel_jobs -lt 1 ]] && parallel_jobs=1
        total_cpu_usage=$((parallel_jobs * FASTP_THREADS))
    fi

    log_message "Trimming configuration:"
    log_message "  Samples to process: $total_samples"
    log_message "  Parallel jobs: $parallel_jobs"
    log_message "  fastp threads per job: $FASTP_THREADS"
    log_message "  Total CPU usage: ~${total_cpu_usage}/${CPUS} cores"
    log_message "  Min read length: $FASTP_LENGTH_REQUIRED"
    log_message "  Extra args: $FASTP_EXTRA_ARGS"
    log_message ""

    # Export functions and variables for GNU parallel
    export -f detect_read_layout trim_single_sample
    export -f log_message log_error log_step_error mark_trim_failed mark_successful
    export ERROR_LOG TRIM_FAILED SUCCESSFUL_SAMPLES
    export FASTP_THREADS FASTP_LENGTH_REQUIRED FASTP_EXTRA_ARGS

    cat "$samples_file" | \
        parallel -j $parallel_jobs --line-buffer --tagstring '[{1}]' \
        '
        sample={1}
        raw_dir="'"$raw_dir"'"
        trimmed_dir="'"$trimmed_dir"'"

        echo "[$(date +%H:%M:%S)] Starting fastp trimming..."

        # Detect layout for logging
        layout=$(detect_read_layout "$raw_dir/$sample")
        echo "[$(date +%H:%M:%S)] Layout: $layout"

        if trim_single_sample \
            "$sample" \
            "$raw_dir" \
            "$trimmed_dir" \
            "'"$FASTP_THREADS"'" \
            "'"$FASTP_LENGTH_REQUIRED"'" \
            "'"$FASTP_EXTRA_ARGS"'" \
            "'"$ERROR_LOG"'"; then

            echo "[$(date +%H:%M:%S)] ✓ Complete"
            mark_successful "$sample"
        else
            echo "[$(date +%H:%M:%S)] ✗ Failed"
            log_step_error "TRIM" "$sample" "fastp trimming failed"
            mark_trim_failed "$sample"
        fi
        '

    # Count results
    local successful=0
    local failed_count=0

    [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
    [[ -f "$TRIM_FAILED" ]] && failed_count=$(wc -l < "$TRIM_FAILED")

    log_message ""
    log_message "Trimming Summary:"
    log_message "  Successful: $successful"
    log_message "  Failed: $failed_count"

    # SMART RETRY: Retry failed samples sequentially with full resources
    if [[ $failed_count -gt 0 ]]; then
        log_message ""
        log_warning "Retrying $failed_count failed sample(s) sequentially..."
        log_message ""

        local retry_list=$(mktemp)
        cp "$TRIM_FAILED" "$retry_list"
        > "$TRIM_FAILED"

        local retry_success=0
        local retry_fail=0

        # Sequential retry with all CPUs for one sample
        local retry_threads=$CPUS
        [[ $retry_threads -gt 16 ]] && retry_threads=16  # fastp caps out around 16

        while IFS= read -r sample; do
            [[ -z "$sample" ]] && continue

            log_message "[RETRY] Processing $sample (threads: $retry_threads)..."

            if trim_single_sample \
                "$sample" \
                "$raw_dir" \
                "$trimmed_dir" \
                "$retry_threads" \
                "$FASTP_LENGTH_REQUIRED" \
                "$FASTP_EXTRA_ARGS" \
                "$ERROR_LOG"; then

                log_message "[RETRY] ✓ $sample succeeded"
                mark_successful "$sample"
                retry_success=$((retry_success + 1))
            else
                log_error "[RETRY] ✗ $sample failed again"
                mark_trim_failed "$sample"
                retry_fail=$((retry_fail + 1))
            fi
        done < "$retry_list"

        rm -f "$retry_list"

        log_message ""
        log_message "Retry Results:"
        log_message "  Succeeded on retry: $retry_success"
        log_message "  Still failed: $retry_fail"

        # Update counts
        successful=0
        failed_count=0
        [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
        [[ -f "$TRIM_FAILED" ]] && failed_count=$(wc -l < "$TRIM_FAILED")
    fi

    if [[ $failed_count -gt 0 ]]; then
        log_warning "Some samples failed trimming - see $TRIM_FAILED"
    else
        log_success "All samples trimmed successfully"
    fi
}

################################################################################
# FastQC (Optional)
################################################################################

run_fastqc() {
    local study_path="$1"
    local samples_file="$2"

    if [[ "$RUN_FASTQC" != true ]]; then
        return 0
    fi

    log_message "========================================="
    log_message "STEP 2: FastQC on Trimmed Files"
    log_message "========================================="

    local trimmed_dir="$study_path/Trimmed_data"
    local fastqc_dir="$trimmed_dir/FastQC"
    mkdir -p "$fastqc_dir"

    # Collect all trimmed FASTQ files
    local fastq_list=$(mktemp)
    find "$trimmed_dir" -name "fastp_*.fastq.gz" -type f | sort > "$fastq_list"

    local total_files=$(wc -l < "$fastq_list")
    if [[ $total_files -eq 0 ]]; then
        log_warning "No trimmed FASTQ files found for FastQC"
        rm -f "$fastq_list"
        return 0
    fi

    local fastqc_jobs=$((CPUS / 2))
    [[ $fastqc_jobs -lt 1 ]] && fastqc_jobs=1

    log_message "Running FastQC on $total_files trimmed files..."
    log_message "  Parallel jobs: $fastqc_jobs"
    log_message "  Output: $fastqc_dir"
    log_message ""

    # Export for parallel
    export -f log_step_error mark_fastqc_failed
    export ERROR_LOG FASTQC_FAILED

    cat "$fastq_list" | \
        parallel -j $fastqc_jobs --line-buffer \
        '
        file={1}
        fname=$(basename "$file")

        if fastqc "$file" --outdir "'"$fastqc_dir"'" --threads 2 --quiet 2>> "'"$ERROR_LOG"'"; then
            echo "[$(date +%H:%M:%S)] ✓ FastQC: $fname"
        else
            echo "[$(date +%H:%M:%S)] ✗ FastQC failed: $fname"
            log_step_error "FASTQC" "$fname" "FastQC analysis failed"
            mark_fastqc_failed "$fname"
        fi
        '

    rm -f "$fastq_list"

    local fastqc_failed_count=0
    [[ -f "$FASTQC_FAILED" ]] && fastqc_failed_count=$(wc -l < "$FASTQC_FAILED")

    if [[ $fastqc_failed_count -gt 0 ]]; then
        log_warning "$fastqc_failed_count FastQC analyses failed"
    else
        log_success "FastQC completed for all files"
    fi
}

################################################################################
# MultiQC (Optional)
################################################################################

run_multiqc() {
    local study_path="$1"

    if [[ "$RUN_MULTIQC" != true ]]; then
        return 0
    fi

    log_message "========================================="
    log_message "STEP 3: MultiQC Aggregation"
    log_message "========================================="

    local trimmed_dir="$study_path/Trimmed_data"
    local multiqc_dir="$trimmed_dir/MultiQC"
    mkdir -p "$multiqc_dir"

    log_message "Generating MultiQC report..."
    log_message "  Scanning: $trimmed_dir"
    log_message "  Output: $multiqc_dir"

    if multiqc "$trimmed_dir" \
        --outdir "$multiqc_dir" \
        --force \
        --no-data-dir \
        2>> "$ERROR_LOG"; then
        log_success "MultiQC report generated: $multiqc_dir/multiqc_report.html"
    else
        log_warning "MultiQC report generation failed - see $ERROR_LOG"
    fi
}

################################################################################
# Summary Generation
################################################################################

generate_summary() {
    local study_path="$1"

    log_message "Generating trimming summary..."

    local trimmed_dir="$study_path/Trimmed_data"

    local sample_count=$(find "$trimmed_dir" -mindepth 1 -maxdepth 1 -type d \
        ! -name "FastQC" ! -name "MultiQC" 2>/dev/null | wc -l)
    local fastq_count=$(find "$trimmed_dir" -name "fastp_*.fastq.gz" 2>/dev/null | wc -l)
    local report_count=$(find "$trimmed_dir" -name "*.html" ! -path "*/FastQC/*" ! -path "*/MultiQC/*" 2>/dev/null | wc -l)
    local total_size=$(du -sh "$trimmed_dir" 2>/dev/null | cut -f1)

    local trim_errors=0
    local fastqc_errors=0
    local successful=0
    local skipped=0

    [[ -f "$TRIM_FAILED" ]] && trim_errors=$(wc -l < "$TRIM_FAILED")
    [[ -f "$FASTQC_FAILED" ]] && fastqc_errors=$(wc -l < "$FASTQC_FAILED")
    [[ -f "$SUCCESSFUL_SAMPLES" ]] && successful=$(wc -l < "$SUCCESSFUL_SAMPLES")
    [[ -f "$SKIPPED_SAMPLES" ]] && skipped=$(wc -l < "$SKIPPED_SAMPLES")

    cat > "$study_path/fastp_trim_summary.txt" << EOF
FASTP Trimming Summary
================================================================================
Study: $(basename "$study_path")
Date: $(date)
System CPUs: $TOTAL_CPUS
CPUs Used: $CPUS
fastp Threads per Job: $FASTP_THREADS

Parameters:
-----------
Minimum Read Length: $FASTP_LENGTH_REQUIRED
Extra Args: $FASTP_EXTRA_ARGS
Force Re-run: $(if [[ "$FORCE_RERUN" == true ]]; then echo "Yes"; else echo "No"; fi)
FastQC: $(if [[ "$RUN_FASTQC" == true ]]; then echo "Yes"; else echo "No"; fi)
MultiQC: $(if [[ "$RUN_MULTIQC" == true ]]; then echo "Yes"; else echo "No"; fi)

Results:
--------
Samples Trimmed: $sample_count
Trimmed FASTQ Files: $fastq_count
fastp Reports: $report_count
Total Size: ${total_size:-Unknown}

Processing Status:
------------------
Successful (trimmed): $successful
Skipped (already trimmed): $skipped
Failed at Trimming: $trim_errors
Failed at FastQC: $fastqc_errors

Error Logs:
-----------
Main error log: $ERROR_LOG
Trim failures: $TRIM_FAILED
FastQC failures: $FASTQC_FAILED

Trimmed FASTQ Files:
---------------------
EOF

    find "$trimmed_dir" -name "fastp_*.fastq.gz" | sort >> "$study_path/fastp_trim_summary.txt"

    log_success "Summary saved to $study_path/fastp_trim_summary.txt"

    echo ""
    log_message "========================================="
    log_message "TRIMMING COMPLETED"
    log_message "========================================="
    log_message "Study: $(basename "$study_path")"
    log_message "Samples: $sample_count"
    log_message "Trimmed FASTQ files: $fastq_count"
    log_message "Total size: ${total_size:-Unknown}"
    log_message ""
    log_message "Processing Summary:"
    log_message "  Successful: $successful"
    log_message "  Skipped (already trimmed): $skipped"
    log_message "  Failed at trimming: $trim_errors"
    log_message "  Failed at FastQC: $fastqc_errors"
    log_message ""
    log_message "Location: $study_path/Trimmed_data/"
    log_message "Error log: $ERROR_LOG"
    log_message "========================================="

    if [[ $((trim_errors + fastqc_errors)) -gt 0 ]]; then
        log_warning "Some files failed processing - check error logs for details"
    fi
}

################################################################################
# Cleanup Handler
################################################################################

cleanup_on_exit() {
    local exit_code=$?

    rm -f /tmp/fastp_samples_*.txt 2>/dev/null

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
    log_message "Sample Filter: ${sample_filter:-all}"
    log_message "========================================="

    if [[ ! -d "$study_path" ]]; then
        log_error "Study directory not found: $study_path"
        return 1
    fi

    if [[ ! -d "$study_path/Raw_data" ]]; then
        log_error "Raw_data directory not found in $study_path"
        log_error "Has geo_download_docker.sh been run for this study?"
        return 1
    fi

    detect_system_resources
    check_dependencies

    # Initialize error tracking inside the study directory
    initialize_error_tracking "$study_path"

    # Discover samples
    local samples_file="$study_path/fastp_samples_to_process.txt"
    if ! discover_samples "$study_path" "$sample_filter" "$samples_file"; then
        log_error "Failed to discover samples"
        return 1
    fi

    # Filter out already-trimmed samples (unless --force)
    filter_already_trimmed "$study_path" "$samples_file"

    local num_samples=$(wc -l < "$samples_file")

    # If all samples already trimmed, skip to summary
    if [[ $num_samples -eq 0 ]]; then
        rm -f "$samples_file"
        generate_summary "$study_path"
        return 0
    fi

    check_disk_space "$study_path" "$num_samples"

    # Run fastp trimming
    run_trimming "$study_path" "$samples_file"

    # Optional: FastQC
    run_fastqc "$study_path" "$samples_file"

    # Optional: MultiQC
    run_multiqc "$study_path"

    # Cleanup and summary
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

    # Extract unique study directories from the input file
    # The input file has: Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC]
    # We aggregate by Output_Directory (study), collecting all GSM samples
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
                # If any line says "all", the whole study is "all"
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
        # Sequential processing
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
        # Parallel study processing
        log_message "Starting parallel processing with $PARALLEL_DATASETS jobs..."

        export -f log_message log_error log_warning log_success log_step_error
        export -f detect_system_resources check_dependencies
        export -f initialize_error_tracking mark_trim_failed mark_fastqc_failed mark_successful mark_skipped
        export -f discover_samples filter_already_trimmed detect_read_layout trim_single_sample
        export -f calculate_optimal_parallelism calculate_fastp_threads
        export -f check_disk_space run_trimming run_fastqc run_multiqc generate_summary
        export -f process_single_study
        export BASE_PATH PARALLEL_DATASETS RUN_FASTQC RUN_MULTIQC FORCE_RERUN
        export FASTP_LENGTH_REQUIRED FASTP_EXTRA_ARGS

        cat "$aggregated" | \
            parallel --colsep '\t' -j $PARALLEL_DATASETS --line-buffer \
            '
            study_dir={1}
            gsm_list={2}

            filter=""
            [[ "$gsm_list" != "all" ]] && filter="$gsm_list"

            bash "'"$0"'" "$study_dir" "$filter" \
                $(if [[ "'"$FORCE_RERUN"'" == true ]]; then echo "--force"; fi) \
                $(if [[ "'"$RUN_FASTQC"'" == true ]]; then echo "--fastqc"; fi) \
                --threads "'"$FASTP_THREADS"'" \
                --length "'"$FASTP_LENGTH_REQUIRED"'" \
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
