#!/bin/bash

################################################################################
# FASTP Trimming Quality Control — Failed Sample Identifier
################################################################################
#
# Purpose:
#   Verifies that all expected samples from a fastp trimming run completed
#   successfully. Identifies failures at any pipeline stage (trimming, FastQC),
#   detects missing samples that were never trimmed, and flags abnormalities
#   in fastp JSON reports (abnormally low pass rate, zero output reads,
#   extreme adapter contamination, etc.). Generates a retry input file for
#   re-trimming only the failed/problematic samples.
#
# Usage:
#   # Cross-reference against an input file (recommended)
#   bash fastp_identify_fails.sh \
#     --input-file /data/<chapter>/config/datasets.txt
#
#   # Scan /data for failures without an input file (legacy mode)
#   bash fastp_identify_fails.sh
#
# Modes:
#   With --input-file:
#     Reads the expected datasets and samples from the input file, then
#     checks whether EVERY expected sample has trimmed FASTQ files on disk,
#     corresponding fastp JSON reports, and no anomalies in the reports.
#     This catches completely missing datasets and samples that were never
#     even attempted.
#
#   Without --input-file (legacy):
#     Scans /data for dataset directories containing GSE/SRP/EGAD in
#     their name, then checks Trimmed_data/ for completeness and
#     anomalies. Cannot detect datasets that failed before creating any
#     directories.
#
# Outputs (written to /data/):
#   fastp_failed_samples_report.txt  — human-readable report of all failures
#   fastp_failed_samples_retry.txt   — input file for re-running trimming
#
# Anomaly checks (per-sample from fastp JSON):
#   1. Zero output reads (trimming produced no data)
#   2. Pass filter rate < 50% (excessive read loss)
#   3. Q20 rate < 70% after trimming (poor base quality)
#   4. Q30 rate < 50% after trimming (poor base quality)
#   5. Adapter contamination > 60% (unusually high)
#   6. Output file is empty (0 bytes)
#   7. JSON report missing (fastp did not produce stats)
#   8. Trimmed FASTQ missing for a given raw input
#
# Expected directory structure (created by fastp_trim_docker.sh):
#
#   /data/
#   └── <Output_Directory>/              e.g. Friedman_2019_GSE125050
#       ├── Raw_data/
#       │   └── <GSM>/                   e.g. GSM3559136
#       │       ├── <SRR>_1.fastq.gz     e.g. SRR8571953_1.fastq.gz
#       │       └── <SRR>_2.fastq.gz     e.g. SRR8571953_2.fastq.gz
#       ├── Trimmed_data/
#       │   └── <GSM>/
#       │       ├── fastp_<SRR>_1.fastq.gz   trimmed R1
#       │       ├── fastp_<SRR>_2.fastq.gz   trimmed R2
#       │       ├── <SRR>.html               HTML report
#       │       └── <SRR>.json               JSON report (parsed for QC)
#       ├── fastp_trim_summary.txt
#       ├── fastp_error_log.txt
#       ├── fastp_trim_failed.txt
#       ├── fastp_fastqc_failed.txt
#       ├── fastp_successful_samples.txt
#       └── fastp_skipped_samples.txt
#
#   This script checks:
#     1. Does /data/<Output_Directory>/Trimmed_data exist?
#     2. Does each expected GSM have a trimmed subdirectory?
#     3. For every raw FASTQ, is there a corresponding fastp_* output?
#     4. Does each trimmed file have a matching JSON report?
#     5. Are there anomalies inside the JSON reports?
#     6. Are there entries in the *_failed.txt error logs?
#     7. Are any trimmed FASTQ files empty (0 bytes)?
#
# Dependencies:
#   Standard bash utilities (awk, grep, find) — no Docker required
#   python3 (for JSON parsing) or awk-based fallback
#
################################################################################

set -u
set -o pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

BASE_PATH="/data"
OUTPUT_REPORT="fastp_failed_samples_report.txt"
OUTPUT_RETRY_FILE="fastp_failed_samples_retry.txt"
INPUT_FILE=""

# Anomaly thresholds (adjustable)
MIN_PASS_RATE=50          # minimum % of reads passing filter
MIN_Q20_RATE=70           # minimum % bases ≥ Q20 after trimming
MIN_Q30_RATE=50           # minimum % bases ≥ Q30 after trimming
MAX_ADAPTER_RATE=60       # maximum % adapter-contaminated reads

# ─────────────────────────────────────────────────────────────────────────────
# Color codes for terminal output
# ─────────────────────────────────────────────────────────────────────────────

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

################################################################################
# Logging Functions
################################################################################

log_info()    { echo -e "${BLUE}[INFO]${NC} $*"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $*"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $*"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $*"; }

################################################################################
# Argument Parsing
################################################################################

show_help() {
    cat << 'EOF'
FASTP Trimming QC — Failed Sample Identifier

USAGE:
  With input file (recommended):
    bash fastp_identify_fails.sh --input-file datasets.txt

  Without input file (legacy scan mode):
    bash fastp_identify_fails.sh

OPTIONS:
  --input-file FILE   Cross-reference expected samples from this file
                      against what was actually trimmed. Same format as
                      geo_download_docker.sh / fastp_trim_docker.sh input.
  --base-path PATH    Override base data path (default: /data)
  --min-pass-rate N   Minimum pass-filter rate % (default: 50)
  --min-q20 N         Minimum Q20 rate % after trimming (default: 70)
  --min-q30 N         Minimum Q30 rate % after trimming (default: 50)
  --max-adapter N     Maximum adapter contamination % (default: 60)
  --help, -h          Show this help message

INPUT FILE FORMAT (tab-separated):
  Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC_file]

  Examples:
    Friedman_2019_GSE125050	GSE125050	all
    Matlik_2023_GSE227729	GSE227729	GSM7106344,GSM7106387

WHAT THIS SCRIPT CHECKS:
  ✓ Trimmed_data directory exists per study
  ✓ Every expected GSM has a trimmed subdirectory
  ✓ Every raw FASTQ has a matching fastp_* trimmed output
  ✓ No trimmed FASTQ files are empty (0 bytes)
  ✓ fastp JSON reports exist for every run
  ✓ JSON reports show no anomalies:
      - Zero output reads
      - Pass rate < threshold
      - Q20 / Q30 below threshold
      - Adapter contamination above threshold
  ✓ fastp_trim_failed.txt entries are flagged
  ✓ fastp_error_log.txt is inspected for errors

EOF
}

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --input-file)
                [[ $# -lt 2 ]] && { log_error "Option --input-file requires a filename"; exit 1; }
                INPUT_FILE="$2"; shift 2 ;;
            --base-path)
                [[ $# -lt 2 ]] && { log_error "Option --base-path requires a path"; exit 1; }
                BASE_PATH="$2"; shift 2 ;;
            --min-pass-rate)
                MIN_PASS_RATE="$2"; shift 2 ;;
            --min-q20)
                MIN_Q20_RATE="$2"; shift 2 ;;
            --min-q30)
                MIN_Q30_RATE="$2"; shift 2 ;;
            --max-adapter)
                MAX_ADAPTER_RATE="$2"; shift 2 ;;
            --help|-h)
                show_help; exit 0 ;;
            *)
                log_error "Unknown option: $1"; show_help; exit 1 ;;
        esac
    done
}

################################################################################
# JSON Parsing Utility
#   Uses python3 if available, otherwise falls back to awk/grep.
#   Extracts key QC metrics from a fastp JSON report.
#
#   Output format (pipe-delimited):
#     total_reads|passed_reads|pass_rate|q20_rate|q30_rate|adapter_rate
#
#   All rates are percentages (0–100). If a field is missing the value
#   is set to "NA".
################################################################################

parse_fastp_json() {
    local json_file="$1"

    if [[ ! -f "$json_file" ]] || [[ ! -s "$json_file" ]]; then
        echo "NA|NA|NA|NA|NA|NA"
        return 1
    fi

    # ── Python path (prefer python3) ──
    if command -v python3 &>/dev/null; then
        python3 << PYEOF
import json, sys

try:
    with open("$json_file") as f:
        data = json.load(f)
except Exception:
    print("NA|NA|NA|NA|NA|NA")
    sys.exit(1)

# ── Total reads before filtering ──
total = data.get("summary", {}).get("before_filtering", {}).get("total_reads", "NA")

# ── Total reads after filtering ──
passed = data.get("summary", {}).get("after_filtering", {}).get("total_reads", "NA")

# ── Pass rate ──
if total not in ("NA", 0, None) and passed not in ("NA", None):
    pass_rate = round(passed / total * 100, 2)
else:
    pass_rate = "NA"

# ── Q20 and Q30 rates (after filtering) ──
q20 = data.get("summary", {}).get("after_filtering", {}).get("q20_rate", "NA")
q30 = data.get("summary", {}).get("after_filtering", {}).get("q30_rate", "NA")

# Convert fraction → percentage if needed (fastp reports as 0–1 fraction)
if isinstance(q20, (int, float)) and q20 <= 1:
    q20 = round(q20 * 100, 2)
if isinstance(q30, (int, float)) and q30 <= 1:
    q30 = round(q30 * 100, 2)

# ── Adapter trimming rate ──
adapter = data.get("adapter_cutting", {}).get("adapter_trimmed_reads", "NA")
if adapter not in ("NA", None) and total not in ("NA", 0, None):
    adapter_rate = round(adapter / total * 100, 2)
else:
    adapter_rate = "NA"

print(f"{total}|{passed}|{pass_rate}|{q20}|{q30}|{adapter_rate}")
PYEOF
    else
        # ── Fallback: awk-based extraction (less robust but functional) ──
        # Extracts values with simple pattern matching on the JSON text
        local total passed q20 q30 adapter

        total=$(grep -oP '"total_reads"\s*:\s*\K[0-9]+' "$json_file" | head -1)
        passed=$(grep -oP '"total_reads"\s*:\s*\K[0-9]+' "$json_file" | tail -1)

        # If before/after are the same line count, we only got one value
        local read_count=$(grep -c '"total_reads"' "$json_file" 2>/dev/null || echo 0)
        if [[ $read_count -lt 2 ]]; then
            passed="NA"
        fi

        q20=$(grep -oP '"q20_rate"\s*:\s*\K[0-9.]+' "$json_file" | tail -1)
        q30=$(grep -oP '"q30_rate"\s*:\s*\K[0-9.]+' "$json_file" | tail -1)
        adapter=$(grep -oP '"adapter_trimmed_reads"\s*:\s*\K[0-9]+' "$json_file" | head -1)

        total=${total:-NA}
        passed=${passed:-NA}
        q20=${q20:-NA}
        q30=${q30:-NA}
        adapter=${adapter:-NA}

        # Calculate rates
        local pass_rate="NA"
        if [[ "$total" != "NA" ]] && [[ "$passed" != "NA" ]] && [[ "$total" -gt 0 ]]; then
            pass_rate=$(awk "BEGIN {printf \"%.2f\", ($passed/$total)*100}")
        fi

        local adapter_rate="NA"
        if [[ "$total" != "NA" ]] && [[ "$adapter" != "NA" ]] && [[ "$total" -gt 0 ]]; then
            adapter_rate=$(awk "BEGIN {printf \"%.2f\", ($adapter/$total)*100}")
        fi

        # Convert q20/q30 from fraction to percentage if needed
        if [[ "$q20" != "NA" ]]; then
            local q20_is_fraction=$(awk "BEGIN {print ($q20 <= 1) ? 1 : 0}")
            [[ "$q20_is_fraction" == "1" ]] && q20=$(awk "BEGIN {printf \"%.2f\", $q20 * 100}")
        fi
        if [[ "$q30" != "NA" ]]; then
            local q30_is_fraction=$(awk "BEGIN {print ($q30 <= 1) ? 1 : 0}")
            [[ "$q30_is_fraction" == "1" ]] && q30=$(awk "BEGIN {printf \"%.2f\", $q30 * 100}")
        fi

        echo "${total}|${passed}|${pass_rate}|${q20}|${q30}|${adapter_rate}"
    fi
}

################################################################################
# Check Anomalies in fastp JSON
#   Compares extracted metrics against thresholds.
#   Returns a comma-separated list of anomaly flags, or "OK" if clean.
################################################################################

check_json_anomalies() {
    local json_metrics="$1"   # pipe-delimited: total|passed|pass_rate|q20|q30|adapter

    local total passed pass_rate q20 q30 adapter_rate
    IFS='|' read -r total passed pass_rate q20 q30 adapter_rate <<< "$json_metrics"

    local anomalies=()

    # 1. JSON couldn't be parsed at all
    if [[ "$total" == "NA" ]]; then
        echo "json_parse_error"
        return 1
    fi

    # 2. Zero reads after filtering
    if [[ "$passed" != "NA" ]] && [[ "$passed" -eq 0 ]]; then
        anomalies+=("zero_output_reads")
    fi

    # 3. Pass rate below threshold
    if [[ "$pass_rate" != "NA" ]]; then
        local below=$(awk "BEGIN {print ($pass_rate < $MIN_PASS_RATE) ? 1 : 0}")
        [[ "$below" == "1" ]] && anomalies+=("low_pass_rate(${pass_rate}%)")
    fi

    # 4. Q20 below threshold
    if [[ "$q20" != "NA" ]]; then
        local below=$(awk "BEGIN {print ($q20 < $MIN_Q20_RATE) ? 1 : 0}")
        [[ "$below" == "1" ]] && anomalies+=("low_q20(${q20}%)")
    fi

    # 5. Q30 below threshold
    if [[ "$q30" != "NA" ]]; then
        local below=$(awk "BEGIN {print ($q30 < $MIN_Q30_RATE) ? 1 : 0}")
        [[ "$below" == "1" ]] && anomalies+=("low_q30(${q30}%)")
    fi

    # 6. Adapter contamination above threshold
    if [[ "$adapter_rate" != "NA" ]]; then
        local above=$(awk "BEGIN {print ($adapter_rate > $MAX_ADAPTER_RATE) ? 1 : 0}")
        [[ "$above" == "1" ]] && anomalies+=("high_adapter(${adapter_rate}%)")
    fi

    if [[ ${#anomalies[@]} -eq 0 ]]; then
        echo "OK"
        return 0
    else
        local IFS=','
        echo "${anomalies[*]}"
        return 1
    fi
}

################################################################################
# Analyze a Single GSM Sample
#   Compares Raw_data/<GSM>/ files against Trimmed_data/<GSM>/ outputs.
#
#   Prints pipe-delimited failure lines to stdout (one per SRR/run):
#     <GSM>|<SRR>|<failure_reason>
#
#   Returns 0 if all runs in this GSM are clean, 1 otherwise.
################################################################################

analyze_gsm_sample() {
    local gsm="$1"
    local study_path="$2"

    local raw_dir="$study_path/Raw_data/$gsm"
    local trim_dir="$study_path/Trimmed_data/$gsm"

    local had_failure=0

    # ── Case 1: Trimmed directory doesn't exist at all ──
    if [[ ! -d "$trim_dir" ]]; then
        echo "${gsm}||trimmed_dir_missing"
        return 1
    fi

    # ── Case 2: Check each raw FASTQ run for a trimmed counterpart ──
    # Discover SRR run prefixes from Raw_data
    #   Paired-end: SRR123_1.fastq.gz + SRR123_2.fastq.gz  → prefix "SRR123"
    #   Single-end: SRR123.fastq.gz                         → prefix "SRR123"
    local -A run_prefixes

    if [[ -d "$raw_dir" ]]; then
        for fq in "$raw_dir"/*.fastq.gz; do
            [[ ! -f "$fq" ]] && continue
            local fname
            fname=$(basename "$fq")

            # Strip _1.fastq.gz, _2.fastq.gz, or .fastq.gz to get the run prefix
            local prefix
            if [[ "$fname" =~ ^(.+)_[12]\.fastq\.gz$ ]]; then
                prefix="${BASH_REMATCH[1]}"
            else
                prefix="${fname%.fastq.gz}"
            fi
            run_prefixes["$prefix"]=1
        done
    fi

    # If Raw_data doesn't exist or is empty, try to infer from Trimmed_data
    if [[ ${#run_prefixes[@]} -eq 0 ]]; then
        for fq in "$trim_dir"/fastp_*.fastq.gz; do
            [[ ! -f "$fq" ]] && continue
            local fname
            fname=$(basename "$fq")
            local prefix
            # Strip fastp_ prefix, then _1/_2 suffix
            fname="${fname#fastp_}"
            if [[ "$fname" =~ ^(.+)_[12]\.fastq\.gz$ ]]; then
                prefix="${BASH_REMATCH[1]}"
            else
                prefix="${fname%.fastq.gz}"
            fi
            run_prefixes["$prefix"]=1
        done
    fi

    # If still nothing, report as empty
    if [[ ${#run_prefixes[@]} -eq 0 ]]; then
        echo "${gsm}||no_fastq_found"
        return 1
    fi

    # ── Case 3: For each run prefix, verify trimmed outputs and JSON ──
    for prefix in "${!run_prefixes[@]}"; do

        # Check for at least one trimmed FASTQ (fastp_<SRR>*.fastq.gz)
        local trimmed_count
        trimmed_count=$(find "$trim_dir" -maxdepth 1 -name "fastp_${prefix}*.fastq.gz" 2>/dev/null | wc -l)

        if [[ $trimmed_count -eq 0 ]]; then
            echo "${gsm}|${prefix}|trimmed_fastq_missing"
            had_failure=1
            continue
        fi

        # Check for empty trimmed files (0 bytes)
        local empty_files
        empty_files=$(find "$trim_dir" -maxdepth 1 -name "fastp_${prefix}*.fastq.gz" -empty 2>/dev/null | wc -l)
        if [[ $empty_files -gt 0 ]]; then
            echo "${gsm}|${prefix}|trimmed_file_empty"
            had_failure=1
            continue
        fi

        # Check for JSON report
        local json_file="$trim_dir/${prefix}.json"
        if [[ ! -f "$json_file" ]]; then
            echo "${gsm}|${prefix}|json_report_missing"
            had_failure=1
            continue
        fi

        # Parse JSON and check for anomalies
        local metrics
        metrics=$(parse_fastp_json "$json_file")
        local anomaly_result
        anomaly_result=$(check_json_anomalies "$metrics")

        if [[ "$anomaly_result" != "OK" ]]; then
            echo "${gsm}|${prefix}|${anomaly_result}"
            had_failure=1
        fi
    done

    return $had_failure
}

################################################################################
# Dataset Detection (legacy mode — no input file)
#   Scans BASE_PATH for directories that look like dataset outputs.
################################################################################

detect_datasets() {
    log_info "Scanning for dataset directories in $BASE_PATH..." >&2

    local temp_datasets
    temp_datasets=$(mktemp)

    while IFS= read -r dir; do
        local dirname
        dirname=$(basename "$dir")

        [[ "$dir" == "$BASE_PATH" ]] && continue

        # Skip special directories
        if [[ "$dirname" == "geo_scripts" ]] || \
           [[ "$dirname" == "docker_images" ]] || \
           [[ "$dirname" == "lost+found" ]]; then
            continue
        fi

        # Only process directories that contain an accession pattern
        if [[ ! "$dirname" =~ GSE[0-9]+ ]] && \
           [[ ! "$dirname" =~ SRP[0-9]+ ]] && \
           [[ ! "$dirname" =~ EGAD[0-9]+ ]]; then
            continue
        fi

        # Check it has Trimmed_data or fastp log files
        if [[ -d "$dir/Trimmed_data" ]] || \
           [[ -f "$dir/fastp_trim_failed.txt" ]] || \
           [[ -f "$dir/fastp_error_log.txt" ]] || \
           [[ -f "$dir/fastp_trim_summary.txt" ]]; then
            echo "$dir" >> "$temp_datasets"
        fi
    done < <(find "$BASE_PATH" -maxdepth 1 -type d)

    echo "$temp_datasets"
}

################################################################################
# Extract Accession From Directory Name
################################################################################

extract_accession() {
    local dirname="$1"

    if [[ "$dirname" =~ (GSE[0-9]+) ]]; then
        echo "${BASH_REMATCH[1]}"
    elif [[ "$dirname" =~ (SRP[0-9]+) ]]; then
        echo "${BASH_REMATCH[1]}"
    elif [[ "$dirname" =~ (EGAD[0-9]+) ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo "UNKNOWN"
    fi
}

################################################################################
# Analyze a Single Dataset (legacy mode)
#   Checks all GSM directories under Trimmed_data/ for completeness and
#   anomalies.
#
#   Stdout output:
#     Line 1: <dataset_name>|<accession>|<failed_count>|<success_count>
#     Lines 2+: <dataset_name>|<accession>|<gsm>|<srr>|<failure_reason>
################################################################################

analyze_dataset() {
    local dataset_dir="$1"
    local dataset_name
    dataset_name=$(basename "$dataset_dir")
    local accession
    accession=$(extract_accession "$dataset_name")

    log_info "Analyzing: $dataset_name" >&2

    local failed_count=0
    local success_count=0
    local temp_failures
    temp_failures=$(mktemp)

    # ── Check fastp_trim_failed.txt for recorded failures ──
    if [[ -f "$dataset_dir/fastp_trim_failed.txt" ]] && \
       [[ -s "$dataset_dir/fastp_trim_failed.txt" ]]; then
        while IFS= read -r gsm; do
            [[ -z "$gsm" ]] && continue
            echo "${gsm}||logged_trim_failure" >> "$temp_failures"
        done < "$dataset_dir/fastp_trim_failed.txt"
    fi

    # ── Check fastp_fastqc_failed.txt for FastQC failures ──
    if [[ -f "$dataset_dir/fastp_fastqc_failed.txt" ]] && \
       [[ -s "$dataset_dir/fastp_fastqc_failed.txt" ]]; then
        while IFS= read -r fname; do
            [[ -z "$fname" ]] && continue
            echo "FASTQC|${fname}|fastqc_failed" >> "$temp_failures"
        done < "$dataset_dir/fastp_fastqc_failed.txt"
    fi

    # ── Check every GSM in Trimmed_data ──
    if [[ -d "$dataset_dir/Trimmed_data" ]]; then
        while IFS= read -r gsm_dir; do
            local gsm
            gsm=$(basename "$gsm_dir")

            # Skip FastQC/MultiQC subdirectories
            [[ "$gsm" == "FastQC" ]] && continue
            [[ "$gsm" == "MultiQC" ]] && continue

            local gsm_output
            gsm_output=$(analyze_gsm_sample "$gsm" "$dataset_dir")

            if [[ -n "$gsm_output" ]]; then
                # Has failures — check if they're real failures or just anomalies
                while IFS= read -r line; do
                    [[ -z "$line" ]] && continue
                    echo "$line" >> "$temp_failures"
                done <<< "$gsm_output"
            fi
        done < <(find "$dataset_dir/Trimmed_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    fi

    # ── Also check Raw_data for GSMs that never got trimmed ──
    if [[ -d "$dataset_dir/Raw_data" ]]; then
        while IFS= read -r gsm_dir; do
            local gsm
            gsm=$(basename "$gsm_dir")

            # Check if this GSM has a Trimmed_data counterpart
            if [[ ! -d "$dataset_dir/Trimmed_data/$gsm" ]]; then
                # Check it's not empty raw (geo download might have failed)
                local raw_fq_count
                raw_fq_count=$(find "$gsm_dir" -maxdepth 1 -name "*.fastq.gz" 2>/dev/null | wc -l)
                if [[ $raw_fq_count -gt 0 ]]; then
                    # Has raw data but no trimming output → missed by fastp
                    if ! grep -q "^${gsm}|" "$temp_failures" 2>/dev/null; then
                        echo "${gsm}||not_trimmed" >> "$temp_failures"
                    fi
                fi
            fi
        done < <(find "$dataset_dir/Raw_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    fi

    # ── Deduplicate and count ──
    local deduped
    deduped=$(mktemp)
    sort -u "$temp_failures" > "$deduped" 2>/dev/null
    failed_count=$(wc -l < "$deduped" 2>/dev/null || echo 0)

    # Count successful GSMs (in Trimmed_data with at least one fastp_*.fastq.gz)
    if [[ -d "$dataset_dir/Trimmed_data" ]]; then
        while IFS= read -r gsm_dir; do
            local gsm
            gsm=$(basename "$gsm_dir")
            [[ "$gsm" == "FastQC" ]] && continue
            [[ "$gsm" == "MultiQC" ]] && continue

            local fq_count
            fq_count=$(find "$gsm_dir" -maxdepth 1 -name "fastp_*.fastq.gz" 2>/dev/null | wc -l)
            if [[ $fq_count -gt 0 ]]; then
                # Only count as success if not in the failure list
                if ! grep -q "^${gsm}|" "$deduped" 2>/dev/null; then
                    ((success_count++))
                fi
            fi
        done < <(find "$dataset_dir/Trimmed_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null)
    fi

    # ── Print summary line ──
    echo "${dataset_name}|${accession}|${failed_count}|${success_count}"

    # ── Print failure detail lines ──
    if [[ -s "$deduped" ]]; then
        while IFS='|' read -r gsm srr reason; do
            echo "${dataset_name}|${accession}|${gsm}|${srr}|${reason}"
        done < "$deduped"
    fi

    # ── Console summary (stderr) ──
    {
        echo ""
        echo "=================================================="
        echo "Dataset: $dataset_name"
        echo "=================================================="
        echo "Accession: ${accession}"
        echo "Successful GSMs: $success_count"
        echo "Failed / anomalous entries: $failed_count"
        echo ""
    } >&2

    rm -f "$temp_failures" "$deduped"
}

################################################################################
# ═══════════════════════════════════════════════════════════════════════════════
#   MODE 1: Input File Cross-Reference
# ═══════════════════════════════════════════════════════════════════════════════
################################################################################

check_expected_samples_from_input() {
    log_info "Cross-referencing expected samples from input file: $INPUT_FILE"
    log_info ""

    if [[ ! -f "$INPUT_FILE" ]]; then
        log_error "Input file not found: $INPUT_FILE"
        exit 1
    fi

    # ── Global counters ──
    local global_total_samples=0
    local global_successful=0
    local global_failed=0
    local global_anomalies=0
    local datasets_with_failures=0
    local total_datasets=0
    local completely_missing_datasets=0

    # Temp files for collecting results
    local temp_all_failed
    temp_all_failed=$(mktemp)
    local temp_missing_all
    temp_missing_all=$(mktemp)

    # ── Step 1: Aggregate input file by Output_Directory ──
    local aggregated
    aggregated=$(mktemp)

    awk -F'\t' '
        /^#/ { next }
        /^[[:space:]]*$/ { next }
        {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $3)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $4)

            if ($1 == "" || $2 == "" || $3 == "") next

            key = $1 "|" $2

            if (key in gsms) {
                if ($3 != "all") {
                    gsms[key] = gsms[key] "," $3
                }
            } else {
                gsms[key] = $3
                order[++n] = key
            }

            if ($4 != "") ngc[key] = $4
        }
        END {
            for (i = 1; i <= n; i++) {
                key = order[i]
                output_dir = key;  sub(/\|.*/, "", output_dir)
                accession = key;   sub(/[^|]*\|/, "", accession)
                gsm_list = gsms[key]
                ngc_file = (key in ngc) ? ngc[key] : ""
                printf "%s\t%s\t%s\t%s\n", output_dir, accession, gsm_list, ngc_file
            }
        }
    ' "$INPUT_FILE" > "$aggregated"

    total_datasets=$(wc -l < "$aggregated")
    log_info "Input file contains $total_datasets unique dataset(s)"
    log_info ""

    # ── Step 2: Initialise report file ──
    cat > "$BASE_PATH/$OUTPUT_REPORT" << EOF
================================================================================
FASTP TRIMMING — FAILED SAMPLES ANALYSIS REPORT
================================================================================
Generated: $(date)
Base Path: $BASE_PATH
Mode: Input file cross-reference ($INPUT_FILE)

Anomaly Thresholds:
  Minimum pass-filter rate: ${MIN_PASS_RATE}%
  Minimum Q20 rate:         ${MIN_Q20_RATE}%
  Minimum Q30 rate:         ${MIN_Q30_RATE}%
  Maximum adapter rate:     ${MAX_ADAPTER_RATE}%

================================================================================
PER-DATASET RESULTS
================================================================================
EOF

    # ── Step 3: Check each expected dataset ──
    while IFS=$'\t' read -r output_dir accession gsm_list ngc_file; do
        [[ -z "$output_dir" ]] && continue
        [[ -z "$accession" ]] && continue

        local dataset_dir="$BASE_PATH/$output_dir"
        local dataset_failed=0
        local dataset_anomalies=0
        local dataset_successful=0
        local dataset_total=0

        log_info "Checking: $output_dir ($accession)" >&2

        # ── Check 1: Dataset directory exists? ──
        if [[ ! -d "$dataset_dir" ]]; then
            log_error "  COMPLETELY MISSING — directory does not exist: $dataset_dir" >&2
            ((completely_missing_datasets++))
            ((datasets_with_failures++))

            echo "$output_dir|$accession|$gsm_list||dataset_missing" >> "$temp_all_failed"

            if [[ "$gsm_list" == "all" ]]; then
                echo "  - $output_dir ($accession) — \"all\" samples" >> "$temp_missing_all"
            else
                local sample_count
                sample_count=$(echo "$gsm_list" | tr ',' '\n' | wc -l)
                global_failed=$((global_failed + sample_count))
                global_total_samples=$((global_total_samples + sample_count))
            fi

            cat >> "$BASE_PATH/$OUTPUT_REPORT" << EOF

--------------------------------------------------
Dataset: $output_dir ($accession)
  STATUS: COMPLETELY MISSING — no directory on disk
--------------------------------------------------
EOF
            continue
        fi

        # ── Check 2: Trimmed_data directory exists? ──
        if [[ ! -d "$dataset_dir/Trimmed_data" ]]; then
            log_error "  MISSING Trimmed_data directory in: $dataset_dir" >&2
            ((datasets_with_failures++))

            echo "$output_dir|$accession|$gsm_list||no_trimmed_data" >> "$temp_all_failed"

            if [[ "$gsm_list" == "all" ]]; then
                echo "  - $output_dir ($accession) — \"all\" samples" >> "$temp_missing_all"
            else
                local sample_count
                sample_count=$(echo "$gsm_list" | tr ',' '\n' | wc -l)
                global_failed=$((global_failed + sample_count))
                global_total_samples=$((global_total_samples + sample_count))
            fi

            cat >> "$BASE_PATH/$OUTPUT_REPORT" << EOF

--------------------------------------------------
Dataset: $output_dir ($accession)
  STATUS: Trimmed_data directory MISSING (fastp not run?)
--------------------------------------------------
EOF
            continue
        fi

        # ── Check 3: Inspect fastp error log files ──
        local logged_trim_fails=0
        if [[ -f "$dataset_dir/fastp_trim_failed.txt" ]] && \
           [[ -s "$dataset_dir/fastp_trim_failed.txt" ]]; then
            logged_trim_fails=$(wc -l < "$dataset_dir/fastp_trim_failed.txt")
        fi

        local logged_fastqc_fails=0
        if [[ -f "$dataset_dir/fastp_fastqc_failed.txt" ]] && \
           [[ -s "$dataset_dir/fastp_fastqc_failed.txt" ]]; then
            logged_fastqc_fails=$(wc -l < "$dataset_dir/fastp_fastqc_failed.txt")
        fi

        # ── Check 4: Per-sample verification ──
        local build_gsm_list=()

        if [[ "$gsm_list" != "all" ]]; then
            # Specific GSMs listed → check each one
            IFS=',' read -ra build_gsm_list <<< "$gsm_list"
        else
            # "all" → check whatever exists in Raw_data + Trimmed_data
            if [[ -d "$dataset_dir/Raw_data" ]]; then
                while IFS= read -r d; do
                    build_gsm_list+=("$(basename "$d")")
                done < <(find "$dataset_dir/Raw_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
            fi
            # Also pick up anything in Trimmed_data not in Raw_data
            if [[ -d "$dataset_dir/Trimmed_data" ]]; then
                while IFS= read -r d; do
                    local g
                    g=$(basename "$d")
                    [[ "$g" == "FastQC" ]] && continue
                    [[ "$g" == "MultiQC" ]] && continue
                    # Add if not already present
                    local found=false
                    for existing in "${build_gsm_list[@]+"${build_gsm_list[@]}"}"; do
                        [[ "$existing" == "$g" ]] && { found=true; break; }
                    done
                    [[ "$found" == false ]] && build_gsm_list+=("$g")
                done < <(find "$dataset_dir/Trimmed_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
            fi
        fi

        for gsm in "${build_gsm_list[@]+"${build_gsm_list[@]}"}"; do
            gsm=$(echo "$gsm" | tr -d '[:space:]')
            [[ -z "$gsm" ]] && continue
            ((dataset_total++))
            ((global_total_samples++))

            # Does the raw sample exist?
            local raw_gsm_dir="$dataset_dir/Raw_data/$gsm"
            local has_raw=true
            if [[ ! -d "$raw_gsm_dir" ]]; then
                has_raw=false
            else
                local raw_fq_count
                raw_fq_count=$(find "$raw_gsm_dir" -maxdepth 1 -name "*.fastq.gz" 2>/dev/null | wc -l)
                [[ $raw_fq_count -eq 0 ]] && has_raw=false
            fi

            # Run the per-GSM analysis
            local gsm_output
            gsm_output=$(analyze_gsm_sample "$gsm" "$dataset_dir" 2>/dev/null)
            local gsm_exit=$?

            if [[ $gsm_exit -ne 0 ]] || [[ -n "$gsm_output" ]]; then
                # Parse each failure line
                while IFS='|' read -r g srr reason; do
                    [[ -z "$g" ]] && continue
                    echo "$output_dir|$accession|$g|$srr|$reason" >> "$temp_all_failed"
                    ((dataset_failed++))
                    ((global_failed++))

                    # Distinguish hard failures from anomalies (quality warnings)
                    case "$reason" in
                        low_pass_rate*|low_q20*|low_q30*|high_adapter*)
                            ((dataset_anomalies++))
                            ((global_anomalies++))
                            log_warning "  ANOMALY: $gsm ($srr) — $reason" >&2
                            ;;
                        *)
                            log_warning "  FAILED: $gsm ($srr) — $reason" >&2
                            ;;
                    esac
                done <<< "$gsm_output"
            else
                ((dataset_successful++))
                ((global_successful++))
            fi
        done

        # ── Also flag logged trim failures not already caught ──
        if [[ -f "$dataset_dir/fastp_trim_failed.txt" ]] && \
           [[ -s "$dataset_dir/fastp_trim_failed.txt" ]]; then
            while IFS= read -r gsm; do
                [[ -z "$gsm" ]] && continue
                if ! grep -q "^${output_dir}|${accession}|${gsm}|" "$temp_all_failed" 2>/dev/null; then
                    echo "$output_dir|$accession|$gsm||logged_trim_failure" >> "$temp_all_failed"
                    ((dataset_failed++))
                    ((global_failed++))
                fi
            done < "$dataset_dir/fastp_trim_failed.txt"
        fi

        # ── Per-dataset report section ──
        {
            echo ""
            echo "--------------------------------------------------"
            echo "Dataset: $output_dir ($accession)"
            echo "--------------------------------------------------"
            echo "  Expected samples:      $dataset_total"
            echo "  Successful:            $dataset_successful"
            echo "  Failed (hard):         $((dataset_failed - dataset_anomalies))"
            echo "  Anomalies (quality):   $dataset_anomalies"
            echo "  Logged trim failures:  $logged_trim_fails"
            echo "  Logged FastQC failures:$logged_fastqc_fails"
        } >> "$BASE_PATH/$OUTPUT_REPORT"

        if [[ $dataset_failed -gt 0 ]]; then
            ((datasets_with_failures++))
        fi

        # Console summary
        if [[ $dataset_failed -eq 0 ]]; then
            log_success "  $output_dir — all $dataset_total sample(s) OK" >&2
        else
            log_warning "  $output_dir — $dataset_failed failure(s) / anomalies out of $dataset_total" >&2
        fi

    done < "$aggregated"

    # ── Step 4: Global summary in report ──
    cat >> "$BASE_PATH/$OUTPUT_REPORT" << EOF

================================================================================
GLOBAL SUMMARY
================================================================================
Total datasets checked:     $total_datasets
Datasets completely missing:$completely_missing_datasets
Datasets with failures:     $datasets_with_failures
Datasets fully successful:  $((total_datasets - datasets_with_failures))

Total samples checked:      $global_total_samples
Successfully trimmed:       $global_successful
Failed (hard errors):       $((global_failed - global_anomalies))
Anomalies (quality flags):  $global_anomalies
Total problematic:          $global_failed
EOF

    if [[ $global_total_samples -gt 0 ]]; then
        echo "Success rate: $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%" \
            >> "$BASE_PATH/$OUTPUT_REPORT"
    else
        echo "Success rate: N/A (no samples checked)" >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # ── Datasets where "all" was specified but entire directory is missing ──
    if [[ -s "$temp_missing_all" ]]; then
        cat >> "$BASE_PATH/$OUTPUT_REPORT" << EOF

NOTE: The following datasets specified "all" samples but are completely
missing. Their sample count could not be determined:
EOF
        cat "$temp_missing_all" >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # ── Step 5: Failed samples detail ──
    cat >> "$BASE_PATH/$OUTPUT_REPORT" << 'EOF'

================================================================================
FAILED / ANOMALOUS SAMPLES DETAIL
================================================================================

Format: Dataset | Accession | GSM | SRR_prefix | Failure_Reason

Failure reasons:
  dataset_missing       — entire study directory absent
  no_trimmed_data       — Trimmed_data/ directory absent (fastp not run?)
  trimmed_dir_missing   — GSM subdirectory absent in Trimmed_data/
  not_trimmed           — raw FASTQ exists but no trimmed output found
  trimmed_fastq_missing — no fastp_*.fastq.gz for this SRR run
  trimmed_file_empty    — trimmed FASTQ is 0 bytes
  json_report_missing   — fastp JSON report not found for this run
  json_parse_error      — JSON file could not be parsed
  zero_output_reads     — 0 reads survived filtering
  low_pass_rate(N%)     — read pass rate below threshold
  low_q20(N%)           — Q20 rate after trimming below threshold
  low_q30(N%)           — Q30 rate after trimming below threshold
  high_adapter(N%)      — adapter contamination above threshold
  logged_trim_failure   — sample appears in fastp_trim_failed.txt

EOF

    if [[ -s "$temp_all_failed" ]]; then
        sort -t'|' -k1,1 -k3,3 "$temp_all_failed" >> "$BASE_PATH/$OUTPUT_REPORT"
    else
        echo "No failed or anomalous samples found! All trimming was successful." \
            >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # ── Step 6: Generate retry input file ──
    local retry_path="$BASE_PATH/$OUTPUT_RETRY_FILE"

    if [[ -s "$temp_all_failed" ]]; then
        cat > "$retry_path" << 'EOF'
# FASTP Trimming Retry Input File
# Generated automatically to retry only failed/anomalous samples
#
# Format: Output_Directory<TAB>GSE_accession<TAB>GSM_samples
#
# Run with:
#   bash fastp_trim_docker.sh --input-file fastp_failed_samples_retry.txt --force
#
EOF

        # Aggregate failed GSMs per dataset
        local temp_agg
        temp_agg=$(mktemp)

        awk -F'|' '
        {
            key = $1 "|" $2
            gsm = $3

            if (!(key in seen)) {
                keys[++n] = key
                seen[key] = 1
            }

            if (gsm != "" && !(key SUBSEP gsm in gsm_seen)) {
                if (key in gsm_list) {
                    gsm_list[key] = gsm_list[key] "," gsm
                } else {
                    gsm_list[key] = gsm
                }
                gsm_seen[key, gsm] = 1
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                key = keys[i]
                split(key, parts, "|")
                gsms = (key in gsm_list) ? gsm_list[key] : "all"
                printf "%s\t%s\t%s\n", parts[1], parts[2], gsms
            }
        }
        ' "$temp_all_failed" > "$temp_agg"

        cat "$temp_agg" >> "$retry_path"
        rm -f "$temp_agg"

        log_success "Generated retry file: $retry_path"
    else
        cat > "$retry_path" << 'EOF'
# No failed samples found — all trimming was successful!
EOF
        log_success "No failed samples — all trimming was successful!"
    fi

    # ── Cleanup ──
    rm -f "$temp_all_failed" "$temp_missing_all" "$aggregated"

    # ── Final console summary ──
    echo ""
    echo "=================================================="
    log_info "FASTP QC Analysis Complete"
    echo "=================================================="
    echo ""
    echo "   Total samples:      $global_total_samples"
    echo "   Successful:         $global_successful"
    echo "   Failed (hard):      $((global_failed - global_anomalies))"
    echo "   Anomalies (quality):$global_anomalies"
    if [[ $global_total_samples -gt 0 ]]; then
        echo "   Success rate:       $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%"
    fi
    echo ""
    echo "   Report:     $BASE_PATH/$OUTPUT_REPORT"
    echo "   Retry file: $retry_path"
    echo ""

    if [[ $global_failed -gt 0 ]]; then
        echo "   Next steps:"
        echo "   1. Review $BASE_PATH/$OUTPUT_REPORT for failure details"
        echo "   2. Re-run trimming with:"
        echo "      docker run --rm -u \"\$(id -u):\$(id -g)\" \\"
        echo "        -v /data:/data \\"
        echo "        stodo3569/fastp_tools:0.0 \\"
        echo "        bash /data/<scripts_dir>/fastp_trim_docker.sh \\"
        echo "        --input-file /data/$OUTPUT_RETRY_FILE --force"
        echo ""
    else
        log_success "All samples completed successfully — no retry needed!"
    fi
}

################################################################################
# ═══════════════════════════════════════════════════════════════════════════════
#   MODE 2: Legacy Scan (no input file)
# ═══════════════════════════════════════════════════════════════════════════════
################################################################################

generate_reports() {
    log_warning "Legacy mode — scanning for dataset directories"
    log_warning "Tip: Use --input-file for comprehensive checking"
    log_info ""

    # Initialise output files
    > "$BASE_PATH/$OUTPUT_RETRY_FILE"

    cat > "$BASE_PATH/$OUTPUT_REPORT" << EOF
================================================================================
FASTP TRIMMING — FAILED SAMPLES ANALYSIS REPORT
================================================================================
Generated: $(date)
Base Path: $BASE_PATH
Mode: Directory scan (no input file)

Anomaly Thresholds:
  Minimum pass-filter rate: ${MIN_PASS_RATE}%
  Minimum Q20 rate:         ${MIN_Q20_RATE}%
  Minimum Q30 rate:         ${MIN_Q30_RATE}%
  Maximum adapter rate:     ${MAX_ADAPTER_RATE}%

WARNING: Without an input file, this report can only check datasets that
have directories on disk. Datasets that never ran will NOT be detected.

For comprehensive checking, use: --input-file datasets.txt

================================================================================
EOF

    local datasets_file
    datasets_file=$(detect_datasets)

    if [[ ! -f "$datasets_file" ]]; then
        log_error "Failed to create datasets list file"
        echo "No datasets found" >> "$BASE_PATH/$OUTPUT_REPORT"
        return 1
    fi

    local total_datasets
    total_datasets=$(wc -l < "$datasets_file" 2>/dev/null || echo 0)

    log_info "Found $total_datasets dataset directories"

    if [[ $total_datasets -eq 0 ]]; then
        log_warning "No dataset directories found in $BASE_PATH"
        echo "No dataset directories found." >> "$BASE_PATH/$OUTPUT_REPORT"
        echo "# No datasets found — nothing to retry" > "$BASE_PATH/$OUTPUT_RETRY_FILE"
        rm -f "$datasets_file"
        return 0
    fi

    echo "" >> "$BASE_PATH/$OUTPUT_REPORT"
    echo "Total dataset directories found: $total_datasets" >> "$BASE_PATH/$OUTPUT_REPORT"
    echo "" >> "$BASE_PATH/$OUTPUT_REPORT"

    local global_total_samples=0
    local global_successful=0
    local global_failed=0
    local datasets_with_failures=0

    local temp_failed
    temp_failed=$(mktemp)

    while IFS= read -r dataset_dir; do
        [[ -z "$dataset_dir" ]] && continue

        local analysis_output
        analysis_output=$(analyze_dataset "$dataset_dir")

        # First line is the summary
        local summary
        summary=$(echo "$analysis_output" | head -1)
        local dataset_name failed_count success_count
        dataset_name=$(echo "$summary" | cut -d'|' -f1)
        failed_count=$(echo "$summary" | cut -d'|' -f3)
        success_count=$(echo "$summary" | cut -d'|' -f4)

        echo "$summary" >> "$BASE_PATH/$OUTPUT_REPORT"

        # Remaining lines are failure details
        local sample_lines
        sample_lines=$(echo "$analysis_output" | tail -n +2)
        if [[ -n "$sample_lines" ]]; then
            echo "$sample_lines" >> "$temp_failed"
            ((datasets_with_failures++))
        fi

        global_total_samples=$((global_total_samples + failed_count + success_count))
        global_successful=$((global_successful + success_count))
        global_failed=$((global_failed + failed_count))
    done < "$datasets_file"

    rm -f "$datasets_file"

    # ── Global summary ──
    cat >> "$BASE_PATH/$OUTPUT_REPORT" << EOF

================================================================================
GLOBAL SUMMARY
================================================================================
Total datasets analyzed:  $total_datasets
Datasets with failures:   $datasets_with_failures
Datasets fully successful:$((total_datasets - datasets_with_failures))

Total samples:            $global_total_samples
Successfully trimmed:     $global_successful
Failed / anomalous:       $global_failed
EOF

    if [[ $global_total_samples -gt 0 ]]; then
        echo "Success rate: $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%" \
            >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # ── Failed sample detail ──
    echo "" >> "$BASE_PATH/$OUTPUT_REPORT"
    echo "FAILED / ANOMALOUS SAMPLES DETAIL" >> "$BASE_PATH/$OUTPUT_REPORT"
    echo "Format: Dataset | Accession | GSM | SRR_prefix | Failure_Reason" >> "$BASE_PATH/$OUTPUT_REPORT"
    echo "" >> "$BASE_PATH/$OUTPUT_REPORT"

    if [[ -s "$temp_failed" ]]; then
        sort -t'|' -k1,1 -k3,3 "$temp_failed" >> "$BASE_PATH/$OUTPUT_REPORT"
    else
        echo "No failed samples found!" >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # ── Retry file ──
    if [[ -s "$temp_failed" ]]; then
        cat > "$BASE_PATH/$OUTPUT_RETRY_FILE" << 'EOF'
# FASTP Trimming Retry Input File
# Run with: bash fastp_trim_docker.sh --input-file fastp_failed_samples_retry.txt --force
#
EOF
        local temp_agg
        temp_agg=$(mktemp)

        awk -F'|' '
        {
            key = $1 "|" $2
            gsm = $3
            if (!(key in seen)) { keys[++n] = key; seen[key] = 1 }
            if (gsm != "" && !(key SUBSEP gsm in gsm_seen)) {
                gsm_list[key] = (key in gsm_list) ? gsm_list[key] "," gsm : gsm
                gsm_seen[key, gsm] = 1
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                key = keys[i]; split(key, p, "|")
                gsms = (key in gsm_list) ? gsm_list[key] : "all"
                printf "%s\t%s\t%s\n", p[1], p[2], gsms
            }
        }
        ' "$temp_failed" > "$temp_agg"

        cat "$temp_agg" >> "$BASE_PATH/$OUTPUT_RETRY_FILE"
        rm -f "$temp_agg"

        log_success "Generated retry file with $(grep -cv "^#" "$BASE_PATH/$OUTPUT_RETRY_FILE") datasets"
    else
        echo "# No failed samples — all trimming was successful!" > "$BASE_PATH/$OUTPUT_RETRY_FILE"
        log_success "No failed samples — all trimming was successful!"
    fi

    rm -f "$temp_failed"

    # Console summary
    echo ""
    echo "=================================================="
    log_info "FASTP QC Analysis Complete"
    echo "=================================================="
    echo ""
    echo "   Total samples:  $global_total_samples"
    echo "   Successful:     $global_successful"
    echo "   Failed:         $global_failed"
    if [[ $global_total_samples -gt 0 ]]; then
        echo "   Success rate:   $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%"
    fi
    echo ""
    echo "   Report:     $BASE_PATH/$OUTPUT_REPORT"
    echo "   Retry file: $BASE_PATH/$OUTPUT_RETRY_FILE"
    echo ""

    if [[ $global_failed -gt 0 ]]; then
        echo "   Next steps:"
        echo "   1. Review $BASE_PATH/$OUTPUT_REPORT for failure details"
        echo "   2. Re-run trimming with:"
        echo "      bash fastp_trim_docker.sh --input-file /data/$OUTPUT_RETRY_FILE --force"
        echo ""
    else
        log_success "All samples completed successfully — no retry needed!"
    fi
}

################################################################################
# Main Entry Point
################################################################################

main() {
    parse_arguments "$@"

    log_info "FASTP Trimming QC — Failed Sample Identifier"
    log_info "Base path: $BASE_PATH"
    log_info "Thresholds: pass≥${MIN_PASS_RATE}%  Q20≥${MIN_Q20_RATE}%  Q30≥${MIN_Q30_RATE}%  adapter≤${MAX_ADAPTER_RATE}%"
    log_info ""

    if [[ -n "$INPUT_FILE" ]]; then
        log_info "Mode: Input file cross-reference"
        log_info "Input file: $INPUT_FILE"
        check_expected_samples_from_input
    else
        log_warning "Mode: Legacy directory scan (no input file)"
        log_warning "Tip: Use --input-file for comprehensive checking"
        generate_reports
    fi

    echo ""
    log_success "QC complete!"
}

main "$@"
