#!/bin/bash

################################################################################
# GEO/SRA Download Quality Control — Failed Sample Identifier
################################################################################
#
# Purpose:
#   Verifies that all expected samples from a download run completed
#   successfully. Identifies failures at any pipeline stage (prefetch,
#   validation, extraction) and detects missing samples that never
#   appeared on disk. Generates a retry input file for re-downloading
#   only the failed/missing samples.
#
# Usage:
#   # Cross-reference against an input file (recommended)
#   bash geo_download_identify_fails.sh \
#     --input-file /data/<chapter>/config/datasets.txt
#
#   # Scan /data for failures without an input file (legacy mode)
#   bash geo_download_identify_fails.sh
#
# Modes:
#   With --input-file:
#     Reads the expected datasets and samples from the input file,
#     then checks whether EVERY expected sample has FASTQ files on disk.
#     This catches completely missing datasets and samples that were
#     never even attempted.
#
#   Without --input-file (legacy):
#     Scans /data for directories containing GSE/SRP in their name
#     and checks error logs. Cannot detect datasets that failed before
#     creating any directories.
#
# Outputs (written to /data/):
#   failed_samples_report.txt  — human-readable report of all failures
#   failed_samples_retry.txt   — input file for re-running downloads
#
# Dependencies:
#   Standard bash utilities (awk, grep, find) — no Docker required
#
# Note on /data:
#   Same mount point as the download script. All dataset directories
#   and error logs are expected to be under /data/<Output_Directory>/.
#
# Expected directory structure (created by geo_download_docker.sh):
#
#   /data/
#   ├── <Output_Directory>/              e.g. Friedman_2019_GSE125050
#   │   ├── Raw_data/
#   │   │   ├── <GSM>/                   e.g. GSM3559136
#   │   │   │   ├── <SRR>.fastq.gz       e.g. SRR8571953_1.fastq.gz
#   │   │   │   └── <SRR>.fastq.gz       e.g. SRR8571953_2.fastq.gz
#   │   │   ├── <GSM>/
#   │   │   │   └── ...
#   │   │   └── ...
#   │   ├── download_summary.txt
#   │   ├── error_log.txt
#   │   ├── prefetch_failed.txt
#   │   ├── validation_failed.txt
#   │   ├── extraction_failed.txt
#   │   └── successful_samples.txt
#   └── <Output_Directory>/
#       └── ...
#
#   This script checks:
#     1. Does /data/<Output_Directory> exist?
#     2. Does /data/<Output_Directory>/Raw_data exist?
#     3. Does /data/<Output_Directory>/Raw_data/<GSM> exist for each expected GSM?
#     4. Does each GSM directory contain at least one .fastq.gz file?
#     5. Are there any entries in the *_failed.txt error logs?
#
################################################################################

set -u
set -o pipefail

# Configuration
BASE_PATH="/data"
OUTPUT_REPORT="failed_samples_report.txt"
OUTPUT_RETRY_FILE="failed_samples_retry.txt"
INPUT_FILE=""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

################################################################################
# Logging Functions
################################################################################

log_info() {
    echo -e "${BLUE}[INFO]${NC} $*"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $*"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $*"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $*"
}

################################################################################
# Argument Parsing
################################################################################

show_help() {
    cat << EOF
GEO/SRA Download QC — Failed Sample Identifier

USAGE:
  With input file (recommended):
    $0 --input-file datasets.txt

  Without input file (legacy scan mode):
    $0

OPTIONS:
  --input-file FILE   Cross-reference expected samples from this file
                      against what was actually downloaded. Same format
                      as geo_download_docker.sh input files.
  --help, -h          Show this help message

INPUT FILE FORMAT (tab-separated):
  Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC_file]

  Examples:
    Friedman_2019	GSE125050	all
    Matlik_2023	GSE227729	GSM7106344
    Matlik_2023	GSE227729	GSM7106387

MODES:
  With --input-file:
    Checks every dataset and sample listed in the input file.
    Catches:
      - Datasets that completely failed (no directory created)
      - Samples that were never attempted
      - Samples that failed at prefetch/validation/extraction
      - Samples with missing FASTQ files (silent failures)

  Without --input-file (legacy):
    Scans /data for existing dataset directories only.
    Cannot detect datasets that failed before creating directories.

EOF
}

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --input-file)
                if [[ $# -lt 2 ]]; then
                    log_error "Option --input-file requires a filename"
                    exit 1
                fi
                INPUT_FILE="$2"
                shift 2
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done
}

################################################################################
# Dataset Detection (legacy mode — no input file)
################################################################################

detect_datasets() {
    log_info "Scanning for dataset directories in $BASE_PATH..." >&2

    local temp_datasets=$(mktemp)

    while IFS= read -r dir; do
        local dirname=$(basename "$dir")

        [[ "$dir" == "$BASE_PATH" ]] && continue

        # Skip special directories
        if [[ "$dirname" == "geo_scripts" ]] || \
           [[ "$dirname" == "docker_images" ]] || \
           [[ "$dirname" == "lost+found" ]]; then
            continue
        fi

        # Only process directories that contain GSE or SRP in their name
        if [[ ! "$dirname" =~ GSE[0-9]+ ]] && [[ ! "$dirname" =~ SRP[0-9]+ ]]; then
            continue
        fi

        # Check if it looks like a dataset directory
        if [[ -d "$dir/Raw_data" ]] || \
           [[ -f "$dir/prefetch_failed.txt" ]] || \
           [[ -f "$dir/validation_failed.txt" ]] || \
           [[ -f "$dir/extraction_failed.txt" ]] || \
           [[ -f "$dir/error_log.txt" ]]; then
            echo "$dir" >> "$temp_datasets"
        fi
    done < <(find "$BASE_PATH" -maxdepth 1 -type d)

    echo "$temp_datasets"
}

################################################################################
# Failed Sample Extraction
################################################################################

extract_gsm_from_path() {
    local raw_data_path=$1
    basename "$(dirname "$raw_data_path")"
}

extract_srr_from_path() {
    local raw_data_path=$1
    basename "$raw_data_path"
}

analyze_dataset() {
    local dataset_dir=$1
    local dataset_name=$(basename "$dataset_dir")

    log_info "Analyzing: $dataset_name" >&2

    # Try to determine GSE accession from directory name
    local gse_accession=""
    if [[ "$dataset_name" =~ (GSE[0-9]+) ]]; then
        gse_accession="${BASH_REMATCH[1]}"
    elif [[ "$dataset_name" =~ (SRP[0-9]+) ]]; then
        gse_accession="${BASH_REMATCH[1]}"
    fi

    # Initialize counters
    local total_samples=0
    local successful_samples=0
    local prefetch_failed=0
    local validation_failed=0
    local extraction_failed=0
    local missing_fastq=0

    # Track failed samples
    local temp_failed_map=$(mktemp)
    local temp_failed_reasons=$(mktemp)

    # Check error logs
    if [[ -f "$dataset_dir/prefetch_failed.txt" ]]; then
        while IFS=$'\t' read -r gsm srr || [[ -n "$gsm" ]]; do
            [[ -z "$gsm" ]] && continue
            echo "${gsm}|${srr}" >> "$temp_failed_map"
            echo "${gsm}|${srr}|prefetch" >> "$temp_failed_reasons"
            ((prefetch_failed++))
        done < "$dataset_dir/prefetch_failed.txt"
    fi

    if [[ -f "$dataset_dir/validation_failed.txt" ]]; then
        while IFS=$'\t' read -r gsm srr || [[ -n "$gsm" ]]; do
            [[ -z "$gsm" ]] && continue
            if ! grep -q "^${gsm}|${srr}$" "$temp_failed_map" 2>/dev/null; then
                echo "${gsm}|${srr}" >> "$temp_failed_map"
                echo "${gsm}|${srr}|validation" >> "$temp_failed_reasons"
                ((validation_failed++))
            fi
        done < "$dataset_dir/validation_failed.txt"
    fi

    if [[ -f "$dataset_dir/extraction_failed.txt" ]]; then
        while IFS=$'\t' read -r gsm srr || [[ -n "$gsm" ]]; do
            [[ -z "$gsm" ]] && continue
            if ! grep -q "^${gsm}|${srr}$" "$temp_failed_map" 2>/dev/null; then
                echo "${gsm}|${srr}" >> "$temp_failed_map"
                echo "${gsm}|${srr}|extraction" >> "$temp_failed_reasons"
                ((extraction_failed++))
            fi
        done < "$dataset_dir/extraction_failed.txt"
    fi

    # Check for samples with no FASTQ files (silent failures)
    if [[ -d "$dataset_dir/Raw_data" ]]; then
        while IFS= read -r gsm_dir; do
            local gsm=$(basename "$gsm_dir")
            ((total_samples++))

            local fastq_count=$(find "$gsm_dir" -name "*.fastq.gz" 2>/dev/null | wc -l)

            if [[ $fastq_count -eq 0 ]]; then
                while IFS= read -r srr_dir; do
                    [[ ! -d "$srr_dir" ]] && continue
                    local srr=$(basename "$srr_dir")

                    if ! grep -q "^${gsm}|${srr}$" "$temp_failed_map" 2>/dev/null; then
                        echo "${gsm}|${srr}" >> "$temp_failed_map"
                        echo "${gsm}|${srr}|missing_fastq" >> "$temp_failed_reasons"
                        ((missing_fastq++))
                    fi
                done < <(find "$gsm_dir" -mindepth 1 -maxdepth 1 -type d 2>/dev/null)
            else
                ((successful_samples++))
            fi
        done < <(find "$dataset_dir/Raw_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null)
    fi

    # Count successfully completed samples
    if [[ -f "$dataset_dir/successful_samples.txt" ]]; then
        local success_count=$(wc -l < "$dataset_dir/successful_samples.txt" 2>/dev/null || echo 0)
        if [[ $success_count -gt $successful_samples ]]; then
            successful_samples=$success_count
        fi
    fi

    # Calculate totals
    local total_failed=$(wc -l < "$temp_failed_map" 2>/dev/null || echo 0)
    if [[ $total_samples -eq 0 ]]; then
        total_samples=$((successful_samples + total_failed))
    fi

    # Print summary for this dataset (to stderr)
    {
        echo ""
        echo "=================================================="
        echo "Dataset: $dataset_name"
        echo "=================================================="
        echo "GSE/SRP Accession: ${gse_accession:-UNKNOWN}"
        echo "Total samples processed: $total_samples"
        echo "Successfully completed: $successful_samples"
        echo "Failed samples: $total_failed"
        if [[ $total_failed -gt 0 ]]; then
            echo "  - Prefetch failures: $prefetch_failed"
            echo "  - Validation failures: $validation_failed"
            echo "  - Extraction failures: $extraction_failed"
            echo "  - Missing FASTQ (silent failure): $missing_fastq"
        fi
        echo ""
    } >&2

    # Return dataset info (to stdout)
    echo "$dataset_name|$gse_accession|$total_failed|$successful_samples"

    # Export failed samples for this dataset (to stdout)
    if [[ -f "$temp_failed_reasons" ]]; then
        while IFS='|' read -r gsm srr reason; do
            echo "$dataset_name|$gse_accession|$gsm|$srr|$reason"
        done < "$temp_failed_reasons"
    fi

    # Cleanup temp files
    rm -f "$temp_failed_map" "$temp_failed_reasons"
}

################################################################################
# Input File Cross-Reference Mode (NEW)
################################################################################

check_expected_samples_from_input() {
    log_info "Cross-referencing expected samples from input file: $INPUT_FILE"
    log_info ""

    if [[ ! -f "$INPUT_FILE" ]]; then
        log_error "Input file not found: $INPUT_FILE"
        exit 1
    fi

    # Global counters
    local global_total_samples=0
    local global_successful=0
    local global_failed=0
    local datasets_with_failures=0
    local total_datasets=0
    local completely_missing_datasets=0

    # Temporary file for collecting all failed samples
    local temp_all_failed=$(mktemp)

    # Temporary file for tracking missing datasets where sample count is unknown ("all")
    local temp_missing_all=$(mktemp)

    # Step 1: Aggregate the input file by Output_Directory + Accession
    # (same logic as the download script)
    local aggregated=$(mktemp)

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

            if ($4 != "") {
                ngc[key] = $4
            }
        }

        END {
            for (i = 1; i <= n; i++) {
                key = order[i]
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

    total_datasets=$(wc -l < "$aggregated")
    log_info "Input file contains $total_datasets unique dataset(s)"
    log_info ""

    # Step 2: Check each expected dataset
    while IFS=$'\t' read -r output_dir accession gsm_list ngc_file; do
        [[ -z "$output_dir" ]] && continue
        [[ -z "$accession" ]] && continue

        local dataset_dir="$BASE_PATH/$output_dir"
        local dataset_failed=0
        local dataset_successful=0
        local dataset_total=0

        log_info "Checking: $output_dir ($accession)" >&2

        # Check 1: Does the dataset directory even exist?
        if [[ ! -d "$dataset_dir" ]]; then
            log_error "  COMPLETELY MISSING — directory does not exist: $dataset_dir" >&2
            ((completely_missing_datasets++))
            ((datasets_with_failures++))

            # Add the whole dataset as failed
            echo "$output_dir|$accession|$gsm_list|$ngc_file|dataset_missing" >> "$temp_all_failed"

            if [[ "$gsm_list" == "all" ]]; then
                # Can't count samples — track separately
                echo "  - $output_dir ($accession) — \"all\" samples" >> "$temp_missing_all"
            else
                local sample_count=$(echo "$gsm_list" | tr ',' '\n' | wc -l)
                global_failed=$((global_failed + sample_count))
                global_total_samples=$((global_total_samples + sample_count))
            fi

            continue
        fi

        # Check 2: Does Raw_data exist?
        if [[ ! -d "$dataset_dir/Raw_data" ]]; then
            log_error "  MISSING Raw_data directory in: $dataset_dir" >&2
            ((datasets_with_failures++))

            echo "$output_dir|$accession|$gsm_list|$ngc_file|no_raw_data" >> "$temp_all_failed"

            if [[ "$gsm_list" == "all" ]]; then
                # Can't count samples — track separately
                echo "  - $output_dir ($accession) — \"all\" samples" >> "$temp_missing_all"
            else
                local sample_count=$(echo "$gsm_list" | tr ',' '\n' | wc -l)
                global_failed=$((global_failed + sample_count))
                global_total_samples=$((global_total_samples + sample_count))
            fi

            continue
        fi

        # Check 3: Check error log files for known failures
        local logged_failures=0
        for fail_file in prefetch_failed.txt validation_failed.txt extraction_failed.txt; do
            if [[ -f "$dataset_dir/$fail_file" ]] && [[ -s "$dataset_dir/$fail_file" ]]; then
                local count=$(wc -l < "$dataset_dir/$fail_file")
                logged_failures=$((logged_failures + count))
            fi
        done

        # Check 4: If specific GSMs are listed, check each one
        if [[ "$gsm_list" != "all" ]]; then
            IFS=',' read -ra EXPECTED_GSMS <<< "$gsm_list"

            for gsm in "${EXPECTED_GSMS[@]}"; do
                gsm=$(echo "$gsm" | tr -d '[:space:]')
                [[ -z "$gsm" ]] && continue
                ((dataset_total++))
                ((global_total_samples++))

                local gsm_dir="$dataset_dir/Raw_data/$gsm"

                if [[ ! -d "$gsm_dir" ]]; then
                    # GSM directory doesn't exist at all
                    log_warning "  MISSING: $gsm — no directory found" >&2
                    echo "$output_dir|$accession|$gsm||sample_missing" >> "$temp_all_failed"
                    ((dataset_failed++))
                    ((global_failed++))
                else
                    # GSM directory exists — check for FASTQ files
                    local fastq_count=$(find "$gsm_dir" -name "*.fastq.gz" 2>/dev/null | wc -l)

                    if [[ $fastq_count -eq 0 ]]; then
                        # Check if it's a known failure
                        local failure_reason="missing_fastq"
                        if [[ -f "$dataset_dir/prefetch_failed.txt" ]] && grep -q "$gsm" "$dataset_dir/prefetch_failed.txt" 2>/dev/null; then
                            failure_reason="prefetch"
                        elif [[ -f "$dataset_dir/validation_failed.txt" ]] && grep -q "$gsm" "$dataset_dir/validation_failed.txt" 2>/dev/null; then
                            failure_reason="validation"
                        elif [[ -f "$dataset_dir/extraction_failed.txt" ]] && grep -q "$gsm" "$dataset_dir/extraction_failed.txt" 2>/dev/null; then
                            failure_reason="extraction"
                        fi

                        log_warning "  FAILED: $gsm — $failure_reason" >&2
                        echo "$output_dir|$accession|$gsm||$failure_reason" >> "$temp_all_failed"
                        ((dataset_failed++))
                        ((global_failed++))
                    else
                        ((dataset_successful++))
                        ((global_successful++))
                    fi
                fi
            done
        else
            # "all" mode — check whatever is on disk + error logs
            # Count GSM directories with FASTQ files
            if [[ -d "$dataset_dir/Raw_data" ]]; then
                while IFS= read -r gsm_dir; do
                    local gsm=$(basename "$gsm_dir")
                    ((dataset_total++))
                    ((global_total_samples++))

                    local fastq_count=$(find "$gsm_dir" -name "*.fastq.gz" 2>/dev/null | wc -l)

                    if [[ $fastq_count -eq 0 ]]; then
                        local failure_reason="missing_fastq"
                        if [[ -f "$dataset_dir/prefetch_failed.txt" ]] && grep -q "$gsm" "$dataset_dir/prefetch_failed.txt" 2>/dev/null; then
                            failure_reason="prefetch"
                        elif [[ -f "$dataset_dir/validation_failed.txt" ]] && grep -q "$gsm" "$dataset_dir/validation_failed.txt" 2>/dev/null; then
                            failure_reason="validation"
                        elif [[ -f "$dataset_dir/extraction_failed.txt" ]] && grep -q "$gsm" "$dataset_dir/extraction_failed.txt" 2>/dev/null; then
                            failure_reason="extraction"
                        fi

                        log_warning "  FAILED: $gsm — $failure_reason" >&2
                        echo "$output_dir|$accession|$gsm||$failure_reason" >> "$temp_all_failed"
                        ((dataset_failed++))
                        ((global_failed++))
                    else
                        ((dataset_successful++))
                        ((global_successful++))
                    fi
                done < <(find "$dataset_dir/Raw_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null)

                # Also check error logs for samples not even on disk
                for fail_file in prefetch_failed.txt validation_failed.txt extraction_failed.txt; do
                    if [[ -f "$dataset_dir/$fail_file" ]]; then
                        while IFS=$'\t' read -r gsm srr || [[ -n "$gsm" ]]; do
                            [[ -z "$gsm" ]] && continue
                            # Check if already counted
                            if ! grep -q "^${output_dir}|${accession}|${gsm}|" "$temp_all_failed" 2>/dev/null; then
                                local reason="${fail_file%%_failed.txt}"
                                echo "$output_dir|$accession|$gsm|$srr|$reason" >> "$temp_all_failed"
                                ((dataset_failed++))
                                ((global_failed++))
                                ((dataset_total++))
                                ((global_total_samples++))
                            fi
                        done < "$dataset_dir/$fail_file"
                    fi
                done
            fi
        fi

        # Report per dataset
        {
            echo ""
            echo "=================================================="
            echo "Dataset: $output_dir"
            echo "=================================================="
            echo "Accession: $accession"
            echo "Expected samples: $dataset_total"
            echo "Successful: $dataset_successful"
            echo "Failed/Missing: $dataset_failed"
            if [[ $dataset_failed -gt 0 ]]; then
                echo "Status: INCOMPLETE"
            else
                echo "Status: COMPLETE"
            fi
            echo ""
        } >&2

        if [[ $dataset_failed -gt 0 ]]; then
            ((datasets_with_failures++))
        fi

    done < "$aggregated"

    rm -f "$aggregated"

    # Step 3: Write the report
    > "$OUTPUT_REPORT"
    cat > "$OUTPUT_REPORT" << EOF
================================================================================
FAILED SAMPLES ANALYSIS REPORT
================================================================================
Generated: $(date)
Base Path: $BASE_PATH
Input File: $INPUT_FILE
Mode: Input file cross-reference

================================================================================
GLOBAL SUMMARY
================================================================================
Total datasets in input file: $total_datasets
Datasets fully successful: $((total_datasets - datasets_with_failures))
Datasets with failures: $datasets_with_failures
Completely missing datasets: $completely_missing_datasets

Total samples checked: $global_total_samples
Successfully completed: $global_successful
Failed/missing: $global_failed
EOF

    if [[ $global_total_samples -gt 0 ]]; then
        echo "Success rate: $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%" >> "$OUTPUT_REPORT"
    else
        echo "Success rate: N/A (no samples processed)" >> "$OUTPUT_REPORT"
    fi

    # List missing "all" datasets separately
    if [[ -s "$temp_missing_all" ]]; then
        echo "" >> "$OUTPUT_REPORT"
        echo "WARNING: The following dataset(s) are completely missing and listed" >> "$OUTPUT_REPORT"
        echo "as \"all\" samples — sample count unknown, NOT included in totals above:" >> "$OUTPUT_REPORT"
        cat "$temp_missing_all" >> "$OUTPUT_REPORT"
    fi

    cat >> "$OUTPUT_REPORT" << 'EOF'

================================================================================
FAILED SAMPLES DETAIL
================================================================================

Format: Dataset | Accession | GSM | SRR | Failure_Reason

Failure reasons:
  dataset_missing  — entire dataset directory not found on disk
  no_raw_data      — dataset directory exists but no Raw_data folder
  sample_missing   — specific GSM directory not found
  prefetch         — SRA download failed
  validation       — file integrity check failed
  extraction       — FASTQ extraction failed
  missing_fastq    — directory exists but no .fastq.gz files found

EOF

    if [[ -s "$temp_all_failed" ]]; then
        sort -t'|' -k1,1 -k3,3 "$temp_all_failed" >> "$OUTPUT_REPORT"
    else
        echo "No failed samples found! All downloads were successful." >> "$OUTPUT_REPORT"
    fi

    # Step 4: Generate retry input file
    > "$OUTPUT_RETRY_FILE"

    if [[ -s "$temp_all_failed" ]]; then
        cat > "$OUTPUT_RETRY_FILE" << 'EOF'
# GEO Download Retry Input File
# Generated automatically to retry only failed/missing samples
#
# Format: Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC_file]
#
# Run with:
#   bash geo_download_docker.sh --input-file failed_samples_retry.txt
#
EOF

        # Aggregate failed GSMs by dataset + accession
        local temp_aggregated=$(mktemp)

        awk -F'|' '
        {
            key = $1 "|" $2
            gsm = $3
            ngc = ""

            if (!(key in seen)) {
                keys[++n] = key
                seen[key] = 1
            }

            # For dataset_missing or no_raw_data, keep original gsm_list
            if ($5 == "dataset_missing" || $5 == "no_raw_data") {
                gsm_list[key] = gsm
                next
            }

            if (gsm == "" || gsm == "all") next

            if (!(key SUBSEP gsm in gsm_seen)) {
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
                printf "%s|%s|%s\n", parts[1], parts[2], gsm_list[key]
            }
        }
        ' "$temp_all_failed" > "$temp_aggregated"

        while IFS='|' read -r dataset_name gse_accession gsm_list; do
            # Try to find NGC file from original input
            local ngc_file=""
            if [[ -f "$INPUT_FILE" ]]; then
                ngc_file=$(awk -F'\t' -v dir="$dataset_name" '$1 == dir && $4 != "" {print $4; exit}' "$INPUT_FILE")
            fi

            printf "%s\t%s\t%s\t%s\n" "$dataset_name" "$gse_accession" "$gsm_list" "$ngc_file" >> "$OUTPUT_RETRY_FILE"
        done < "$temp_aggregated"

        rm -f "$temp_aggregated"

        local retry_count=$(grep -v "^#" "$OUTPUT_RETRY_FILE" | grep -v "^$" | wc -l)
        log_success "Generated retry input file with $retry_count dataset(s) to retry"
    else
        echo "# No failed samples found — all downloads were successful!" >> "$OUTPUT_RETRY_FILE"
        log_success "No failed samples found — all downloads were successful!"
    fi

    rm -f "$temp_all_failed"

    # Final summary
    echo ""
    echo "=================================================="
    log_info "QUALITY CONTROL SUMMARY"
    echo "=================================================="
    echo ""
    echo "   Input file: $INPUT_FILE"
    echo "   Total datasets: $total_datasets"
    echo "   Datasets fully successful: $((total_datasets - datasets_with_failures))"
    echo "   Datasets with failures: $datasets_with_failures"
    echo "   Completely missing: $completely_missing_datasets"
    echo ""
    echo "   Total samples (verifiable): $global_total_samples"
    echo "   Successful: $global_successful"
    echo "   Failed/Missing: $global_failed"
    if [[ $global_total_samples -gt 0 ]]; then
        echo "   Success rate: $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%"
    else
        echo "   Success rate: N/A"
    fi

    if [[ -s "$temp_missing_all" ]]; then
        echo ""
        echo "   ⚠ Missing dataset(s) with unknown sample count (not in totals above):"
        cat "$temp_missing_all"
    fi

    echo ""
    echo "   Report: $OUTPUT_REPORT"
    echo "   Retry file: $OUTPUT_RETRY_FILE"
    echo ""

    if [[ $global_failed -gt 0 ]] || [[ -s "$temp_missing_all" ]]; then
        echo "   Next steps:"
        echo "   1. Review $OUTPUT_REPORT for failure details"
        echo "   2. Retry with:"
        echo "      docker run --rm -u \"\$(id -u):\$(id -g)\" \\"
        echo "        --shm-size=60g \\"
        echo "        -v /data:/data \\"
        echo "        stodo3569/geo-tools:0.0 \\"
        echo "        bash /data/<scripts_dir>/geo_download_docker.sh \\"
        echo "        --input-file /data/$OUTPUT_RETRY_FILE"
        echo ""
    else
        log_success "All expected samples are present — no retry needed!"
    fi

    rm -f "$temp_missing_all"
}

################################################################################
# Legacy Report Generation (no input file)
################################################################################

generate_reports() {
    log_info "Generating reports (legacy scan mode — no input file)..."

    > "$OUTPUT_REPORT"
    > "$OUTPUT_RETRY_FILE"

    cat > "$OUTPUT_REPORT" << 'EOF'
================================================================================
FAILED SAMPLES ANALYSIS REPORT
================================================================================
Generated: $(date)
Base Path: $BASE_PATH
Mode: Directory scan (no input file)

WARNING: Without an input file, this report can only check datasets that
created directories on disk. Datasets that completely failed before
creating any directories will NOT be detected.

For comprehensive checking, use: --input-file datasets.txt

================================================================================
EOF

    local datasets_file=$(detect_datasets)

    if [[ ! -f "$datasets_file" ]]; then
        log_error "Failed to create datasets list file"
        echo "No datasets found" >> "$OUTPUT_REPORT"
        return 1
    fi

    local total_datasets=$(wc -l < "$datasets_file" 2>/dev/null || echo 0)

    log_info "Found $total_datasets dataset directories"

    if [[ $total_datasets -eq 0 ]]; then
        log_warning "No dataset directories found in $BASE_PATH"
        log_warning "Looking for directories with GSE[numbers] or SRP[numbers] in their name"
        log_warning ""
        log_warning "Available directories in $BASE_PATH:"
        ls -1 "$BASE_PATH" 2>/dev/null | head -20 | while read dir; do
            log_warning "  - $dir"
        done

        echo "No dataset directories found with GSE or SRP accessions." >> "$OUTPUT_REPORT"
        echo "# No datasets found — nothing to retry" > "$OUTPUT_RETRY_FILE"

        rm -f "$datasets_file"

        log_warning "No datasets to analyze — check directory names contain GSE or SRP"
        return 0
    fi

    echo "" >> "$OUTPUT_REPORT"
    echo "Total dataset directories found: $total_datasets" >> "$OUTPUT_REPORT"
    echo "" >> "$OUTPUT_REPORT"

    local global_total_samples=0
    local global_successful=0
    local global_failed=0
    local datasets_with_failures=0

    local temp_failed=$(mktemp)

    while IFS= read -r dataset_dir; do
        [[ -z "$dataset_dir" ]] && continue

        local analysis_output=$(analyze_dataset "$dataset_dir")

        local summary=$(echo "$analysis_output" | head -1)
        local dataset_name=$(echo "$summary" | cut -d'|' -f1)
        local gse_accession=$(echo "$summary" | cut -d'|' -f2)
        local failed_count=$(echo "$summary" | cut -d'|' -f3)
        local success_count=$(echo "$summary" | cut -d'|' -f4)

        echo "$analysis_output" | head -1 >> "$OUTPUT_REPORT"

        local sample_lines=$(echo "$analysis_output" | tail -n +2)
        if [[ -n "$sample_lines" ]]; then
            echo "$sample_lines" >> "$temp_failed"
            ((datasets_with_failures++))
        fi

        global_total_samples=$((global_total_samples + failed_count + success_count))
        global_successful=$((global_successful + success_count))
        global_failed=$((global_failed + failed_count))
    done < "$datasets_file"

    rm -f "$datasets_file"

    cat >> "$OUTPUT_REPORT" << EOF

================================================================================
GLOBAL SUMMARY
================================================================================
Total datasets analyzed: $total_datasets
Datasets with failures: $datasets_with_failures
Datasets fully successful: $((total_datasets - datasets_with_failures))

Total samples across all datasets: $global_total_samples
Successfully completed: $global_successful
Failed samples: $global_failed
EOF

    if [[ $global_total_samples -gt 0 ]]; then
        echo "Success rate: $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%" >> "$OUTPUT_REPORT"
    else
        echo "Success rate: N/A (no samples processed)" >> "$OUTPUT_REPORT"
    fi

    cat >> "$OUTPUT_REPORT" << 'EOF'

================================================================================
FAILED SAMPLES DETAIL
================================================================================

Format: Dataset | GSE/SRP | GSM | SRR | Failure_Reason

EOF

    if [[ -s "$temp_failed" ]]; then
        sort -t'|' -k1,1 -k3,3 "$temp_failed" >> "$OUTPUT_REPORT"
    else
        echo "No failed samples found! All downloads were successful." >> "$OUTPUT_REPORT"
    fi

    # Generate retry input file
    if [[ -s "$temp_failed" ]]; then
        cat >> "$OUTPUT_RETRY_FILE" << 'EOF'
# GEO Download Retry Input File
# Generated automatically to retry only failed samples
#
# Format: Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC_file]
#
# Run with: bash geo_download_docker.sh --input-file failed_samples_retry.txt
#
EOF

        local temp_aggregated=$(mktemp)

        awk -F'|' '
        {
            key = $1 "|" $2
            gsm = $3

            if (!(key in seen)) {
                keys[++n] = key
                seen[key] = 1
            }

            if (!(key SUBSEP gsm in gsm_seen)) {
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
                printf "%s|%s|%s\n", parts[1], parts[2], gsm_list[key]
            }
        }
        ' "$temp_failed" > "$temp_aggregated"

        while IFS='|' read -r dataset_name gse_accession gsm_list; do
            local ngc_file=""
            if [[ -f "$BASE_PATH/$dataset_name/download_summary.txt" ]]; then
                ngc_file=$(grep -oP "NGC file: \K.*" "$BASE_PATH/$dataset_name/download_summary.txt" 2>/dev/null || echo "")
            fi
            printf "%s\t%s\t%s\t%s\n" "$dataset_name" "$gse_accession" "$gsm_list" "$ngc_file" >> "$OUTPUT_RETRY_FILE"
        done < "$temp_aggregated"

        rm -f "$temp_aggregated"

        log_success "Generated retry input file with $(grep -v "^#" "$OUTPUT_RETRY_FILE" | wc -l) datasets to retry"
    else
        echo "# No failed samples found — all downloads were successful!" >> "$OUTPUT_RETRY_FILE"
        log_success "No failed samples found — all downloads were successful!"
    fi

    rm -f "$temp_failed"

    # Final summary
    echo ""
    echo "=================================================="
    log_success "Analysis complete!"
    echo "=================================================="
    echo ""
    echo "   Total samples: $global_total_samples"
    echo "   Successful: $global_successful"
    echo "   Failed: $global_failed"
    if [[ $global_total_samples -gt 0 ]]; then
        echo "   Success rate: $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%"
    else
        echo "   Success rate: N/A (no samples processed)"
    fi
    echo ""
    echo "   Report: $OUTPUT_REPORT"
    echo "   Retry file: $OUTPUT_RETRY_FILE"
    echo ""

    if [[ $global_failed -gt 0 ]]; then
        echo "   Next steps:"
        echo "   1. Review $OUTPUT_REPORT for failure details"
        echo "   2. Run retry with:"
        echo "      docker run --rm -u \"\$(id -u):\$(id -g)\" \\"
        echo "        --shm-size=60g \\"
        echo "        -v /data:/data \\"
        echo "        stodo3569/geo-tools:0.0 \\"
        echo "        bash /data/<scripts_dir>/geo_download_docker.sh \\"
        echo "        --input-file /data/$OUTPUT_RETRY_FILE"
        echo ""
    else
        log_success "All samples completed successfully — no retry needed!"
    fi
}

################################################################################
# Main
################################################################################

main() {
    parse_arguments "$@"

    log_info "GEO/SRA Download QC — Failed Sample Identifier"
    log_info "Base path: $BASE_PATH"

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
