#!/bin/bash

################################################################################
# Salmon Quantification Quality Control — Failed Sample Identifier
################################################################################
#
# Purpose:
#   Verifies that all expected samples from a salmon quantification run
#   completed successfully. Identifies failures at any pipeline stage
#   (missing quantification output, salmon errors), detects missing samples
#   that were never quantified, and flags abnormalities in salmon's
#   meta_info.json reports (abnormally low mapping rate, zero mapped reads,
#   excessive decoy mapping, library type inconsistencies, etc.). Generates
#   a retry input file for re-quantifying only the failed/problematic samples.
#
# Usage:
#   # Cross-reference against an input file (recommended)
#   bash salmon_identify_fails.sh \
#     --input-file /data/<chapter>/config/datasets.txt
#
#   # Scan /data for failures without an input file (legacy mode)
#   bash salmon_identify_fails.sh
#
# Modes:
#   With --input-file:
#     Reads the expected datasets and samples from the input file, then
#     checks whether EVERY expected sample has salmon quantification
#     output on disk, corresponding meta_info.json reports, and no
#     anomalies in the reports. This catches completely missing datasets
#     and samples that were never even attempted.
#
#   Without --input-file (legacy):
#     Scans /data for dataset directories containing GSE/SRP/EGAD in
#     their name, then checks Aligned_data/ for completeness and
#     anomalies. Cannot detect datasets that failed before creating any
#     directories.
#
# Outputs (written to /data/):
#   salmon_failed_samples_report.txt  — human-readable report of all failures
#   salmon_failed_samples_retry.txt   — input file for re-running quantification
#
# Anomaly checks (per-sample from salmon output files):
#   1. quant.sf missing or empty (salmon produced no quantification)
#   2. aux_info/ directory missing
#   3. meta_info.json missing (salmon did not produce mapping stats)
#   4. Zero processed reads (salmon processed no input)
#   5. Zero mapped reads (no reads mapped to transcriptome)
#   6. Mapping rate < 20% (unusually low alignment)
#   7. Decoy mapping rate > 40% (excessive genomic/decoy contamination)
#   8. Quantified transcripts with non-zero counts < 1000 (suspiciously few)
#   9. Library type inconsistency via lib_format_counts.json
#  10. Entries in salmon_quant_failed.txt error logs
#
# Expected directory structure (created by salmon_quant_docker.sh):
#
#   Salmon quantifies ALL runs for a given GSM sample together into a
#   single output directory. This means quantification output is per-GSM,
#   not per-SRR — even when a GSM has multiple runs (e.g. technical
#   replicates or split lanes), salmon merges them at quantification time.
#
#   /data/
#   └── <Output_Directory>/              e.g. Friedman_2019_GSE125050
#       ├── Trimmed_data/                (input — from fastp_trim_docker.sh)
#       │   └── <GSM>/                   e.g. GSM3559136
#       │       ├── fastp_<SRR>_1.fastq.gz   (may have multiple SRR runs)
#       │       └── fastp_<SRR>_2.fastq.gz
#       ├── Aligned_data/                (salmon quantification output)
#       │   └── <GSM>/                   one output dir per SAMPLE (not per run)
#       │       ├── quant.sf                  transcript-level quantification
#       │       ├── cmd_info.json             salmon command used
#       │       ├── lib_format_counts.json    library type fragment counts
#       │       ├── libParams/                library parameter estimates
#       │       ├── aux_info/
#       │       │   └── meta_info.json        mapping stats (parsed for QC)
#       │       └── logs/
#       │           └── quant.log             salmon log output
#       ├── salmon_quant_summary.txt
#       ├── salmon_error_log.txt
#       ├── salmon_quant_failed.txt
#       ├── salmon_successful_samples.txt
#       └── salmon_skipped_samples.txt
#
#   This script checks:
#     1. Does /data/<Output_Directory>/Aligned_data exist?
#     2. Does each expected GSM have a quantified subdirectory?
#     3. Does each GSM output contain quant.sf, aux_info/, meta_info.json?
#     4. Is quant.sf non-empty and has sufficient quantified transcripts?
#     5. Are there anomalies inside meta_info.json (mapping rate, decoy rate)?
#     6. Is library type consistent (lib_format_counts.json)?
#     7. Are there entries in the *_failed.txt error logs?
#
# Dependencies:
#   Standard bash utilities (awk, grep, find) — no Docker required
#   python3 (for JSON parsing) or awk-based fallback
#
################################################################################

set -u
set -o pipefail

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

BASE_PATH="/data"
OUTPUT_REPORT="salmon_failed_samples_report.txt"
OUTPUT_RETRY_FILE="salmon_failed_samples_retry.txt"
INPUT_FILE=""

# Anomaly thresholds (adjustable)
MIN_MAPPING_RATE=20       # minimum % of reads mapping to transcriptome
MAX_DECOY_RATE=40         # maximum % of reads mapping to decoy sequences
MIN_QUANTIFIED_TX=1000    # minimum non-zero transcripts in quant.sf
MAX_LIBTYPE_DISCORD=50    # maximum % reads disagreeing with expected lib type

# ═══════════════════════════════════════════════════════════════════════════════
# Color codes for terminal output
# ═══════════════════════════════════════════════════════════════════════════════

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
Salmon Quantification QC — Failed Sample Identifier

USAGE:
  With input file (recommended):
    bash salmon_identify_fails.sh --input-file datasets.txt

  Without input file (legacy scan mode):
    bash salmon_identify_fails.sh

OPTIONS:
  --input-file FILE      Cross-reference expected samples from this file
                         against what was actually quantified. Same format as
                         geo_download_docker.sh / salmon_quant_docker.sh input.
  --base-path PATH       Override base data path (default: /data)
  --min-mapping-rate N   Minimum mapping rate % (default: 20)
  --max-decoy-rate N     Maximum decoy mapping rate % (default: 40)
  --min-quant-tx N       Minimum non-zero quantified transcripts (default: 1000)
  --max-libtype-discord N Maximum library type discordance % (default: 50)
  --help, -h             Show this help message

INPUT FILE FORMAT (tab-separated):
  Output_Directory<TAB>GSE_accession<TAB>GSM_samples<TAB>[NGC_file]

  Examples:
    Friedman_2019_GSE125050	GSE125050	all
    Matlik_2023_GSE227729	GSE227729	GSM7106344,GSM7106387

WHAT THIS SCRIPT CHECKS:
  ✔ Aligned_data directory exists per study
  ✔ Every expected GSM has a quantified subdirectory
  ✔ Every expected GSM has a salmon output directory with quant.sf
  ✔ No quant.sf files are empty (0 bytes)
  ✔ meta_info.json reports exist for every run
  ✔ meta_info.json reports show no anomalies:
      - Zero processed reads
      - Zero mapped reads
      - Mapping rate < threshold
      - Decoy mapping rate above threshold
  ✔ quant.sf has sufficient non-zero transcripts
  ✔ Library type consistency (lib_format_counts.json)
  ✔ salmon_quant_failed.txt entries are flagged

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
            --min-mapping-rate)
                MIN_MAPPING_RATE="$2"; shift 2 ;;
            --max-decoy-rate)
                MAX_DECOY_RATE="$2"; shift 2 ;;
            --min-quant-tx)
                MIN_QUANTIFIED_TX="$2"; shift 2 ;;
            --max-libtype-discord)
                MAX_LIBTYPE_DISCORD="$2"; shift 2 ;;
            --help|-h)
                show_help; exit 0 ;;
            *)
                log_error "Unknown option: $1"; show_help; exit 1 ;;
        esac
    done
}

################################################################################
# JSON Parsing Utility — meta_info.json
#   Uses python3 if available, otherwise falls back to awk/grep.
#   Extracts key QC metrics from salmon's meta_info.json.
#
#   Output format (pipe-delimited):
#     num_processed|num_mapped|mapping_rate|num_decoy|decoy_rate|num_valid_targets
#
#   All rates are percentages (0–100). If a field is missing the value
#   is set to "NA".
################################################################################

parse_salmon_meta_json() {
    local json_file="$1"

    if [[ ! -f "$json_file" ]] || [[ ! -s "$json_file" ]]; then
        echo "NA|NA|NA|NA|NA|NA"
        return 1
    fi

    # —— Python path (prefer python3) ——
    if command -v python3 &>/dev/null; then
        python3 << PYEOF
import json, sys

try:
    with open("$json_file") as f:
        data = json.load(f)
except Exception:
    print("NA|NA|NA|NA|NA|NA")
    sys.exit(1)

# —— Reads processed (total input fragments) ——
num_processed = data.get("num_processed", "NA")

# —— Reads mapped (mapped to transcriptome + decoys) ——
num_mapped = data.get("num_mapped", "NA")

# —— Mapping rate ——
# salmon provides percent_mapped directly, but we compute as fallback
percent_mapped = data.get("percent_mapped", "NA")
if percent_mapped == "NA" and num_processed not in ("NA", 0, None) and num_mapped not in ("NA", None):
    percent_mapped = round(num_mapped / num_processed * 100, 2)
elif isinstance(percent_mapped, (int, float)):
    percent_mapped = round(percent_mapped, 2)

# —— Decoy fragments ——
num_decoy = data.get("num_decoy_fragments", "NA")
if num_decoy in (None,):
    num_decoy = "NA"

# —— Decoy rate ——
decoy_rate = "NA"
if num_decoy not in ("NA", None) and num_processed not in ("NA", 0, None):
    decoy_rate = round(num_decoy / num_processed * 100, 2)

# —— Number of valid targets (transcripts with non-zero counts) ——
# This is not in meta_info.json; we output NA here — quant.sf is checked separately
num_valid_targets = "NA"

print(f"{num_processed}|{num_mapped}|{percent_mapped}|{num_decoy}|{decoy_rate}|{num_valid_targets}")
PYEOF
    else
        # —— Fallback: awk/grep-based extraction ——
        local num_processed num_mapped percent_mapped num_decoy

        num_processed=$(grep -oP '"num_processed"\s*:\s*\K[0-9]+' "$json_file" | head -1)
        num_mapped=$(grep -oP '"num_mapped"\s*:\s*\K[0-9]+' "$json_file" | head -1)
        percent_mapped=$(grep -oP '"percent_mapped"\s*:\s*\K[0-9.]+' "$json_file" | head -1)
        num_decoy=$(grep -oP '"num_decoy_fragments"\s*:\s*\K[0-9]+' "$json_file" | head -1)

        num_processed=${num_processed:-NA}
        num_mapped=${num_mapped:-NA}
        percent_mapped=${percent_mapped:-NA}
        num_decoy=${num_decoy:-NA}

        # Calculate decoy rate
        local decoy_rate="NA"
        if [[ "$num_processed" != "NA" ]] && [[ "$num_decoy" != "NA" ]] && [[ "$num_processed" -gt 0 ]]; then
            decoy_rate=$(awk "BEGIN {printf \"%.2f\", ($num_decoy/$num_processed)*100}")
        fi

        # Calculate mapping rate if not directly available
        if [[ "$percent_mapped" == "NA" ]] && \
           [[ "$num_processed" != "NA" ]] && [[ "$num_mapped" != "NA" ]] && \
           [[ "$num_processed" -gt 0 ]]; then
            percent_mapped=$(awk "BEGIN {printf \"%.2f\", ($num_mapped/$num_processed)*100}")
        fi

        echo "${num_processed}|${num_mapped}|${percent_mapped}|${num_decoy}|${decoy_rate}|NA"
    fi
}

################################################################################
# Count Non-Zero Transcripts in quant.sf
#   Counts rows where NumReads (column 5) > 0.
#   Returns the count or "NA" on error.
################################################################################

count_quantified_transcripts() {
    local quant_file="$1"

    if [[ ! -f "$quant_file" ]] || [[ ! -s "$quant_file" ]]; then
        echo "NA"
        return 1
    fi

    # quant.sf format: Name  Length  EffectiveLength  TPM  NumReads
    # Skip header line, count rows where NumReads (col 5) > 0
    local count
    count=$(awk 'NR > 1 && $5 > 0 { n++ } END { print n+0 }' "$quant_file" 2>/dev/null)

    if [[ -z "$count" ]]; then
        echo "NA"
        return 1
    fi

    echo "$count"
}

################################################################################
# Check Library Type Consistency (lib_format_counts.json)
#   Reads lib_format_counts.json and checks whether the dominant library
#   type accounts for at least (100 - MAX_LIBTYPE_DISCORD)% of reads.
#
#   Returns "OK", a discordance flag, or "NA" if file missing.
################################################################################

check_libtype_consistency() {
    local libcount_file="$1"

    if [[ ! -f "$libcount_file" ]] || [[ ! -s "$libcount_file" ]]; then
        echo "NA"
        return 0
    fi

    if command -v python3 &>/dev/null; then
        python3 << PYEOF
import json, sys

try:
    with open("$libcount_file") as f:
        data = json.load(f)
except Exception:
    print("NA")
    sys.exit(0)

# The expected_format and compatible_fragment_ratio fields
# lib_format_counts contains: {read_files, expected_format, compatible_fragment_ratio, ...}
# "compatible_fragment_ratio" is the fraction of reads consistent with the expected lib type
compat_ratio = data.get("compatible_fragment_ratio", None)

if compat_ratio is None:
    # Older salmon versions: count concordant vs discordant
    # Keys are library type strings like "IU", "ISF", "ISR", etc.
    counts = {}
    for key, val in data.items():
        if isinstance(val, (int, float)) and key not in ("read_files",):
            counts[key] = val

    if not counts:
        print("NA")
        sys.exit(0)

    total = sum(counts.values())
    if total == 0:
        print("NA")
        sys.exit(0)

    dominant = max(counts.values())
    dominant_pct = round(dominant / total * 100, 2)
    discord_pct = round(100 - dominant_pct, 2)

    if discord_pct > $MAX_LIBTYPE_DISCORD:
        print(f"libtype_discord({discord_pct}%)")
    else:
        print("OK")
else:
    compat_pct = round(compat_ratio * 100, 2)
    discord_pct = round(100 - compat_pct, 2)

    if discord_pct > $MAX_LIBTYPE_DISCORD:
        print(f"libtype_discord({discord_pct}%)")
    else:
        print("OK")
PYEOF
    else
        # Fallback: grep for compatible_fragment_ratio
        local compat_ratio
        compat_ratio=$(grep -oP '"compatible_fragment_ratio"\s*:\s*\K[0-9.]+' "$libcount_file" | head -1)

        if [[ -z "$compat_ratio" ]]; then
            echo "NA"
            return 0
        fi

        local discord_pct
        discord_pct=$(awk "BEGIN {printf \"%.2f\", (1 - $compat_ratio) * 100}")
        local above
        above=$(awk "BEGIN {print ($discord_pct > $MAX_LIBTYPE_DISCORD) ? 1 : 0}")

        if [[ "$above" == "1" ]]; then
            echo "libtype_discord(${discord_pct}%)"
        else
            echo "OK"
        fi
    fi
}

################################################################################
# Check Anomalies in Salmon Quantification
#   Compares extracted metrics against thresholds.
#   Returns a comma-separated list of anomaly flags, or "OK" if clean.
################################################################################

check_salmon_anomalies() {
    local meta_metrics="$1"   # pipe-delimited: num_processed|num_mapped|mapping_rate|num_decoy|decoy_rate|num_valid_targets
    local quant_tx_count="$2" # integer: non-zero transcript count from quant.sf
    local libtype_result="$3" # string: "OK", "libtype_discord(...)", or "NA"

    local num_processed num_mapped mapping_rate num_decoy decoy_rate _unused
    IFS='|' read -r num_processed num_mapped mapping_rate num_decoy decoy_rate _unused <<< "$meta_metrics"

    local anomalies=()

    # 1. meta_info.json couldn't be parsed at all
    if [[ "$num_processed" == "NA" ]]; then
        echo "meta_json_parse_error"
        return 1
    fi

    # 2. Zero reads processed
    if [[ "$num_processed" -eq 0 ]]; then
        anomalies+=("zero_processed_reads")
    fi

    # 3. Zero mapped reads
    if [[ "$num_mapped" != "NA" ]] && [[ "$num_mapped" -eq 0 ]]; then
        anomalies+=("zero_mapped_reads")
    fi

    # 4. Mapping rate below threshold
    if [[ "$mapping_rate" != "NA" ]]; then
        local below
        below=$(awk "BEGIN {print ($mapping_rate < $MIN_MAPPING_RATE) ? 1 : 0}")
        [[ "$below" == "1" ]] && anomalies+=("low_mapping_rate(${mapping_rate}%)")
    fi

    # 5. Decoy rate above threshold
    if [[ "$decoy_rate" != "NA" ]]; then
        local above
        above=$(awk "BEGIN {print ($decoy_rate > $MAX_DECOY_RATE) ? 1 : 0}")
        [[ "$above" == "1" ]] && anomalies+=("high_decoy_rate(${decoy_rate}%)")
    fi

    # 6. Too few quantified transcripts
    if [[ "$quant_tx_count" != "NA" ]] && [[ "$quant_tx_count" -lt "$MIN_QUANTIFIED_TX" ]]; then
        anomalies+=("few_quant_transcripts(${quant_tx_count})")
    fi

    # 7. Library type discordance
    if [[ "$libtype_result" != "OK" ]] && [[ "$libtype_result" != "NA" ]]; then
        anomalies+=("$libtype_result")
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
#   Checks Aligned_data/<GSM>/ for salmon quantification output.
#   Salmon quantifies all runs for a GSM together, so the output is
#   per-sample (per-GSM), not per-SRR run.
#
#   Prints pipe-delimited failure lines to stdout:
#     <GSM>||<failure_reason>
#
#   Returns 0 if the GSM is clean, 1 otherwise.
################################################################################

analyze_gsm_sample() {
    local gsm="$1"
    local study_path="$2"

    local align_dir="$study_path/Aligned_data/$gsm"

    # —— Case 1: Aligned directory doesn't exist at all ——
    if [[ ! -d "$align_dir" ]]; then
        echo "${gsm}||aligned_dir_missing"
        return 1
    fi

    # —— Case 2: Check quant.sf exists ——
    local quant_file="$align_dir/quant.sf"
    if [[ ! -f "$quant_file" ]]; then
        echo "${gsm}||quant_sf_missing"
        return 1
    fi

    # —— Case 3: Check quant.sf is not empty ——
    if [[ ! -s "$quant_file" ]]; then
        echo "${gsm}||quant_sf_empty"
        return 1
    fi

    # —— Case 4: Check aux_info directory exists ——
    local aux_dir="$align_dir/aux_info"
    if [[ ! -d "$aux_dir" ]]; then
        echo "${gsm}||aux_info_missing"
        return 1
    fi

    # —— Case 5: Check meta_info.json exists ——
    local meta_file="$aux_dir/meta_info.json"
    if [[ ! -f "$meta_file" ]]; then
        echo "${gsm}||meta_json_missing"
        return 1
    fi

    # —— Case 6: Parse meta_info.json and check for anomalies ——
    local meta_metrics
    meta_metrics=$(parse_salmon_meta_json "$meta_file")

    # Count non-zero transcripts in quant.sf
    local quant_tx_count
    quant_tx_count=$(count_quantified_transcripts "$quant_file")

    # Check library type consistency
    # lib_format_counts.json lives at the GSM level, not inside aux_info/
    local libtype_result
    libtype_result=$(check_libtype_consistency "$align_dir/lib_format_counts.json")

    local anomaly_result
    anomaly_result=$(check_salmon_anomalies "$meta_metrics" "$quant_tx_count" "$libtype_result")

    if [[ "$anomaly_result" != "OK" ]]; then
        echo "${gsm}||${anomaly_result}"
        return 1
    fi

    return 0
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
           [[ "$dirname" == "salmon_transcriptome" ]] || \
           [[ "$dirname" == "lost+found" ]]; then
            continue
        fi

        # Only process directories that contain an accession pattern
        if [[ ! "$dirname" =~ GSE[0-9]+ ]] && \
           [[ ! "$dirname" =~ SRP[0-9]+ ]] && \
           [[ ! "$dirname" =~ EGAD[0-9]+ ]]; then
            continue
        fi

        # Check it has Aligned_data or salmon log files
        if [[ -d "$dir/Aligned_data" ]] || \
           [[ -f "$dir/salmon_quant_failed.txt" ]] || \
           [[ -f "$dir/salmon_error_log.txt" ]] || \
           [[ -f "$dir/salmon_quant_summary.txt" ]]; then
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
#   Checks all GSM directories under Aligned_data/ for completeness and
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

    # —— Check salmon_quant_failed.txt for recorded failures ——
    if [[ -f "$dataset_dir/salmon_quant_failed.txt" ]] && \
       [[ -s "$dataset_dir/salmon_quant_failed.txt" ]]; then
        while IFS= read -r gsm; do
            [[ -z "$gsm" ]] && continue
            echo "${gsm}||logged_quant_failure" >> "$temp_failures"
        done < "$dataset_dir/salmon_quant_failed.txt"
    fi

    # —— Check every GSM in Aligned_data ——
    if [[ -d "$dataset_dir/Aligned_data" ]]; then
        while IFS= read -r gsm_dir; do
            local gsm
            gsm=$(basename "$gsm_dir")

            local gsm_output
            gsm_output=$(analyze_gsm_sample "$gsm" "$dataset_dir")

            if [[ -n "$gsm_output" ]]; then
                while IFS= read -r line; do
                    [[ -z "$line" ]] && continue
                    echo "$line" >> "$temp_failures"
                done <<< "$gsm_output"
            fi
        done < <(find "$dataset_dir/Aligned_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    fi

    # —— Also check Trimmed_data for GSMs that never got quantified ——
    if [[ -d "$dataset_dir/Trimmed_data" ]]; then
        while IFS= read -r gsm_dir; do
            local gsm
            gsm=$(basename "$gsm_dir")

            # Skip FastQC/MultiQC subdirectories
            [[ "$gsm" == "FastQC" ]] && continue
            [[ "$gsm" == "MultiQC" ]] && continue

            # Check if this GSM has an Aligned_data counterpart
            if [[ ! -d "$dataset_dir/Aligned_data/$gsm" ]]; then
                # Has trimmed data but no quantification output → missed by salmon
                local trim_fq_count
                trim_fq_count=$(find "$gsm_dir" -maxdepth 1 -name "fastp_*.fastq.gz" 2>/dev/null | wc -l)
                if [[ $trim_fq_count -gt 0 ]]; then
                    if ! grep -q "^${gsm}|" "$temp_failures" 2>/dev/null; then
                        echo "${gsm}||not_quantified" >> "$temp_failures"
                    fi
                fi
            fi
        done < <(find "$dataset_dir/Trimmed_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
    fi

    # —— Deduplicate and count ——
    local deduped
    deduped=$(mktemp)
    sort -u "$temp_failures" > "$deduped" 2>/dev/null
    failed_count=$(wc -l < "$deduped" 2>/dev/null || echo 0)

    # Count successful GSMs (in Aligned_data with quant.sf directly present)
    if [[ -d "$dataset_dir/Aligned_data" ]]; then
        while IFS= read -r gsm_dir; do
            local gsm
            gsm=$(basename "$gsm_dir")

            if [[ -f "$gsm_dir/quant.sf" ]] && [[ -s "$gsm_dir/quant.sf" ]]; then
                # Only count as success if not in the failure list
                if ! grep -q "^${gsm}|" "$deduped" 2>/dev/null; then
                    ((success_count++))
                fi
            fi
        done < <(find "$dataset_dir/Aligned_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null)
    fi

    # —— Print summary line ——
    echo "${dataset_name}|${accession}|${failed_count}|${success_count}"

    # —— Print failure detail lines ——
    if [[ -s "$deduped" ]]; then
        while IFS='|' read -r gsm _unused reason; do
            echo "${dataset_name}|${accession}|${gsm}||${reason}"
        done < "$deduped"
    fi

    # —— Console summary (stderr) ——
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
# ╔══════════════════════════════════════════════════════════════════════════════
#   MODE 1: Input File Cross-Reference
# ╚══════════════════════════════════════════════════════════════════════════════
################################################################################

check_expected_samples_from_input() {
    log_info "Cross-referencing expected samples from input file: $INPUT_FILE"
    log_info ""

    if [[ ! -f "$INPUT_FILE" ]]; then
        log_error "Input file not found: $INPUT_FILE"
        exit 1
    fi

    # —— Global counters ——
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

    # —— Step 1: Aggregate input file by Output_Directory ——
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

    # —— Step 2: Initialise report file ——
    cat > "$BASE_PATH/$OUTPUT_REPORT" << EOF
================================================================================
SALMON QUANTIFICATION — FAILED SAMPLES ANALYSIS REPORT
================================================================================
Generated: $(date)
Base Path: $BASE_PATH
Mode: Input file cross-reference ($INPUT_FILE)

Anomaly Thresholds:
  Minimum mapping rate:       ${MIN_MAPPING_RATE}%
  Maximum decoy rate:         ${MAX_DECOY_RATE}%
  Minimum quantified tx:      ${MIN_QUANTIFIED_TX}
  Maximum libtype discordance:${MAX_LIBTYPE_DISCORD}%

================================================================================
PER-DATASET RESULTS
================================================================================
EOF

    # —— Step 3: Check each expected dataset ——
    while IFS=$'\t' read -r output_dir accession gsm_list ngc_file; do
        [[ -z "$output_dir" ]] && continue
        [[ -z "$accession" ]] && continue

        local dataset_dir="$BASE_PATH/$output_dir"
        local dataset_failed=0
        local dataset_anomalies=0
        local dataset_successful=0
        local dataset_total=0

        log_info "Checking: $output_dir ($accession)" >&2

        # —— Check 1: Dataset directory exists? ——
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

        # —— Check 2: Aligned_data directory exists? ——
        if [[ ! -d "$dataset_dir/Aligned_data" ]]; then
            log_error "  MISSING Aligned_data directory in: $dataset_dir" >&2
            ((datasets_with_failures++))

            echo "$output_dir|$accession|$gsm_list||no_aligned_data" >> "$temp_all_failed"

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
  STATUS: Aligned_data directory MISSING (salmon not run?)
--------------------------------------------------
EOF
            continue
        fi

        # —— Check 3: Inspect salmon error log files ——
        local logged_quant_fails=0
        if [[ -f "$dataset_dir/salmon_quant_failed.txt" ]] && \
           [[ -s "$dataset_dir/salmon_quant_failed.txt" ]]; then
            logged_quant_fails=$(wc -l < "$dataset_dir/salmon_quant_failed.txt")
        fi

        # —— Check 4: Per-sample verification ——
        local build_gsm_list=()

        if [[ "$gsm_list" != "all" ]]; then
            # Specific GSMs listed → check each one
            IFS=',' read -ra build_gsm_list <<< "$gsm_list"
        else
            # "all" → check whatever exists in Trimmed_data + Aligned_data
            if [[ -d "$dataset_dir/Trimmed_data" ]]; then
                while IFS= read -r d; do
                    local g
                    g=$(basename "$d")
                    [[ "$g" == "FastQC" ]] && continue
                    [[ "$g" == "MultiQC" ]] && continue
                    build_gsm_list+=("$g")
                done < <(find "$dataset_dir/Trimmed_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
            fi
            # Also pick up anything in Aligned_data not in Trimmed_data
            if [[ -d "$dataset_dir/Aligned_data" ]]; then
                while IFS= read -r d; do
                    local g
                    g=$(basename "$d")
                    # Add if not already present
                    local found=false
                    for existing in "${build_gsm_list[@]+"${build_gsm_list[@]}"}"; do
                        [[ "$existing" == "$g" ]] && { found=true; break; }
                    done
                    [[ "$found" == false ]] && build_gsm_list+=("$g")
                done < <(find "$dataset_dir/Aligned_data" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | sort)
            fi
        fi

        for gsm in "${build_gsm_list[@]+"${build_gsm_list[@]}"}"; do
            gsm=$(echo "$gsm" | tr -d '[:space:]')
            [[ -z "$gsm" ]] && continue
            ((dataset_total++))
            ((global_total_samples++))

            # Run the per-GSM analysis
            local gsm_output
            gsm_output=$(analyze_gsm_sample "$gsm" "$dataset_dir" 2>/dev/null)
            local gsm_exit=$?

            if [[ $gsm_exit -ne 0 ]] || [[ -n "$gsm_output" ]]; then
                # Parse each failure line (format: GSM||reason)
                while IFS='|' read -r g _unused reason; do
                    [[ -z "$g" ]] && continue
                    echo "$output_dir|$accession|$g||$reason" >> "$temp_all_failed"
                    ((dataset_failed++))
                    ((global_failed++))

                    # Distinguish hard failures from anomalies (quality warnings)
                    case "$reason" in
                        low_mapping_rate*|high_decoy_rate*|few_quant_transcripts*|libtype_discord*)
                            ((dataset_anomalies++))
                            ((global_anomalies++))
                            log_warning "  ANOMALY: $gsm — $reason" >&2
                            ;;
                        *)
                            log_warning "  FAILED: $gsm — $reason" >&2
                            ;;
                    esac
                done <<< "$gsm_output"
            else
                ((dataset_successful++))
                ((global_successful++))
            fi
        done

        # —— Also flag logged quant failures not already caught ——
        if [[ -f "$dataset_dir/salmon_quant_failed.txt" ]] && \
           [[ -s "$dataset_dir/salmon_quant_failed.txt" ]]; then
            while IFS= read -r gsm; do
                [[ -z "$gsm" ]] && continue
                if ! grep -q "^${output_dir}|${accession}|${gsm}|" "$temp_all_failed" 2>/dev/null; then
                    echo "$output_dir|$accession|$gsm||logged_quant_failure" >> "$temp_all_failed"
                    ((dataset_failed++))
                    ((global_failed++))
                fi
            done < "$dataset_dir/salmon_quant_failed.txt"
        fi

        # —— Per-dataset report section ——
        {
            echo ""
            echo "--------------------------------------------------"
            echo "Dataset: $output_dir ($accession)"
            echo "--------------------------------------------------"
            echo "  Expected samples:       $dataset_total"
            echo "  Successful:             $dataset_successful"
            echo "  Failed (hard):          $((dataset_failed - dataset_anomalies))"
            echo "  Anomalies (quality):    $dataset_anomalies"
            echo "  Logged quant failures:  $logged_quant_fails"
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

    # —— Step 4: Global summary in report ——
    cat >> "$BASE_PATH/$OUTPUT_REPORT" << EOF

================================================================================
GLOBAL SUMMARY
================================================================================
Total datasets checked:     $total_datasets
Datasets completely missing:$completely_missing_datasets
Datasets with failures:     $datasets_with_failures
Datasets fully successful:  $((total_datasets - datasets_with_failures))

Total samples checked:      $global_total_samples
Successfully quantified:    $global_successful
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

    # —— Datasets where "all" was specified but entire directory is missing ——
    if [[ -s "$temp_missing_all" ]]; then
        cat >> "$BASE_PATH/$OUTPUT_REPORT" << EOF

NOTE: The following datasets specified "all" samples but are completely
missing. Their sample count could not be determined:
EOF
        cat "$temp_missing_all" >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # —— Step 5: Failed samples detail ——
    cat >> "$BASE_PATH/$OUTPUT_REPORT" << 'EOF'

================================================================================
FAILED / ANOMALOUS SAMPLES DETAIL
================================================================================

Format: Dataset | Accession | GSM | (unused) | Failure_Reason

Failure reasons:
  dataset_missing           — entire study directory absent
  no_aligned_data           — Aligned_data/ directory absent (salmon not run?)
  aligned_dir_missing       — GSM subdirectory absent in Aligned_data/
  not_quantified            — trimmed FASTQs exist but no salmon output found
  quant_sf_missing          — quant.sf not found in salmon output
  quant_sf_empty            — quant.sf is 0 bytes
  aux_info_missing          — aux_info/ directory not found
  meta_json_missing         — meta_info.json not found for this sample
  meta_json_parse_error     — meta_info.json could not be parsed
  zero_processed_reads      — 0 reads processed by salmon
  zero_mapped_reads         — 0 reads mapped to transcriptome
  low_mapping_rate(N%)      — mapping rate below threshold
  high_decoy_rate(N%)       — decoy mapping rate above threshold
  few_quant_transcripts(N)  — fewer than threshold non-zero transcripts
  libtype_discord(N%)       — library type discordance above threshold
  logged_quant_failure      — sample appears in salmon_quant_failed.txt

EOF

    if [[ -s "$temp_all_failed" ]]; then
        sort -t'|' -k1,1 -k3,3 "$temp_all_failed" >> "$BASE_PATH/$OUTPUT_REPORT"
    else
        echo "No failed or anomalous samples found! All quantification was successful." \
            >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # —— Step 6: Generate retry input file ——
    local retry_path="$BASE_PATH/$OUTPUT_RETRY_FILE"

    if [[ -s "$temp_all_failed" ]]; then
        cat > "$retry_path" << 'EOF'
# Salmon Quantification Retry Input File
# Generated automatically to retry only failed/anomalous samples
#
# Format: Output_Directory<TAB>GSE_accession<TAB>GSM_samples
#
# Run with:
#   bash salmon_quant_docker.sh --input-file salmon_failed_samples_retry.txt --force
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
# No failed samples found — all quantification was successful!
EOF
        log_success "No failed samples — all quantification was successful!"
    fi

    # —— Cleanup ——
    rm -f "$temp_all_failed" "$temp_missing_all" "$aggregated"

    # —— Final console summary ——
    echo ""
    echo "=================================================="
    log_info "Salmon QC Analysis Complete"
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
        echo "   2. Re-run quantification with:"
        echo "      docker run --rm -u \"\$(id -u):\$(id -g)\" \\"
        echo "        -v /data:/data \\"
        echo "        stodo3569/salmon-tools:0.1 \\"
        echo "        bash /data/<scripts_dir>/salmon_quant_docker.sh \\"
        echo "        --input-file /data/$OUTPUT_RETRY_FILE --force"
        echo ""
    else
        log_success "All samples completed successfully — no retry needed!"
    fi
}

################################################################################
# ╔══════════════════════════════════════════════════════════════════════════════
#   MODE 2: Legacy Scan (no input file)
# ╚══════════════════════════════════════════════════════════════════════════════
################################################################################

generate_reports() {
    log_warning "Legacy mode — scanning for dataset directories"
    log_warning "Tip: Use --input-file for comprehensive checking"
    log_info ""

    # Initialise output files
    > "$BASE_PATH/$OUTPUT_RETRY_FILE"

    cat > "$BASE_PATH/$OUTPUT_REPORT" << EOF
================================================================================
SALMON QUANTIFICATION — FAILED SAMPLES ANALYSIS REPORT
================================================================================
Generated: $(date)
Base Path: $BASE_PATH
Mode: Directory scan (no input file)

Anomaly Thresholds:
  Minimum mapping rate:        ${MIN_MAPPING_RATE}%
  Maximum decoy rate:          ${MAX_DECOY_RATE}%
  Minimum quantified tx:       ${MIN_QUANTIFIED_TX}
  Maximum libtype discordance: ${MAX_LIBTYPE_DISCORD}%

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

    # —— Global summary ——
    cat >> "$BASE_PATH/$OUTPUT_REPORT" << EOF

================================================================================
GLOBAL SUMMARY
================================================================================
Total datasets analyzed:  $total_datasets
Datasets with failures:   $datasets_with_failures
Datasets fully successful:$((total_datasets - datasets_with_failures))

Total samples:            $global_total_samples
Successfully quantified:  $global_successful
Failed / anomalous:       $global_failed
EOF

    if [[ $global_total_samples -gt 0 ]]; then
        echo "Success rate: $(awk "BEGIN {printf \"%.1f\", ($global_successful/$global_total_samples)*100}")%" \
            >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # —— Failed sample detail ——
    echo "" >> "$BASE_PATH/$OUTPUT_REPORT"
    echo "FAILED / ANOMALOUS SAMPLES DETAIL" >> "$BASE_PATH/$OUTPUT_REPORT"
    echo "Format: Dataset | Accession | GSM | (unused) | Failure_Reason" >> "$BASE_PATH/$OUTPUT_REPORT"
    echo "" >> "$BASE_PATH/$OUTPUT_REPORT"

    if [[ -s "$temp_failed" ]]; then
        sort -t'|' -k1,1 -k3,3 "$temp_failed" >> "$BASE_PATH/$OUTPUT_REPORT"
    else
        echo "No failed samples found!" >> "$BASE_PATH/$OUTPUT_REPORT"
    fi

    # —— Retry file ——
    if [[ -s "$temp_failed" ]]; then
        cat > "$BASE_PATH/$OUTPUT_RETRY_FILE" << 'EOF'
# Salmon Quantification Retry Input File
# Run with: bash salmon_quant_docker.sh --input-file salmon_failed_samples_retry.txt --force
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
        echo "# No failed samples — all quantification was successful!" > "$BASE_PATH/$OUTPUT_RETRY_FILE"
        log_success "No failed samples — all quantification was successful!"
    fi

    rm -f "$temp_failed"

    # Console summary
    echo ""
    echo "=================================================="
    log_info "Salmon QC Analysis Complete"
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
        echo "   2. Re-run quantification with:"
        echo "      bash salmon_quant_docker.sh --input-file /data/$OUTPUT_RETRY_FILE --force"
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

    log_info "Salmon Quantification QC — Failed Sample Identifier"
    log_info "Base path: $BASE_PATH"
    log_info "Thresholds: mapping≥${MIN_MAPPING_RATE}%  decoy≤${MAX_DECOY_RATE}%  tx≥${MIN_QUANTIFIED_TX}  libtype_discord≤${MAX_LIBTYPE_DISCORD}%"
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
