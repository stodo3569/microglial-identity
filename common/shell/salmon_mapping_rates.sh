#!/bin/bash

################################################################################
# Salmon Mapping Rate Summary
################################################################################
#
# Purpose:
#   Extracts the mapping rate from every salmon quantification output and
#   prints a tab-separated table: Study  GSM  mapping_rate(%)
#
# Usage:
#   bash salmon_mapping_rates.sh [--base-path /data]
#
# Output:
#   Tab-separated to stdout (pipe to file or column as needed):
#     bash salmon_mapping_rates.sh | column -t
#     bash salmon_mapping_rates.sh > /data/salmon_mapping_rates.tsv
#
# Expected directory structure:
#   /data/<Study>/Aligned_data/<GSM>/aux_info/meta_info.json
#
################################################################################

BASE_PATH="/data"

while [[ $# -gt 0 ]]; do
    case $1 in
        --base-path) BASE_PATH="$2"; shift 2 ;;
        --help|-h)
            echo "Usage: bash salmon_mapping_rates.sh [--base-path /data]"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Header
printf "Study\tGSM\tMapping_Rate(%%)\n"

for meta_file in "$BASE_PATH"/*/Aligned_data/*/aux_info/meta_info.json; do
    [[ ! -f "$meta_file" ]] && continue

    # Extract Study and GSM from path
    # .../Study/Aligned_data/GSM/aux_info/meta_info.json
    gsm=$(basename "$(dirname "$(dirname "$meta_file")")")
    study=$(basename "$(dirname "$(dirname "$(dirname "$(dirname "$meta_file")")")")")

    # Extract mapping rate
    if command -v python3 &>/dev/null; then
        rate=$(python3 -c "
import json, sys
try:
    data = json.load(open('$meta_file'))
    pm = data.get('percent_mapped', None)
    if pm is not None:
        print(f'{pm:.2f}')
    else:
        np = data.get('num_processed', 0)
        nm = data.get('num_mapped', 0)
        print(f'{nm/np*100:.2f}' if np > 0 else 'NA')
except Exception:
    print('NA')
" 2>/dev/null)
    else
        rate=$(grep -oP '"percent_mapped"\s*:\s*\K[0-9.]+' "$meta_file" | head -1)
        rate=${rate:-NA}
    fi

    printf "%s\t%s\t%s\n" "$study" "$gsm" "$rate"
done | sort -t$'\t' -k1,1 -k2,2
