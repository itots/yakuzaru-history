#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../bam/depth/samples"
SAMPLE_LIST="../list/sample_sex_raw.txt"
OUTPUT_DIR="../mask/detect_nonpar"

# Create output directories
mkdir -p "${OUTPUT_DIR}/union" "${OUTPUT_DIR}/autosome_mean" "${OUTPUT_DIR}/sexchrom"

# Function to process samples by sex
process_sex() {
    local sex=$1
    local union_file="${OUTPUT_DIR}/union/${sex}.bed.gz"

    # Create union file for the given sex
    bedtools unionbedg -i $(mawk -v sex=$sex '{if ($2 == sex) print $1}' "${SAMPLE_LIST}" | xargs -I {} echo "${INPUT_DIR}/{}.regions.bed.gz") | gzip > "$union_file"

    # Calculate autosome mean depth
    autosome_mean=$(zgrep -P "^NC_0417(5[4-9]|6[0-9]|7[0-3]).1" "$union_file" | \
    mawk '{sum=0; for(i=4; i<=NF; i++) sum+=$i; print sum}' | \
    datamash mean 1)

    # Process sex chromosomes
    declare -A CHROMOSOMES=(
        ["^NC_041774.1"]="xchrom"
        ["^NC_027914.1"]="ychrom"
    )

    for CHR_PATTERN in "${!CHROMOSOMES[@]}"; do
        zgrep -P "$CHR_PATTERN" "$union_file" | \
        mawk -v x=$autosome_mean '{sum=0; for(i=4; i<=NF; i++) sum+=$i; print $1, $2, $3, sum/x}' | \
        gzip > "${OUTPUT_DIR}/sexchrom/${CHROMOSOMES[$CHR_PATTERN]}_${sex}.bed.gz"
    done
}

export -f process_sex
export INPUT_DIR SAMPLE_LIST OUTPUT_DIR

# Process female and male samples
parallel process_sex ::: "female" "male"
