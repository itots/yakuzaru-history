#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../bam/depth/samples"
OUTPUT_DIR="../mask/accessibility"
UNION_FILE="${OUTPUT_DIR}/union.bed.gz"

# Create output directories
mkdir -p "${OUTPUT_DIR}"

# Concatenate per-base.bed.gz files
bedtools unionbedg -i "${INPUT_DIR}"/*.per-base.bed.gz | bgzip > "${UNION_FILE}"
tabix -p bed "${UNION_FILE}"

# Chromosome identifiers and output files
declare -A CHROMOSOMES=(
    ["^NC_0417(5[4-9]|6[0-9]|7[0-3])\\.1"]="autosome"
    ["^NC_041774\\.1"]="xchrom"
    ["^NC_027914\\.1"]="ychrom"
)

# Process each chromosome
for CHR_PATTERN in "${!CHROMOSOMES[@]}"; do
    echo "Processing chromosome pattern: $CHR_PATTERN"
    sumbedfile="${OUTPUT_DIR}/${CHROMOSOMES[$CHR_PATTERN]}_sum.bed.gz"
    fltbedfile="${OUTPUT_DIR}/${CHROMOSOMES[$CHR_PATTERN]}_flt.bed.gz"
    statsfile="${OUTPUT_DIR}/${CHROMOSOMES[$CHR_PATTERN]}_stats.txt"

    zgrep -E "$CHR_PATTERN" "$UNION_FILE" | \
    mawk 'BEGIN {OFS="\t"} {sum=0; for(i=4; i<=NF; i++) sum+=$i; print $1, $2, $3, sum}' | \
    bgzip > "${sumbedfile}"
    tabix -p bed "${sumbedfile}"

    zcat "${sumbedfile}" | \
    mawk '{for(i=$2; i<$3; i++) print $4}' | \
    datamash mode 1 median 1 mean 1 > "${statsfile}"

    median_dp=$(awk '{print $2}' "${statsfile}")
    min_dp=$((median_dp / 2))
    max_dp=$((median_dp * 2))

    zcat "${sumbedfile}" | \
    mawk -v min_dp="$min_dp" -v max_dp="$max_dp" 'BEGIN {OFS="\t"} {if ($4 >= min_dp && $4 <= max_dp) print $1, $2, $3, $4}' | \
    sort -k1,1 -k2,2n | bgzip > "${fltbedfile}"
    tabix -p bed "${fltbedfile}"
done