#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../bam/depth/samples"
OUTPUT_DIR="../mask/depth_stats"
MASTER_LIST="../list/master/master_raw.txt"
SAMPLE_LIST="../list/sample_sex_raw.txt"
PARALLEL_JOBS=4

# Create output directories
mkdir -p "${OUTPUT_DIR}"

# Make sample list with sex
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $6}' "${MASTER_LIST}" > "${SAMPLE_LIST}"

# Function to calculate mode, median, and mean of depth
process_sample () {
    local sample=$1
    local input_file="${INPUT_DIR}/${sample}.per-base.bed.gz"

    rm -f "${OUTPUT_DIR}/${sample}.txt"

    # Chromosome identifiers and output files
    declare -A CHROMOSOMES=(
        ["^NC_0417(5[4-9]|6[0-9]|7[0-3])\\.1"]="autosome"
        ["^NC_041774\\.1"]="xchrom"
        ["^NC_027914\\.1"]="ychrom"
        ["^NC_005943\\.1"]="mt"
    )

    # Process each chromosome
    for CHR_PATTERN in "${!CHROMOSOMES[@]}"; do
        echo "Processing ${CHROMOSOMES[$CHR_PATTERN]} for sample ${sample}" 
        stats=$(zgrep -E "${CHR_PATTERN}" "${input_file}" | \
        mawk '{for(i=$2; i<$3; i++) print $4}' | \
        datamash mode 1 median 1 mean 1)
        echo -e "${CHROMOSOMES[$CHR_PATTERN]}\t${stats}" >> "${OUTPUT_DIR}/${sample}.txt"
    done
}

export -f process_sample
export INPUT_DIR OUTPUT_DIR

parallel -j"${PARALLEL_JOBS}" --colsep '\t' -a "$SAMPLE_LIST" process_sample {1}
