#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/bial_snp/all"
OUTPUT_DIR="../stats/pixy"
CHROM_LIST="../list/chromosome.txt"
MASTER_LIST="../list/master/master_qc_filtered.txt"
PIIND_LIST="../list/sample_sample.txt"
PARALLEL_JOBS=20
THREADS=3

awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $1}' "${MASTER_LIST}" > "${PIIND_LIST}"

mkdir -p "${OUTPUT_DIR}/pi_ind/"

# Function to calculate per site per individual pi
process_chrom () {
    local chrom=$1

    echo "Running pixy for ${chrom}..."

    # Run pixy
    pixy \
        --stats pi \
        --vcf "${INPUT_DIR}/${chrom}.vcf.gz" \
        --populations "${PIIND_LIST}" \
        --window_size 1000000000 \
        --n_cores "${THREADS}" \
        --output_folder "${OUTPUT_DIR}/pi_ind/" \
        --output_prefix "${chrom}" \
        &> "${OUTPUT_DIR}/pi_ind/${chrom}.log"     

}

export -f process_chrom
export INPUT_DIR OUTPUT_DIR THREADS PIIND_LIST

parallel --halt now,fail=1 -j "${PARALLEL_JOBS}" process_chrom ::: \
    $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1")
