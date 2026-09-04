#!/bin/bash

set -e
set -o pipefail

# Define variables
RAWFASTQ_DIR="../fastq/rawdata/ncbi_sra"
SPLITFASTQ_DIR="../fastq/rg_split"
OUTPUT_DIR="../fastq/fastp"
SAMPLE_LIST="../list/sample_rg.txt"
THREADS=15
PARALLEL_JOBS=4

mkdir -p "$OUTPUT_DIR"

# Function to run fastp
do_fastp() {
    local file_id=$1
    local input_dir=$2
    local output_dir=$3

    fastp \
        -i "${input_dir}/${file_id}_1.fastq.gz" \
        -I "${input_dir}/${file_id}_2.fastq.gz" \
        -o "${output_dir}/${file_id}_1.fastq.gz" \
        -O "${output_dir}/${file_id}_2.fastq.gz" \
        -j "${output_dir}/${file_id}.json" \
        -h "${output_dir}/${file_id}.html" \
        -w "${THREADS}" \
        &> "${output_dir}/${file_id}.log"
}

# Function to process a single sample
process_sample() {
    local sample=$1
    local rg_existence=$2

    if [ "${rg_existence}" == "rg_exists" ] ; then
        mkdir -p "${OUTPUT_DIR}/${sample}"
        find "${SPLITFASTQ_DIR}/${sample}" -name '*_1.fastq.gz' | while read -r fastq_fwd; do
            local fastq_fwd="${fastq_fwd}"
            local file_id=$(basename "$fastq_fwd" _1.fastq.gz)
            local input_dir="${SPLITFASTQ_DIR}/${sample}"
            local output_dir="${OUTPUT_DIR}/${sample}"
            
            do_fastp "${file_id}" "$input_dir" "$output_dir"
        done
    else
        do_fastp "${sample}" "${RAWFASTQ_DIR}" "${OUTPUT_DIR}"
    fi
}

export -f process_sample do_fastp
export RAWFASTQ_DIR SPLITFASTQ_DIR OUTPUT_DIR THREADS

parallel -j"${PARALLEL_JOBS}" --colsep '\t' -a "$SAMPLE_LIST" process_sample {1} {2}