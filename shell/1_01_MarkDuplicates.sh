#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../bam/bwa_mem"
OUTPUT_DIR="../bam/mark_dp"
SAMPLE_LIST="../list/sample_rg.txt"
PARALLEL_JOBS=6
MEMORY="80G"

mkdir -p "${OUTPUT_DIR}"

# Function to process BAM files with MarkDuplicates
process_sample() {
    local sample=$1
    local rg_existence=$2
    local input_bam_list=()

    if [ "${rg_existence}" == "rg_exists" ]; then
        while read -r input_bam; do
            input_bam_list+=("-I" "${input_bam}")
        done < <(find "${INPUT_DIR}/${sample}" -name '*.bam') # Use process substitution for safety
    else
        input_bam_list=("-I" "${INPUT_DIR}/${sample}.bam")
    fi

    echo "Processing ${sample}"
    
    if ! gatk --java-options "-Xmx${MEMORY}" MarkDuplicates \
        "${input_bam_list[@]}" \
        -O "${OUTPUT_DIR}/${sample}.bam" \
        -M "${OUTPUT_DIR}/${sample}.metrics" &> "${OUTPUT_DIR}/${sample}.log"; then
        echo -e "\033[31mError: MarkDuplicates failed for ${sample}" >&2
        return 1
    fi

    samtools index "${OUTPUT_DIR}/${sample}.bam"
}

export -f process_sample
export INPUT_DIR OUTPUT_DIR MEMORY

# Run MarkDuplicates in parallel using GNU Parallel
parallel -j"${PARALLEL_JOBS}" --colsep '\t' -a "${SAMPLE_LIST}" process_sample {1} {2}
