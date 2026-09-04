#!/bin/bash

set -e
set -o pipefail

# Define variables
OUTPUT_DIR="../bam/depth"
INPUT_DIR="../bam/mark_dp"
PARALLEL_JOBS=2
THREADS=1

# Create an output directory for depth files
mkdir -p "${OUTPUT_DIR}/samples"

# Define a function to calculate the average depth for a single BAM file
process_sample() {
    local bam_file="$1"
    local fname=$(basename "$bam_file" .bam)  # Extract filename without extension
    local out_prefix="${OUTPUT_DIR}/samples/${fname}"  # Define output prefix for mosdepth

    # Run mosdepth
    mosdepth --by 10000 --use-median -t"${THREADS}" "$out_prefix" "$bam_file"
}

export -f process_sample
export THREADS OUTPUT_DIR

# Process BAM files in parallel
find "${INPUT_DIR}" -name "*.bam" | parallel -j"${PARALLEL_JOBS}" process_sample {}


# Once all BAM files are processed, extract TOTAL mean depth
rm -f "${OUTPUT_DIR}/mean_depth.txt"

find "${OUTPUT_DIR}/samples/" -name "*.mosdepth.summary.txt" | while read -r summary_file; do
    fname=$(basename "${summary_file}" .mosdepth.summary.txt)  # Extract sample name
    mean_depth=$(awk '$1=="total" {print $4}' "${summary_file}")  # Extract TOTAL mean depth
    echo -e "${fname}\t${mean_depth}" >> "${OUTPUT_DIR}/mean_depth.txt"
done