#!/bin/bash

set -euo pipefail

# Define variables
INPUT_DIR="../bam/mark_dp"
OUTPUT_DIR="../bam/sexing"
PARALLEL_JOBS=43

# Create required directories
mkdir -p "${OUTPUT_DIR}/idxstats"

export OUTPUT_DIR

# Generate idxstats after all indexing is complete
find "${INPUT_DIR}" -name "*.bam" | parallel -j"${PARALLEL_JOBS}" \
    'samtools idxstats "{}" > "${OUTPUT_DIR}/idxstats/{/.}.txt"'

# Run the Python script
python ../python/sam_sexing/sam_sexing.py -i "${OUTPUT_DIR}/idxstats" -y NC_027914.1 -x NC_041774.1 -p -o "${OUTPUT_DIR}"
