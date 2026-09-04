#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../stats/data/treemix/data"
OUTPUT_DIR="../stats/treemix"

mkdir -p "${OUTPUT_DIR}"

run_treemix() {
    local i=$1

    echo "Running TreeMix with m=${i}..."

    treemix -i "${INPUT_DIR}/treemix.frq.gz" \
        -m "${i}" \
        -o "${OUTPUT_DIR}/treemix.${i}" \
        -root fas_Indonesia \
        -k 50 \
        -noss \
        &> "${OUTPUT_DIR}/treemix_${i}.log"
}

export -f run_treemix
export INPUT_DIR OUTPUT_DIR

parallel --halt now,fail=1 -j 6 run_treemix ::: $(seq 0 5)
