#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../stats/data/iqtree/phylip"
OUTPUT_DIR="../stats/iqtree"

mkdir -p "${OUTPUT_DIR}/pre" "${OUTPUT_DIR}/main"

# Create variant sites for ASC model
iqtree \
    -s "${INPUT_DIR}/thined.min4.fasta" \
    -m MFP+ASC \
    -T 10 \
    --prefix "${OUTPUT_DIR}/pre/thined.min4" \
    --redo || true

# Main run
iqtree \
    -s "${OUTPUT_DIR}/pre/thined.min4.varsites.phy" \
    -m MFP+ASC \
    --mrate E,G \
    -b 100 \
    -T 10 \
    --prefix "${OUTPUT_DIR}/main/thined.min4.variants" \
    --redo
