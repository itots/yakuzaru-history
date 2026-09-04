#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../stats/data/plink/out"
OUTPUT_DIR="../stats/data/ibs_mdist"

mkdir -p "${OUTPUT_DIR}"

SAMPLE_COUNT=$(wc -l < "${INPUT_DIR}/ibs.mdist.id" | awk '{print $1}')

echo "${SAMPLE_COUNT}" > "${OUTPUT_DIR}/ibs.mdist"

awk 'NR==1 {print $2}' "${INPUT_DIR}/ibs.mdist.id" >> "${OUTPUT_DIR}/ibs.mdist"

paste <(awk 'NR>1 {print $2}' "${INPUT_DIR}/ibs.mdist.id") "${INPUT_DIR}/ibs.mdist" >> "${OUTPUT_DIR}/ibs.mdist"

