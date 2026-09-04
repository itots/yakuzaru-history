#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DATA="../stats/data/plink/out/king_pruned"
OUTPUT_DIR="../stats/data/admixtools"

mkdir -p "${OUTPUT_DIR}"

plink --bfile "${INPUT_DATA}" \
      --chr-set 20 --allow-extra-chr \
      --autosome \
      --make-bed \
      --out "${OUTPUT_DIR}/A" \
      2> "${OUTPUT_DIR}/A.log"

plink --bfile "${INPUT_DATA}" \
      --chr-set 20 --allow-extra-chr \
      --chr 21 \
      --make-bed \
      --out "${OUTPUT_DIR}/X" \
      2> "${OUTPUT_DIR}/X.log"
