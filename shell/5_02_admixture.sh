#!/bin/bash

set -e 
set -o pipefail

# Define variables
INPUT_FILE="../stats/data/plink/out/jm_ld_pruned"
OUTPUT_DIR="../stats/admixture"
THREADS=62

mkdir -p "${OUTPUT_DIR}"

# Prepare PLINK files for ADMIXTURE
plink \
    --bfile "${INPUT_FILE}" \
    --autosome \
    --make-bed \
    --out "${OUTPUT_DIR}/jm_ld_pruned_autosome" \
    &> "${OUTPUT_DIR}/jm_ld_pruned_autosome.log"

# Run ADMIXTURE for K=1 to K=5
(
    cd "${OUTPUT_DIR}"
    
    for K in $(seq 1 5); do 
        admixture -j"${THREADS}" --cv "jm_ld_pruned_autosome.bed" "$K" | tee "K${K}.log"
    done

    grep -h "CV error" K*.log > "cross_validation.txt"
)
