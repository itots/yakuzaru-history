#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_FILE="../refseq/GCF_003339765.1_Mmul_10_rm.out"
OUTPUT_DIR="../mask/repeatmask"

mkdir -p "${OUTPUT_DIR}"

mawk 'BEGIN{OFS="\t"}{if(NR>3){print $5, $6-1, $7}}' "${INPUT_FILE}" | \
    sort -k1,1 -k2,2n -S 10G --parallel=10 | \
    bedtools merge -i - | \
    bgzip -c > "${OUTPUT_DIR}/repeatmask.bed.gz"

tabix -p bed "${OUTPUT_DIR}/repeatmask.bed.gz"
