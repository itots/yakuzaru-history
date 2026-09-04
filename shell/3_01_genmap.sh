#!/bin/bash

set -e
set -o pipefail

# Define variables
OUTPUT_DIR="../mask/genmap"
REFSEQ="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"
THREADS=60

mkdir -p "${OUTPUT_DIR}"

#genmap index -F "${REFSEQ}" -I "${OUTPUT_DIR}/index" &> "${OUTPUT_DIR}/index.log"

genmap map -K 100 -E 2 -T "${THREADS}" -I "${OUTPUT_DIR}/index" -O "${OUTPUT_DIR}/map" -t -w -bg -v &> "${OUTPUT_DIR}/map.log"


mawk '{if ($4 == 1) print}' "${OUTPUT_DIR}/map.bedgraph" | sort -k1,1 -k2,2n | bgzip > "${OUTPUT_DIR}/mappable.bed.gz"
tabix -p bed "${OUTPUT_DIR}/mappable.bed.gz"