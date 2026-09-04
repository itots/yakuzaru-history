#!/bin/bash

set -e 
set -o pipefail

# Define variables
INPUT_DIR="../stats/data/plink/out"
OUTPUT_DIR="../stats/fst"
MASTER_LIST="../list/master/master_kinship_filtered.txt"

mkdir -p "${OUTPUT_DIR}/info" "${OUTPUT_DIR}/out"

# Create cluster file for FST calculation
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $7=="fuscata" {
    print "0", $1, $10;
}' "${MASTER_LIST}" > "${OUTPUT_DIR}/info/cluster.txt"

# Calculate FST
plink2 \
    --bfile ../stats/data/plink/out/jm_ld_pruned \
    --chr-set 20 --allow-extra-chr \
    --autosome \
    --within ../stats/fst/info/cluster.txt \
    --fst CATPHENO method=hudson \
    --out ../stats/fst/out/fst_hudson \
      &> "${OUTPUT_DIR}/out/fst_hudson.log"

plink2 \
    --bfile ../stats/data/plink/out/jm_ld_pruned \
    --chr-set 20 --allow-extra-chr \
    --autosome \
    --within ../stats/fst/info/cluster.txt \
    --fst CATPHENO method=wc \
    --out ../stats/fst/out/fst_wc \
      &> "${OUTPUT_DIR}/out/fst_wc.log"

