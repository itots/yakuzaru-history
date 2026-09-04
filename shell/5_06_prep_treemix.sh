#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DATA="../stats/data/plink/out/all_bp_space"
OUTPUT_DIR="../stats/data/treemix"
MASTER_LIST="../list/master/master_kinship_filtered.txt"
SAMPLE_LIST="../list/sample_treemix.txt"

SCRIPT_PLINK2TREEMIX="../python/scripts/plink2treemix.py"
PYTHON2_EXEC="/home/itots/anaconda3/envs/py2/bin/python"

# Create sample list for TreeMix
awk -F'\t' 'BEGIN {OFS="\t"}
NR > 1 {
        pop_name = substr($7, 1, 3) "_" $9
        print $1, pop_name
}' "${MASTER_LIST}" > "${SAMPLE_LIST}"


# Create output directory
mkdir -p "${OUTPUT_DIR}/tmp" "${OUTPUT_DIR}/data"

awk 'BEGIN{OFS="\t"} {print "0", $1, $2}' "${SAMPLE_LIST}" > "${OUTPUT_DIR}/tmp/clust.txt"

plink \
    --bfile "${INPUT_DATA}" \
    --chr-set 20 --allow-extra-chr \
    --autosome \
    --freq \
    --missing \
    --within "${OUTPUT_DIR}/tmp/clust.txt" \
    --out "${OUTPUT_DIR}/tmp/treemix" \
    &> "${OUTPUT_DIR}/tmp/treemix.log"

gzip -f "${OUTPUT_DIR}/tmp/treemix.frq.strat"

# This script should be run using Python 2 environment
$PYTHON2_EXEC "$SCRIPT_PLINK2TREEMIX" \
    "${OUTPUT_DIR}/tmp/treemix.frq.strat.gz" \
    "${OUTPUT_DIR}/data/treemix.frq.gz"

