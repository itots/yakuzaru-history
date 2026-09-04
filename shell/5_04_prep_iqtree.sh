#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../stats/data/plink/out/all_bp_space"
OUTPUT_DIR="../stats/data/iqtree"
CHROM_LIST="../list/chromosome.txt"
VCF2PHYLIP=$(realpath "../python/vcf2phylip/vcf2phylip.py")

mkdir -p "${OUTPUT_DIR}/vcf" "${OUTPUT_DIR}/phylip"

plink \
    --bfile "${INPUT_DIR}" \
    --chr-set 20 --allow-extra-chr \
    --autosome \
    --mac 1 \
    --bp-space 5000 \
    --recode vcf-iid bgz \
    --out "${OUTPUT_DIR}/vcf/thined" \
    &> "${OUTPUT_DIR}/vcf/thined.log"

bcftools index --threads 8 -f -t "${OUTPUT_DIR}/vcf/thined.vcf.gz"

python "${VCF2PHYLIP}" -i "${OUTPUT_DIR}/vcf/thined.vcf.gz" \
    --output-folder "${OUTPUT_DIR}/phylip" \
    --fasta \
    &> "${OUTPUT_DIR}/phylip/thined.log"

