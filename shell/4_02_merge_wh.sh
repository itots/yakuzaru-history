#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_WH_DIR="../vcf/whatshap/merged"
INPUT_TRIO_DIR="../vcf/trio_phase"
OUTPUT_DIR="../vcf/merge_wh"
CHROM_LIST="../list/chromosome.txt"

PARALLEL_JOBS=20

mkdir -p "${OUTPUT_DIR}"

# Function to merge phased VCF files for a chromosome
merging_vcfs () {
    local chrom=$1

    echo "Merging phased VCF files for chromosome ${chrom}"

    bcftools merge --threads 3 -Oz -o "${OUTPUT_DIR}/${chrom}.vcf.gz" \
        "${INPUT_WH_DIR}/${chrom}.vcf.gz" \
        "${INPUT_TRIO_DIR}/${chrom}.vcf.gz" \
        &> "${OUTPUT_DIR}/${chrom}.log"

    bcftools index --threads 3 -f -t "${OUTPUT_DIR}/${chrom}.vcf.gz"
}

export -f merging_vcfs
export INPUT_WH_DIR INPUT_TRIO_DIR OUTPUT_DIR

parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" merging_vcfs ::: \
    $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1")

