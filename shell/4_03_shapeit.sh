#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/merge_wh"
OUTPUT_DIR="../vcf/shapeit"
CHROM_LIST="../list/chromosome.txt"
PARALLEL_JOBS=20
THREADS=3

mkdir -p "${OUTPUT_DIR}/bial_snp"
mkdir -p "${OUTPUT_DIR}/bial_snp_phased"

phasing_shapeit () {
    local chrom=$1

    echo "Extracting biallelic SNPs for ${chrom}"

    bcftools view -m 2 -M 2 -v snps \
        --threads "${THREADS}" \
        -Oz -o "${OUTPUT_DIR}/bial_snp/${chrom}.vcf.gz" \
        "${INPUT_DIR}/${chrom}.vcf.gz" \
        2> "${OUTPUT_DIR}/bial_snp/${chrom}.log"

    bcftools index --threads "${THREADS}" -t -f "${OUTPUT_DIR}/bial_snp/${chrom}.vcf.gz"

    echo "Phasing ${chrom}"

    # Run SHAPEIT
    shapeit4 \
        --input "${OUTPUT_DIR}/bial_snp/${chrom}.vcf.gz" \
        --output "${OUTPUT_DIR}/bial_snp_phased/${chrom}.vcf.gz" \
        --region "${chrom}" \
        --thread "${THREADS}" \
        --effective-size 50000 \
        --use-PS 0.0001 \
        --sequencing \
        --log "${OUTPUT_DIR}/bial_snp_phased/${chrom}.log" 

    bcftools index --threads "${THREADS}" -t -f "${OUTPUT_DIR}/bial_snp_phased/${chrom}.vcf.gz"

}

export -f phasing_shapeit
export INPUT_DIR OUTPUT_DIR THREADS

parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" phasing_shapeit ::: \
    $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1")
