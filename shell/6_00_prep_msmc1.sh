#!/bin/bash

set -e 
set -o pipefail

# Define variables
INPUT_VAR_DIR="../vcf/shapeit/bial_snp_phased"
INPUT_ALL_DIR="../vcf/bial_snp/all"
OUTPUT_DIR="../stats/data/msmc"
CHROM_LIST="../list/chromosome.txt"
JM_LIST="../list/jm.txt"
PARALLEL_JOBS=60

# Create output directories
mkdir -p "${OUTPUT_DIR}/mask_ind"
mkdir -p "${OUTPUT_DIR}/vcf_msmc"

process_sample_chrom () {
    local sample=$1
    local chrom=$2

    echo "Processing ${sample} ${chrom}"

    mkdir -p "${OUTPUT_DIR}/mask_ind/${chrom}"
    mkdir -p "${OUTPUT_DIR}/vcf_msmc/${chrom}"

    # Make mask bed file
    bcftools view -Ou --samples "${sample}" \
        "${INPUT_ALL_DIR}/${chrom}.vcf.gz" | \
    bcftools +fill-tags -Ou -- -t F_MISSING | \
    bcftools view -Ou -i 'F_MISSING=0' | \
    bcftools query -f '%CHROM\t%POS0\t%POS\n' | \
    bedtools merge | \
    bgzip > "${OUTPUT_DIR}/mask_ind/${chrom}/${sample}.bed.gz" \
        2> "${OUTPUT_DIR}/mask_ind/${chrom}/${sample}.log"
    
    tabix -f -p bed "${OUTPUT_DIR}/mask_ind/${chrom}/${sample}.bed.gz"

    # Extract a sample
    bcftools view -Ou \
        --samples "${sample}" \
        "${INPUT_VAR_DIR}/${chrom}.vcf.gz" | \
    bcftools +fill-tags -Ou -- -t F_MISSING | \
    bcftools view -Oz -o "${OUTPUT_DIR}/vcf_msmc/${chrom}/${sample}.vcf.gz" \
        -i 'F_MISSING=0' \
        2> "${OUTPUT_DIR}/vcf_msmc/${chrom}/${sample}.log"
        
    bcftools index -t -f "${OUTPUT_DIR}/vcf_msmc/${chrom}/${sample}.vcf.gz"
}

export -f process_sample_chrom
export INPUT_VAR_DIR INPUT_ALL_DIR OUTPUT_DIR CHROM_LIST JM_LIST

TARGET_CHROMS=$(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1")

parallel -j"${PARALLEL_JOBS}" process_sample_chrom {1} {2} \
    :::: <(awk '{print $1}' "${JM_LIST}") \
    ::: ${TARGET_CHROMS}
