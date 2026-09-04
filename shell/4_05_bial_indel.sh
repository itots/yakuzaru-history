#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/masked/all"
OUTPUT_DIR="../vcf/bial_indel"
CHROM_NAME="../vcf/annotated/chromAlias/refseq2ncbi.txt"
CHROM_LIST="../list/chromosome.txt"
JM_LIST="../list/jm.txt"
PARALLEL_JOBS=21
THREADS=3

# Make sure the output directory exists
mkdir -p "${OUTPUT_DIR}/indel"
mkdir -p "${OUTPUT_DIR}/indel_chr"
mkdir -p "${OUTPUT_DIR}/indel_chr_jm"

process_chrom () {
    local chrom=$1

    echo "Processing ${chrom}"

    # Step 1: Extract biallelic INDELs with F_MISSING < 0.2
    bcftools view --threads "${THREADS}" -Oz -o "${OUTPUT_DIR}/indel/${chrom}.vcf.gz" "${INPUT_DIR}/${chrom}.vcf.gz" \
        -v indels \
        -m 2 -M 2 \
        -i 'F_MISSING < 0.2' \
        &> "${OUTPUT_DIR}/indel/${chrom}.log"
    
    bcftools index --threads "${THREADS}" -t "${OUTPUT_DIR}/indel/${chrom}.vcf.gz"

    # Step 2: Change chromosome name
    bcftools annotate --threads "${THREADS}" -Oz -o "${OUTPUT_DIR}/indel_chr/${chrom}.vcf.gz" \
        --rename-chrs "${CHROM_NAME}" \
        "${OUTPUT_DIR}/indel/${chrom}.vcf.gz" \
        &> "${OUTPUT_DIR}/indel_chr/${chrom}.log"
    
    bcftools index --threads "${THREADS}" -t "${OUTPUT_DIR}/indel_chr/${chrom}.vcf.gz"

    # Step 3: Extract Japanese macaques
    bcftools view \
        -S "${JM_LIST}" --force-samples \
        "${OUTPUT_DIR}/indel_chr/${chrom}.vcf.gz" |\
        bcftools +fill-tags - -- -t AN,AC,AF,F_MISSING,MAF,HWE | \
        bcftools view -c1 --threads "${THREADS}" -Oz -o "${OUTPUT_DIR}/indel_chr_jm/${chrom}.vcf.gz" \
        &> "${OUTPUT_DIR}/indel_chr_jm/${chrom}.log"
    
    bcftools index --threads "${THREADS}" -t "${OUTPUT_DIR}/indel_chr_jm/${chrom}.vcf.gz"
}

export -f process_chrom
export INPUT_DIR OUTPUT_DIR CHROM_NAME THREADS JM_LIST

parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" process_chrom ::: \
    $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_027914.1")
