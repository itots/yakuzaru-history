#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/masked/var"
OUTPUT_DIR="../vcf/trio_phase"
BAMDIR="../bam/mark_dp"
CHROM_LIST="../list/chromosome.txt"
PEDFILE="../list/trio.ped"
REFSEQ="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"
PARALLEL_JOBS=20
RECOMBRATE=0.433 # Xue et al. 2016

# Create output directories
mkdir -p "${OUTPUT_DIR}/tmp"

phasing_trio () {
    local chrom=$1

    local -a SAMPLES_WH=()
    local -a SAMPLES=()
    local -a BAMS=()

    while read -r _ child father mother _ ; do
        SAMPLES_WH+=("--sample" "${child}" "--sample" "${father}" "--sample" "${mother}")
        SAMPLES+=("${child}" "${father}" "${mother}")
        BAMS+=("${BAMDIR}/${child}.bam" "${BAMDIR}/${father}.bam" "${BAMDIR}/${mother}.bam")
    done < "${PEDFILE}"

    local SAMPLES_BT=$(IFS=,; echo "${SAMPLES[*]}")

    bcftools view \
        --samples "${SAMPLES_BT}" \
        --threads 3 -Oz -o "${OUTPUT_DIR}/tmp/${chrom}.vcf.gz" \
        "${INPUT_DIR}/${chrom}.vcf.gz" \
        &> "${OUTPUT_DIR}/tmp/${chrom}.log"
    
    bcftools index --threads 3 -f -t "${OUTPUT_DIR}/tmp/${chrom}.vcf.gz"

    echo "Phasing chromosome ${chrom}"
    
    if ! whatshap phase \
        --ped "${PEDFILE}" \
        --reference="${REFSEQ}" \
        --chromosome="${chrom}" \
        --recombrate "${RECOMBRATE}" \
        "${SAMPLES_WH[@]}" \
        -o "${OUTPUT_DIR}/${chrom}.vcf.gz" \
        "${OUTPUT_DIR}/tmp/${chrom}.vcf.gz" \
        "${BAMS[@]}" &> "${OUTPUT_DIR}/${chrom}.log"; then
        echo "Error: whatshap failed for ${chrom}" >&2
        exit 1
    fi
    
    bcftools index --threads 3 -f -t "${OUTPUT_DIR}/${chrom}.vcf.gz"
}

export -f phasing_trio
export INPUT_DIR OUTPUT_DIR REFSEQ BAMDIR PEDFILE RECOMBRATE

parallel --halt now,fail=1 -j "${PARALLEL_JOBS}" phasing_trio ::: \
    $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1")

rm -r "${OUTPUT_DIR}/tmp"