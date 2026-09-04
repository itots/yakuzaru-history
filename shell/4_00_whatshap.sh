#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/masked/var"
OUTPUT_DIR="../vcf/whatshap"
BAMDIR="../bam/mark_dp"
CHROM_LIST="../list/chromosome.txt"
MASTER_LIST="../list/master/master_qc_filtered.txt"
JM_LIST="../list/jm.txt"
PEDFILE="../list/trio.ped"
REFSEQ="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"

PARALLEL_JOBS=30

mkdir -p "${OUTPUT_DIR}/merged"
mkdir -p "${OUTPUT_DIR}/tmp"
mkdir -p "${OUTPUT_DIR}/info"

# Prepare sample lists
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $7=="fuscata" {print $1}' "${MASTER_LIST}" > "${JM_LIST}"

awk '{print $2"\n"$3"\n"$4}' "${PEDFILE}" | sort | uniq > "${OUTPUT_DIR}/info/trio_samples.txt"

grep -v -f "${OUTPUT_DIR}/info/trio_samples.txt" "${JM_LIST}" > "${OUTPUT_DIR}/info/isolate_samples.txt"

# Chromosomes to process (excluding Y and MT for phasing)
CHROMS=$(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_027914.1|NC_041774.1")

for chrom in ${CHROMS}; do
    mkdir -p "${OUTPUT_DIR}/tmp/${chrom}"
    mkdir -p "${OUTPUT_DIR}/phased/${chrom}"
done

# Function to phase variants for a sample-chromosome
phasing_whatshap () {
    local sample=$1
    local chrom=$2

    local input_vcf="${INPUT_DIR}/${chrom}.vcf.gz"
    local input_bam="${BAMDIR}/${sample}.bam"
    local tmp_vcf="${OUTPUT_DIR}/tmp/${chrom}/${sample}.vcf.gz"
    local out_vcf="${OUTPUT_DIR}/phased/${chrom}/${sample}.vcf.gz"

    if [[ ! -f "${input_vcf}" || ! -f "${input_bam}" ]]; then
        echo "Warning: Missing VCF or BAM for ${sample} on ${chrom}. Skipping."
        return 0
    fi

    echo "Processing ${sample} for chromosome ${chrom}"

    # Extract sample-specific variants
    bcftools view -Oz -o "${tmp_vcf}" \
        --samples "${sample}" \
        --threads 2 \
        "${input_vcf}"
    bcftools index --threads 2 -f -t "${tmp_vcf}"

    # Phasing via Whatshap
    if ! whatshap phase \
        -o "${out_vcf}" \
        --reference="${REFSEQ}" \
        --chromosome="${chrom}" \
        "${tmp_vcf}" \
        "${input_bam}" &> "${OUTPUT_DIR}/phased/${chrom}/${sample}.log"; then
        echo "Error: whatshap failed for ${sample} on ${chrom}" >&2
        exit 1
    fi
    
    bcftools index --threads 2 -f -t "${out_vcf}"
}

export -f phasing_whatshap
export INPUT_DIR OUTPUT_DIR BAMDIR REFSEQ

parallel --colsep '\t' --halt now,fail=1 -j"${PARALLEL_JOBS}" phasing_whatshap {1} {2} \
    :::: "${OUTPUT_DIR}/info/isolate_samples.txt" \
    ::: ${CHROMS}

# Merge phased VCFs for each chromosome
echo "Merging phased VCFs for each chromosome..."

merge_chroms() {
    local chrom=$1
        
    echo "Merging ${chrom}..."
    bcftools merge --threads 3 -Oz -o "${OUTPUT_DIR}/merged/${chrom}.vcf.gz" \
        "${OUTPUT_DIR}/phased/${chrom}"/*.vcf.gz
        
    bcftools index --threads 3 -f -t "${OUTPUT_DIR}/merged/${chrom}.vcf.gz"
}

export -f merge_chroms
export OUTPUT_DIR

PARALLEL_JOBS=20
parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" merge_chroms ::: ${CHROMS}

rm -r "${OUTPUT_DIR}/tmp"

