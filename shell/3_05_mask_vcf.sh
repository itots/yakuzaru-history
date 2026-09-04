#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_VAR_DIR="../vcf/annotated"
INPUT_NONVAR_DIR="../vcf/classified/NO_VARIATION" 

BED_DIR="../mask/merged"
OUTPUT_DIR="../vcf/masked"
CHROM_LIST="../list/chromosome.txt"

PARALLEL_JOBS=30
THREADS=2

# Create output directories for each type
mkdir -p "${OUTPUT_DIR}/split/SNP"
mkdir -p "${OUTPUT_DIR}/split/INDEL"
mkdir -p "${OUTPUT_DIR}/split/NO_VARIATION"
mkdir -p "${OUTPUT_DIR}/split/tmp"
mkdir -p "${OUTPUT_DIR}/var"
mkdir -p "${OUTPUT_DIR}/all"

# Function to mask VCF files
mask_vcf() {
    local type=$1
    local chrom=$2
    local input_vcf=""
    local output_vcf="${OUTPUT_DIR}/split/${type}/${chrom}.vcf.gz"
    local log_file="${OUTPUT_DIR}/split/${type}/${chrom}.log"

    # Set input path based on type
    if [ "${type}" == "NO_VARIATION" ]; then
        input_vcf="${INPUT_NONVAR_DIR}/${chrom}.vcf.gz"
    else
        # SNP / INDEL
        input_vcf="${INPUT_VAR_DIR}/${type}/${chrom}.vcf.gz"
    fi

    # Skip if input doesn't exist (e.g., Y chromosome missing)
    if [ ! -f "${input_vcf}" ]; then
        echo "Warning: Input file not found for ${type} ${chrom}. Skipping."
        return 0
    fi

    # Set bedfile
    if [ "${chrom}" == "NC_027914.1" ]; then
        local bedfile="${BED_DIR}/ychrom.bed.gz"
    elif [ "${chrom}" == "NC_041774.1" ]; then
        local bedfile="${BED_DIR}/xchrom.bed.gz"
    else
        local bedfile="${BED_DIR}/autosome.bed.gz"
    fi

    echo "Masking ${type} - ${chrom}..."
    
    # Mask the VCF file using bcftools
    bcftools filter --threads "${THREADS}" -R "${bedfile}" -Oz -o "${output_vcf}" \
        "${input_vcf}" &> "${log_file}"
        
    bcftools index --threads "${THREADS}" -t "${output_vcf}"
}

export -f mask_vcf
export INPUT_VAR_DIR INPUT_NONVAR_DIR BED_DIR OUTPUT_DIR THREADS

# Target chromosomes (exclude MT)
TARGET_CHROMS=$(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -v "NC_005943.1")

parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" mask_vcf {1} {2} ::: \
    "SNP" "INDEL" "NO_VARIATION" ::: \
    ${TARGET_CHROMS}

# Concatenate masked VCFs
concat_vars() {
    local chrom=$1
    local snp_vcf="${OUTPUT_DIR}/split/SNP/${chrom}.vcf.gz"
    local indel_vcf="${OUTPUT_DIR}/split/INDEL/${chrom}.vcf.gz"
    local novar_vcf="${OUTPUT_DIR}/split/NO_VARIATION/${chrom}.vcf.gz"
    
    local target_indel="${indel_vcf}"
    local target_novar="${novar_vcf}"

    # To prevent header mismatch for Y and X chromosome
    if [[ "${chrom}" == "NC_041774.1" || "${chrom}" == "NC_027914.1" ]]; then
        echo "Aligning sample headers for ${chrom}..."
        local sample_list="${OUTPUT_DIR}/split/tmp/samples_${chrom}.txt"

        target_indel="${OUTPUT_DIR}/split/tmp/indel_${chrom}.vcf.gz"
        target_novar="${OUTPUT_DIR}/split/tmp/novar_${chrom}.vcf.gz"
        
        bcftools query -l "${snp_vcf}" > "${sample_list}"
        
        bcftools view -Oz -o "${target_indel}" -S "${sample_list}" --force-samples "${indel_vcf}"
        bcftools index --threads 3 -t -f "${target_indel}"
        
        bcftools view -Oz -o "${target_novar}" -S "${sample_list}" --force-samples "${novar_vcf}"
        bcftools index --threads 3 -t -f "${target_novar}"
    else
        echo "Concatenating ${chrom}..."
    fi

    # --- Variants (SNP + INDEL) ---
    bcftools concat --threads 3 -Ou --allow-overlaps \
        "${snp_vcf}" \
        "${target_indel}" | \
        bcftools sort --max-mem 25G -Oz -o "${OUTPUT_DIR}/var/${chrom}.vcf.gz" \
        2> "${OUTPUT_DIR}/var/${chrom}.log"
        
    bcftools index --threads 3 -f -t "${OUTPUT_DIR}/var/${chrom}.vcf.gz"

    # --- All (SNP + INDEL + NO_VARIATION) ---
    bcftools concat --threads 3 -Ou --allow-overlaps \
        "${snp_vcf}" \
        "${target_indel}" \
        "${target_novar}" | \
        bcftools sort --max-mem 25G -Oz -o "${OUTPUT_DIR}/all/${chrom}.vcf.gz" \
        2> "${OUTPUT_DIR}/all/${chrom}.log"
        
    bcftools index --threads 3 -f -t "${OUTPUT_DIR}/all/${chrom}.vcf.gz"

}

export -f concat_vars
export OUTPUT_DIR

PARALLEL_JOBS=22

parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" concat_vars ::: ${TARGET_CHROMS}

rm -rf "${OUTPUT_DIR}/split/tmp"
