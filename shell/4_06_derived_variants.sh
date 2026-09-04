#!/bin/bash
set -euo pipefail

# Define variables
INPUT_SNP_DIR="../vcf/bial_snp/snp"
INPUT_INDEL_DIR="../vcf/bial_indel/indel"
INPUT_ALL_DIR="../vcf/bial_snp/all"
OUTPUT_DIR="../vcf/derived"
MASTER_LIST="../list/master/master_kinship_filtered.txt"
OUTGROUP_LIST="../list/outgroup.txt"
JM_LIST="../list/jm.txt"
CHROM_LIST="../list/chromosome.txt"
PARALLEL_JOBS=21

mkdir -p "${OUTPUT_DIR}"

# Prepare outgroup sample list
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $7!="fuscata" {print $1}' "${MASTER_LIST}" > "${OUTGROUP_LIST}"

# Function to process each chromosome and variant type
extract_derived_variants() {
    local chrom=$1
    local type=$2

    local input_dir output_dir
    
    if [ "${type}" == "SNP" ]; then
        input_dir="${INPUT_SNP_DIR}"
        output_dir="${OUTPUT_DIR}/snp"
    elif [ "${type}" == "INDEL" ]; then
        input_dir="${INPUT_INDEL_DIR}"
        output_dir="${OUTPUT_DIR}/indel"
    elif [ "${type}" == "ALL" ]; then
        input_dir="${INPUT_ALL_DIR}"
        output_dir="${OUTPUT_DIR}/all"
    else
        echo "Unknown type: ${type}" >&2
        exit 1
    fi

    echo "Processing ${chrom} ${type}"
    mkdir -p "${output_dir}"

    # Extract ancestral SNP positions based on outgroup samples
    bcftools view -Ou -S "${OUTGROUP_LIST}" \
        "${input_dir}/${chrom}.vcf.gz" | \
        bcftools +fill-tags - -- -t AN,AC,AF,F_MISSING,MAF,HWE | \
        bcftools view -i 'F_MISSING<0.2 && (ALT="." || AF<0.2)' | \
        bcftools query -f '%CHROM\t%POS\n' > "${output_dir}/${chrom}_ancestral_positions.txt"

    # Extract derived SNPs
    bcftools view -Ou -S "${JM_LIST}" \
        -T "${output_dir}/${chrom}_ancestral_positions.txt" \
        "${input_dir}/${chrom}.vcf.gz" | \
        bcftools +fill-tags - -- -t AN,AC,AF,F_MISSING,MAF,HWE | \
        bcftools view -Oz -o "${output_dir}/${chrom}.vcf.gz" \
        &> "${output_dir}/${chrom}.log"
        
    bcftools index -f -t "${output_dir}/${chrom}.vcf.gz"
}

export -f extract_derived_variants
export INPUT_SNP_DIR INPUT_INDEL_DIR INPUT_ALL_DIR OUTPUT_DIR OUTGROUP_LIST JM_LIST

parallel -j "${PARALLEL_JOBS}" extract_derived_variants ::: \
    $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_027914.1") ::: SNP INDEL ALL
