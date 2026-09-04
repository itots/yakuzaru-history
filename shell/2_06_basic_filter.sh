#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/classified"
OUTPUT_DIR="../vcf/basic_filtered"
CHROM_LIST="../list/chromosome.txt"
REFSEQ="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"
PARALLEL_JOBS=62
MEMORY="8G"

# Ensure output directories exist
for type in "SNP" "INDEL" "NO_VARIATION"; do
    mkdir -p "${OUTPUT_DIR}/filtered/${type}" "${OUTPUT_DIR}/selected/${type}"
done

# Function to filter variants
basic_filtering() {
    local chrom=$1
    local type=$2
    local input_vcf="${INPUT_DIR}/${type}/${chrom}.vcf.gz"

    if [ ! -f "${input_vcf}" ]; then
        echo "Warning: Input file not found for ${type} ${chrom}. Skipping."
        return 0
    fi

    echo "Filtering ${type} variants for chromosome: ${chrom}"

    # SNP
    if [ "${type}" == "SNP" ]; then
        if ! gatk --java-options "-Xmx${MEMORY}" VariantFiltration \
            -V "${input_vcf}" \
            -O "${OUTPUT_DIR}/filtered/${type}/${chrom}.vcf.gz" \
            -R "${REFSEQ}" \
            --filter-name "QD2" --filter-expression "QD < 2.0" \
            --filter-name "SOR3" --filter-expression "SOR > 3.0" \
            --filter-name "FS60" --filter-expression "FS > 60.0" \
            --filter-name "MQ40" --filter-expression "MQ < 40.0" \
            --filter-name "MQRS12.5" --filter-expression "MQRankSum < -12.5" \
            --filter-name "RPRS-8" --filter-expression "ReadPosRankSum < -8.0" \
            &> "${OUTPUT_DIR}/filtered/${type}/${chrom}.log"; then
            echo "Error: Failed to filter ${type} for chromosome ${chrom}" >&2
            exit 1
        fi
    # INDEL
    elif [ "${type}" == "INDEL" ]; then
        if ! gatk --java-options "-Xmx${MEMORY}" VariantFiltration \
            -V "${input_vcf}" \
            -O "${OUTPUT_DIR}/filtered/${type}/${chrom}.vcf.gz" \
            -R "${REFSEQ}" \
            --filter-name "QD2" --filter-expression "QD < 2.0" \
            --filter-name "FS200" --filter-expression "FS > 200.0" \
            --filter-name "SOR10" --filter-expression "SOR > 10.0" \
            --filter-name "RPRS-20" --filter-expression "ReadPosRankSum < -20.0" \
            &> "${OUTPUT_DIR}/filtered/${type}/${chrom}.log"; then
            echo "Error: Failed to filter ${type} for chromosome ${chrom}" >&2
            exit 1
        fi   
    # NO_VARIATION (copy)
    elif [ "${type}" == "NO_VARIATION" ]; then
        echo "Skipping VariantFiltration for NO_VARIATION. Copying directly to selected."
        cp "${input_vcf}" "${OUTPUT_DIR}/selected/${type}/${chrom}.vcf.gz"
        cp "${input_vcf}.tbi" "${OUTPUT_DIR}/selected/${type}/${chrom}.vcf.gz.tbi" 2>/dev/null || true
        return 0 
    fi

    # Select filtered variants
    if ! gatk --java-options "-Xmx${MEMORY}" SelectVariants \
        -V "${OUTPUT_DIR}/filtered/${type}/${chrom}.vcf.gz" \
        -O "${OUTPUT_DIR}/selected/${type}/${chrom}.vcf.gz" \
        -R "${REFSEQ}" \
        --exclude-filtered \
        &> "${OUTPUT_DIR}/selected/${type}/${chrom}.log"; then
        echo "Error: Failed to select filtered ${type} for chromosome ${chrom}" >&2
        exit 1
    fi
}

export -f basic_filtering
export INPUT_DIR OUTPUT_DIR REFSEQ MEMORY

# Run variant filtration in parallel
parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" basic_filtering \
    ::: $(awk 'NR>1 {print $7}' "${CHROM_LIST}") ::: "SNP" "INDEL" "NO_VARIATION"
