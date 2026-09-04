#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/basic_filtered/selected"
OUTPUT_DIR="../vcf/annotated"
DATABASE="Mmul_10.105"
CHROMTOUCSC="../python/chromToUcsc/chromToUcsc"
PARALLEL_JOBS=25
CHROM_LIST="../list/chromosome.txt"

# Create output directories
mkdir -p "${OUTPUT_DIR}/chromAlias"
for type in "SNP" "INDEL"; do
    mkdir -p "${OUTPUT_DIR}/tmp/${type}" "${OUTPUT_DIR}/${type}"
done

# Generate chromosome alias file
echo "Generating chromosome alias file..."
${CHROMTOUCSC} --get GCF_003339765.1 -o "${OUTPUT_DIR}/chromAlias/"

# Convert alias file to refseq-to-genbank mapping
awk 'NR>1 {print $1, $5}' "${OUTPUT_DIR}/chromAlias/GCF_003339765.1.chromAlias.tsv" \
> "${OUTPUT_DIR}/chromAlias/refseq2ncbi.txt"

awk 'NR>1 {print $5, $1}' "${OUTPUT_DIR}/chromAlias/GCF_003339765.1.chromAlias.tsv" \
> "${OUTPUT_DIR}/chromAlias/ncbi2refseq.txt"


# Function to annotate VCF files with snpEff
annotating_vcf () {
    local chrom=$1
    local type=$2
    local input_vcf="${INPUT_DIR}/${type}/${chrom}.vcf.gz"

    if [ ! -f "${input_vcf}" ]; then
        echo "Warning: Input file not found for ${type} ${chrom}. Skipping."
        return 0
    fi

    # Rename chromosome names to match the Mmul_10.105 database
    echo "Converting chromosome names for ${chrom} ${type}"
    bcftools annotate --rename-chrs "${OUTPUT_DIR}/chromAlias/refseq2ncbi.txt" \
        -Oz -o "${OUTPUT_DIR}/tmp/${type}/${chrom}.vcf.gz" \
        "${input_vcf}"

    # Annotate with snpEff
    echo "Annotating ${chrom} ${type} with snpEff..."
    snpEff -v \
        -stats "${OUTPUT_DIR}/${type}/${chrom}_summary.html" \
        -csvStats "${OUTPUT_DIR}/${type}/${chrom}_stats.csv" \
        "${DATABASE}" "${OUTPUT_DIR}/tmp/${type}/${chrom}.vcf.gz" \
        2> "${OUTPUT_DIR}/${type}/${chrom}.log" | \
    bcftools annotate --rename-chrs "${OUTPUT_DIR}/chromAlias/ncbi2refseq.txt" \
        -Oz -o "${OUTPUT_DIR}/${type}/${chrom}.vcf.gz" -

    bcftools index -f -t "${OUTPUT_DIR}/${type}/${chrom}.vcf.gz"
}

export -f annotating_vcf
export INPUT_DIR OUTPUT_DIR DATABASE CHROMTOUCSC

parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" annotating_vcf ::: $(awk 'NR>1 {print $7}' "${CHROM_LIST}") ::: "SNP" "INDEL"

rm -r "${OUTPUT_DIR}/tmp"

