#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/masked/all"
OUTPUT_DIR="../vcf/bial_snp"
CHROM_NAME="../vcf/annotated/chromAlias/refseq2ncbi.txt"
CHROM_LIST="../list/chromosome.txt"
JM_LIST="../list/jm.txt"
PARALLEL_JOBS=21
THREADS=3

# Make sure the output directory exists
mkdir -p "${OUTPUT_DIR}/all"
mkdir -p "${OUTPUT_DIR}/snp"
mkdir -p "${OUTPUT_DIR}/snp_chr"
mkdir -p "${OUTPUT_DIR}/snp_chr_jm"
mkdir -p "${OUTPUT_DIR}/snp_chr_jm_concat"
mkdir -p "${OUTPUT_DIR}/snp_chr_concat"

process_chrom () {
    local chrom=$1

    echo "Processing ${chrom}"

    # Step 1: Extract biallelic SNPs + RefHom with low missingness
    bcftools view --threads "${THREADS}" -Oz -o "${OUTPUT_DIR}/all/${chrom}.vcf.gz" "${INPUT_DIR}/${chrom}.vcf.gz" \
        -i '(TYPE="snp" || TYPE="ref") && F_MISSING<0.2' \
        -m 1 -M 2 \
        &> "${OUTPUT_DIR}/all/${chrom}.log"
    bcftools index --threads "${THREADS}" -t "${OUTPUT_DIR}/all/${chrom}.vcf.gz"

    # Step 2: Extract biallelic SNPs
    bcftools view --threads "${THREADS}" -Oz -o "${OUTPUT_DIR}/snp/${chrom}.vcf.gz" "${OUTPUT_DIR}/all/${chrom}.vcf.gz" \
        -m 2 -M 2 \
        &> "${OUTPUT_DIR}/snp/${chrom}.log"
    bcftools index --threads "${THREADS}" -t "${OUTPUT_DIR}/snp/${chrom}.vcf.gz"

    # Step 3: Rename chromosome (RefSeq -> Chr1, etc.)
    bcftools annotate --threads "${THREADS}" -Oz -o "${OUTPUT_DIR}/snp_chr/${chrom}.vcf.gz" \
        --rename-chrs "${CHROM_NAME}" \
        "${OUTPUT_DIR}/snp/${chrom}.vcf.gz" \
        &> "${OUTPUT_DIR}/snp_chr/${chrom}.log"
    bcftools index --threads "${THREADS}" -t "${OUTPUT_DIR}/snp_chr/${chrom}.vcf.gz"

    # Step 4: Extract Japanese macaques and Recalculate Tags
    bcftools view \
        -S "${JM_LIST}" --force-samples \
        "${OUTPUT_DIR}/snp_chr/${chrom}.vcf.gz" | \
        bcftools +fill-tags - -- -t AN,AC,AF,F_MISSING,MAF,HWE | \
        bcftools view -c1 --threads "${THREADS}" -Oz -o "${OUTPUT_DIR}/snp_chr_jm/${chrom}.vcf.gz" \
        &> "${OUTPUT_DIR}/snp_chr_jm/${chrom}.log"
    bcftools index --threads "${THREADS}" -t "${OUTPUT_DIR}/snp_chr_jm/${chrom}.vcf.gz"
}

export -f process_chrom
export INPUT_DIR OUTPUT_DIR CHROM_NAME JM_LIST THREADS

TARGET_CHROMS=$(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_027914.1")

parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" process_chrom ::: ${TARGET_CHROMS}

echo "Concatenating all chromosomes..."

JM_VCFS=()
ALL_VCFS=()

for chrom in ${TARGET_CHROMS}; do
    JM_VCFS+=("${OUTPUT_DIR}/snp_chr_jm/${chrom}.vcf.gz")
    ALL_VCFS+=("${OUTPUT_DIR}/snp_chr/${chrom}.vcf.gz")
done

bcftools concat --threads 8 -Oz -o "${OUTPUT_DIR}/snp_chr_jm_concat/concatenated.vcf.gz" "${JM_VCFS[@]}" \
    &> "${OUTPUT_DIR}/snp_chr_jm_concat/concatenated.log"
bcftools index --threads 8 -t "${OUTPUT_DIR}/snp_chr_jm_concat/concatenated.vcf.gz"

bcftools concat --threads 8 -Oz -o "${OUTPUT_DIR}/snp_chr_concat/concatenated.vcf.gz" "${ALL_VCFS[@]}" \
    &> "${OUTPUT_DIR}/snp_chr_concat/concatenated.log"
bcftools index --threads 8 -t "${OUTPUT_DIR}/snp_chr_concat/concatenated.vcf.gz"
