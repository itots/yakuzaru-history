#!/bin/bash
set -euo pipefail

# Define variables
INPUT_DIR="../vcf/derived_intergenic/vcf"
OUTPUT_DIR="../stats/data/fsc"
CHROM_LIST="../list/chromosome.txt"
YAKU_LIST="../list/yaku.txt"
EASYSFS="../python/easySFS/easySFS.py"
PARALLEL_JOBS=20

# Create output directories
mkdir -p "${OUTPUT_DIR}/vcf/chroms" "${OUTPUT_DIR}/sfs/easysfs" "${OUTPUT_DIR}/list"

awk '{print $1, "Yakushima"}' "${YAKU_LIST}" > "${OUTPUT_DIR}/list/yakushima.pop"

# Function to process each chromosome
process_chrom() {
    local chrom=$1

    echo "Processing ${chrom}"
    mkdir -p "${OUTPUT_DIR}/vcf/chroms"

    # Extract Yakushima samples at these sites
    bcftools view --threads 3 -S "${YAKU_LIST}" \
        "${INPUT_DIR}/${chrom}.vcf.gz" | \
        bcftools +fill-tags - -- -t AN,AC,AF,F_MISSING,MAF,HWE | \
        bcftools view -i 'F_MISSING=0' -Oz -o "${OUTPUT_DIR}/vcf/chroms/${chrom}.vcf.gz" \
        &> "${OUTPUT_DIR}/vcf/chroms/${chrom}.log"
    
    bcftools index --threads 3 -t -f "${OUTPUT_DIR}/vcf/chroms/${chrom}.vcf.gz"
}

export -f process_chrom
export INPUT_DIR OUTPUT_DIR YAKU_LIST

parallel -j "${PARALLEL_JOBS}" process_chrom ::: \
  $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_027914.1|NC_041774.1")

# Concat
echo "Concatenating VCFs..."
find "${OUTPUT_DIR}/vcf/chroms/" -name "*.vcf.gz" | sort -V > "${OUTPUT_DIR}/vcf_list.txt"

bcftools concat --threads 8 -f "${OUTPUT_DIR}/vcf_list.txt" -Oz -o "${OUTPUT_DIR}/vcf/concatenated.vcf.gz" \
  &> "${OUTPUT_DIR}/vcf/concatenated.log"
bcftools index -f -t "${OUTPUT_DIR}/vcf/concatenated.vcf.gz"

# Calculate Total Length (seq length)
total_length=$(bcftools index -n "${OUTPUT_DIR}/vcf/concatenated.vcf.gz")
echo "Total length (including monomorphic): ${total_length}" > "${OUTPUT_DIR}/vcf/total_length.txt"

# Extract Polymorphic SNPs
bcftools view --threads 8 -Oz -o "${OUTPUT_DIR}/vcf/concatenated_snp_noFixed.vcf.gz" \
    -m2 -M2 -c 1 -e 'INFO/AC=INFO/AN' \
    "${OUTPUT_DIR}/vcf/concatenated.vcf.gz" \
    &> "${OUTPUT_DIR}/vcf/concatenated_snp.log"
bcftools index --threads 8 -f -t "${OUTPUT_DIR}/vcf/concatenated_snp_noFixed.vcf.gz"

# easySFS
python "${EASYSFS}" \
  -i "${OUTPUT_DIR}/vcf/concatenated_snp_noFixed.vcf.gz" \
  -p "${OUTPUT_DIR}/list/yakushima.pop" \
  --proj 14 \
  --total-length "${total_length}" \
  -a \
  -f \
  --unfolded \
  -o "${OUTPUT_DIR}/sfs/easysfs"


