#!/bin/bash
set -euo pipefail

# Define variables
INPUT_DIR="../stats/data/fsc/vcf"
OUTPUT_DIR="../stats/data/fsc/bs"
YAKU_LIST="../stats/data/fsc/list/yakushima.pop"
EASYSFS="../python/easySFS/easySFS.py"

snp_length=$(bcftools view -H "${INPUT_DIR}/concatenated_snp_noFixed.vcf.gz" | wc -l)
total_length=$(bcftools index -n "${INPUT_DIR}/concatenated.vcf.gz")

mkdir -p "${OUTPUT_DIR}/tmp" "${OUTPUT_DIR}/vcf" "${OUTPUT_DIR}/sfs"

bcftools view -h "${INPUT_DIR}/concatenated_snp_noFixed.vcf.gz" > "${OUTPUT_DIR}/tmp/header.txt"
bcftools view -H "${INPUT_DIR}/concatenated_snp_noFixed.vcf.gz" > "${OUTPUT_DIR}/tmp/data.txt"

block_size=$(( snp_length / 100 ))
remainder=$(( snp_length % block_size ))

split -d -a 3 -l "${block_size}" "${OUTPUT_DIR}/tmp/data.txt" "${OUTPUT_DIR}/tmp/block_"
rm -f "${OUTPUT_DIR}/tmp/block_100" # Discard last block if it has remainder

echo "Discard last block (${remainder})"

{
  echo "Total sites: ${snp_length}"
  echo "Block size: ${block_size}"
  echo "Remainder: ${remainder}"
} > "${OUTPUT_DIR}/vcf/block_size.txt"

# Create bootstrap samples

process_rep (){
    local rep=$1

    echo "Creating bootstrap replicate ${rep}"

    cat "${OUTPUT_DIR}/tmp/header.txt" > "${OUTPUT_DIR}/vcf/bs_${rep}.vcf"

    for i in {1..100}; do
        cat `shuf -n1 -e ${OUTPUT_DIR}/tmp/block_*` >> "${OUTPUT_DIR}/vcf/bs_${rep}.vcf"
    done

    bgzip "${OUTPUT_DIR}/vcf/bs_${rep}.vcf"

}

# Calculate SFS for each bootstrap replicate
process_rep_sfs (){
    local rep=$1

    echo "Calculating SFS for bootstrap replicate ${rep}"

    python "${EASYSFS}" \
        -i "${OUTPUT_DIR}/vcf/bs_${rep}.vcf.gz" \
        -p "${YAKU_LIST}" \
        --proj 14 \
        --total-length "${total_length}" \
        -a \
        -f \
        --unfolded \
        -o "${OUTPUT_DIR}/sfs/bs_${rep}"
}

export -f process_rep process_rep_sfs
export INPUT_DIR OUTPUT_DIR YAKU_LIST total_length EASYSFS

#parallel -j62 process_rep ::: {1..100}
parallel -j6 process_rep_sfs ::: {1..100} # To avoid memory issue

