#!/bin/bash

set -e
set -o pipefail

# Define directories and files
INPUT_DIR="../vcf/gvcf"
OUTPUT_DIR="../vcf/combined_gvcf"
MASTER_LIST="../list/master/master_qc_filtered.txt"
SAMPLE_LIST="../list/sample_sex_qc_filtered.txt"
CHROM_LIST="../list/chromosome.txt"
REFSEQ="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"
PARALLEL_JOBS=24
MEMORY="20G"

# Create a directory for the joint genotyping database
mkdir -p "${OUTPUT_DIR}"

# Make sample list with sex
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $6}' "${MASTER_LIST}" > "${SAMPLE_LIST}"

# Function to run CombineGVCFs
combine_gvcfs() {
    local chrom="$1"
    local sex="${2:-}"  # Default to empty (for autosomes and MT)

    echo "Processing: ${chrom} ${sex}"

    # Create a list file for input gVCFs
    local list_file="${OUTPUT_DIR}/${chrom}${sex:+_${sex}}_inputs.list"
    > "${list_file}" # Initialize empty file

    # Read sample list and collect relevant gVCF files
    while read -r sample_id sample_sex; do
        case "${chrom}" in
            "NC_041774.1") # X chromosome (separate for male and female)
                if [ "${sample_sex}" == "${sex}" ]; then
                    echo "${INPUT_DIR}/${sample_id}/${chrom}.g.vcf.gz" >> "${list_file}"
                fi
                ;;
            "NC_027914.1") # Y chromosome (only for males)
                if [ "${sample_sex}" == "male" ]; then
                    echo "${INPUT_DIR}/${sample_id}/${chrom}.g.vcf.gz" >> "${list_file}"
                fi
                ;;
            *) # Autosomes & Mitochondrial DNA (MT)
                echo "${INPUT_DIR}/${sample_id}/${chrom}.g.vcf.gz" >> "${list_file}"
                ;;
        esac
    done < <(awk '{print $1, $2}' "${SAMPLE_LIST}")

    # Check if there are any files to process
    if [ ! -s "${list_file}" ]; then
        echo "No samples found for ${chrom} ${sex}. Skipping."
        return 0
    fi

    # Run CombineGVCFs
    if ! gatk --java-options "-Xmx${MEMORY}" CombineGVCFs \
        -R "${REFSEQ}" \
        -L "${chrom}" \
        -V "${list_file}" \
        -O "${OUTPUT_DIR}/${chrom}${sex:+_${sex}}.g.vcf.gz" \
        &> "${OUTPUT_DIR}/${chrom}${sex:+_${sex}}.log"
    then
        echo "Error: CombineGVCFs failed for ${chrom} ${sex}" >&2
        exit 1
    fi
}

export -f combine_gvcfs
export SAMPLE_LIST REFSEQ OUTPUT_DIR INPUT_DIR MEMORY

{
    # Autosomes & MT (no sex specification)
    awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_041774.1|NC_027914.1"
    
    # X chromosome (male and female)
    echo "NC_041774.1 male"
    echo "NC_041774.1 female"
    
    # Y chromosome (male only)
    echo "NC_027914.1 male"
} | parallel --halt now,fail=1 --colsep ' ' -j"${PARALLEL_JOBS}" combine_gvcfs {1} {2}
