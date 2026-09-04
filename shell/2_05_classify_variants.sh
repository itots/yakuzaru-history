#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/joint"
OUTPUT_DIR="../vcf/classified"
REFSEQ="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"
PARALLEL_JOBS=62
MEMORY="8G"

# Function to classify variants
classify_variants() {
    local inputfile=$1
    local select_type=$2

    local chrom=$(basename "${inputfile}" .vcf.gz)
    local output_dir="${OUTPUT_DIR}/${select_type}"

    # --- X chromosome Special Handling ---
    if [[ "${chrom}" == "NC_041774.1_female" || "${chrom}" == "NC_041774.1_male" ]]; then
        output_dir="${OUTPUT_DIR}/${select_type}/NC_041774.1"
    fi

    # --- Y chromosome Special Handling ---
    if [[ "${chrom}" == "NC_027914.1_male" ]]; then
        local output_chrom="NC_027914.1"
    else
        local output_chrom="${chrom}"
    fi

    # Ensure the output directory exists
    mkdir -p "${output_dir}"

    echo "Selection: Type=${select_type}, Input=${chrom} -> Output=${output_chrom}.vcf.gz"

    # Run GATK SelectVariants
    if ! gatk --java-options "-Xmx${MEMORY}" SelectVariants \
        -R "${REFSEQ}" \
        -V "${inputfile}" \
        -O "${output_dir}/${output_chrom}.vcf.gz" \
        --select-type-to-include "${select_type}" \
        &> "${output_dir}/${output_chrom}.log"; then
        
        echo "Error: Failed to select variants for ${output_chrom}" >&2
        return 1
    fi

}

export -f classify_variants
export INPUT_DIR OUTPUT_DIR REFSEQ MEMORY

# Run parallel processing
echo "Starting parallel classification..."
parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" classify_variants ::: \
    $(find "${INPUT_DIR}" -name "*.vcf.gz") ::: \
    "SNP" "INDEL" "NO_VARIATION"


# --- Post-Processing: Merge X Chromosomes ---
echo "Merging X chromosome components..."

for type in "SNP" "INDEL" "NO_VARIATION"; do
    X_DIR="${OUTPUT_DIR}/${type}/NC_041774.1"
    
    if [[ -f "${X_DIR}/NC_041774.1_female.vcf.gz" && -f "${X_DIR}/NC_041774.1_male.vcf.gz" ]]; then
        echo "Merging ${type} X chromosome..."
        
        MERGED_OUT="${OUTPUT_DIR}/${type}/NC_041774.1.vcf.gz"
        
        bcftools merge --threads 4 -Oz -o "${MERGED_OUT}" \
            "${X_DIR}/NC_041774.1_female.vcf.gz" \
            "${X_DIR}/NC_041774.1_male.vcf.gz"
            
        bcftools index -f -t "${MERGED_OUT}"
        echo "Created: ${MERGED_OUT}"
    else
        echo "Warning: X chromosome files missing for ${type}, skipping merge."
    fi
done

