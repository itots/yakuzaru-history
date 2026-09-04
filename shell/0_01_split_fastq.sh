#!/bin/bash

# This script splits fastq files by read group
# Note that FASTQ retreaved from NCBI SRA may not have read group information or have extra information in the header
# In the former case, the script will not split the fastq files
# In the latter case, the script will remove the extra information and split the fastq files
# The script is designed for Illumina and DNBSEQ

# References
# https://github.com/bsmn/fastq-split-readgroup-tool
# https://github.com/stevekm/fastq-split

set -e 
set -o pipefail

# Define variables
INPUT_DIR="../fastq/rawdata"
OUTPUT_DIR="../fastq/rg_split"
PARALLEL_JOBS=43

mkdir -p ${OUTPUT_DIR}

LOG_FILE="${OUTPUT_DIR}/split_fastq.log"

# Function to split fastq by read group
split_fastq() {
    local file=$1
    local sample_name=$2
    local f_or_r=$3
    local delimiter=$4
    local flowcell_index=$5
    local lane_index=$6

    zcat ${file} | paste - - - - | \
    mawk -F"\t" -v output_dir="${OUTPUT_DIR}" \
                -v sample_name="${sample_name}" \
                -v f_or_r="${f_or_r}" \
                -v delimiter="${delimiter}" \
                -v flowcell_index="${flowcell_index}" \
                -v lane_index="${lane_index}" '
        {
            read_header=$1
            split(read_header, arr_rg, delimiter)
            read_group=arr_rg[flowcell_index] "." arr_rg[lane_index]

            # To retreave original read header for SRR and DRR
            if (index(read_header, "@SRR") || index(read_header, "@DRR")) {
                split(read_header, arr_hd, " ")
                sub(/^@/, "", arr_hd[1])
                read_header = "@" arr_hd[2] " " arr_hd[1] " " arr_hd[3] # SRA information appended after supposed original readname
            }

            print read_header "\n" $2 "\n" $3 "\n" $4 > output_dir "/" sample_name "/" read_group "_" f_or_r ".fastq"
        }'
}

# Function to process a single sample
process_sample() {
    local fastq_fwd=$1
    local fastq_rev=${fastq_fwd/_1.fastq.gz/_2.fastq.gz}
    local sample_name=$(basename "$fastq_fwd" _1.fastq.gz)
    local header_fwd=$(zcat ${fastq_fwd} | head -n1)

    # Illumina
    if [[ "${header_fwd}" == *":"* ]]; then
        mkdir -p ${OUTPUT_DIR}/${sample_name}
        echo "${sample_name}: split by read group"
        echo -e "${sample_name}\trg_exists\tILLUMINA" >> ${LOG_FILE}
        split_fastq ${fastq_fwd} ${sample_name} 1 ":" 3 4 # forward
        split_fastq ${fastq_rev} ${sample_name} 2 ":" 3 4 # reverse
        pigz "${OUTPUT_DIR}/${sample_name}"/*.fastq
    # DNBSEQ
    elif  [[ "${header_fwd}" == *"V"* ]] && [[ "${header_fwd}" == *"L"* ]] && [[ "${header_fwd}" == *"C"* ]]; then
        mkdir -p ${OUTPUT_DIR}/${sample_name}
        echo "${sample_name}: split by read group"
        echo -e "${sample_name}\trg_exists\tDNBSEQ" >> ${LOG_FILE}
        split_fastq ${fastq_fwd} ${sample_name} 1 "V|L|C" 2 3 # forward
        split_fastq ${fastq_rev} ${sample_name} 2 "V|L|C" 2 3 # reverse
        pigz "${OUTPUT_DIR}/${sample_name}"/*.fastq
    # Missing read group
    else
        echo "${sample_name}: read group not found"
        echo -e "${sample_name}\trg_missing\tNULL" >> ${LOG_FILE}
    fi

    echo "${sample_name}: done"
}

export -f process_sample split_fastq
export OUTPUT_DIR export LOG_FILE

find ${INPUT_DIR} -name "*_1.fastq.gz" | parallel -j"${PARALLEL_JOBS}" process_sample

