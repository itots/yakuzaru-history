#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../fastq/fastp"
OUTPUT_DIR="../bam/bwa_mem"
REF_GENOME="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"
SAMPLE_LIST="../list/sample_rg.txt"
THREADS=15
PARALLEL_JOBS=4

# Create output directories if they do not exist
mkdir -p "${OUTPUT_DIR}"

# Check if the reference genome index exists; if not, generate it
if [[ ! -f "${REF_GENOME}.bwt.2bit.64" ]]; then
    echo "Indexing reference genome..."
    bwa-mem2 index "${REF_GENOME}"
else
    echo "Reference genome index already exists. Skipping indexing."
fi

do_bwamem() {
    local fastq_fwd=$1
    local fastq_rev=$2
    local rg_info=$3
    local output_prefix=$4

    bwa-mem2 mem -t"${THREADS}" -R "${rg_info}" "${REF_GENOME}" "${fastq_fwd}" "${fastq_rev}" 2> "${output_prefix}.log" |
    samtools sort -@"${THREADS}" -o "${output_prefix}.bam"
}

# Function to process a single sample
process_sample() {
    local sample=$1
    local rg_existence=$2
    local platform=$3
    if [ "${platform}" == NULL ]; then
        platform="ILLUMINA" # Check SRA metadata for platform information
    fi

    echo "Processing ${sample}"

    if [ "${rg_existence}" == "rg_exists" ]; then
        mkdir -p "${OUTPUT_DIR}/${sample}"
        find "${INPUT_DIR}/${sample}" -name '*_1.fastq.gz' | while read -r fastq_fwd; do
            local fastq_fwd="${fastq_fwd}"
            local fastq_rev="${fastq_fwd/_1.fastq.gz/_2.fastq.gz}"
            local file_id=$(basename "${fastq_fwd}" _1.fastq.gz)
            local rg_id=$(basename "${fastq_fwd}" _1.fastq.gz) # File name without _1.fastq.gz is read group ID (flowcell . lane)
            local rg_info="@RG\tID:${rg_id}\tPL:${platform}\tSM:${sample}"
            local output_prefix="${OUTPUT_DIR}/${sample}/${file_id}"
            
            do_bwamem "${fastq_fwd}" "${fastq_rev}" "${rg_info}" "${output_prefix}"
        done
    # If read group information is not available, process the sample as a single library
    else
        local fastq_fwd="${INPUT_DIR}/${sample}_1.fastq.gz"
        local fastq_rev="${INPUT_DIR}/${sample}_2.fastq.gz"
        local rg_info="@RG\tID:${sample}\tPL:${platform}\tSM:${sample}" # ID is set to sample name
        local output_prefix="${OUTPUT_DIR}/${sample}"
        do_bwamem "${fastq_fwd}" "${fastq_rev}" "${rg_info}" "${output_prefix}"
    fi 
}

export -f process_sample do_bwamem
export INPUT_DIR OUTPUT_DIR REF_GENOME THREADS

# Find all paired-end read files and run BWA-MEM2 in parallel
parallel -j"${PARALLEL_JOBS}" --colsep '\t' -a "${SAMPLE_LIST}" process_sample {1} {2} {3}

