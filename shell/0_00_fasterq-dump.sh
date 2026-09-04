#!/bin/bash

set -eo pipefail

# Define directories and files
OUTPUT_DIR="../fastq/rawdata/ncbi_sra"
SRAID_LIST="../../../list/sra.txt"

cd ${OUTPUT_DIR}

while read sra_id; do
    echo "${sra_id}"

    if ! prefetch "${sra_id}" --max-size u -L info &> ${sra_id}_prefetch.log; then
        echo "Failed to prefetch ${sra_id}"
        exit 1
    fi
    if ! fasterq-dump -e62 "${sra_id}" -L info &> ${sra_id}_fasterq-dump.log; then
        echo "Failed to run fasterq-dump for ${sra_id}"
        exit 1
    fi
    pigz -p62 "${sra_id}"*.fastq

done < ${SRAID_LIST}

cd ../../../shell


