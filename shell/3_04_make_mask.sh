#!/bin/bash

set -e
set -o pipefail

# Define variables
OUTPUT_DIR="../mask/merged"
GENMAP_BED="../mask/genmap/mappable.bed.gz"
RMSK_BED="../mask/repeatmask/repeatmask.bed.gz"
NONPAR_BED="../mask/detect_nonpar/nonpar.bed.gz"
A_BED="../mask/accessibility/autosome_flt.bed.gz"
X_BED="../mask/accessibility/xchrom_flt.bed.gz"
Y_BED="../mask/accessibility/ychrom_flt.bed.gz"

# Create output directories
mkdir -p "${OUTPUT_DIR}"

bedtools multiinter -i "${GENMAP_BED}" \
    "${A_BED}" | mawk '$4 == 2' | \
    bedtools subtract -a - -b "${RMSK_BED}" | bgzip > "${OUTPUT_DIR}/autosome.bed.gz"
tabix -p bed "${OUTPUT_DIR}/autosome.bed.gz"

bedtools multiinter -i "${GENMAP_BED}" \
    "${Y_BED}" | mawk '$4 == 2' | \
    bedtools subtract -a - -b "${RMSK_BED}" | bgzip > "${OUTPUT_DIR}/ychrom.bed.gz"
tabix -p bed "${OUTPUT_DIR}/ychrom.bed.gz"

bedtools multiinter -i "${GENMAP_BED}" \
    "${X_BED}" \
    "${NONPAR_BED}" | mawk '$4 == 3' | \
    bedtools subtract -a - -b "${RMSK_BED}" | bgzip > "${OUTPUT_DIR}/xchrom.bed.gz"
tabix -p bed "${OUTPUT_DIR}/xchrom.bed.gz"