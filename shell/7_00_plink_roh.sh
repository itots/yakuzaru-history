set -e 
set -o pipefail

# Define variables
INPUT_DIR="../stats/data/plink/out"
OUTPUT_DIR="../stats/plink_roh"

mkdir -p "${OUTPUT_DIR}"

plink --bfile "${INPUT_DIR}/king_pruned" \
    --chr-set 20 --allow-extra-chr \
    --autosome \
    --homozyg \
    --homozyg-density 1000 \
    --homozyg-kb 100 \
    --out "${OUTPUT_DIR}/roh" \
    &> "${OUTPUT_DIR}/roh.log"
