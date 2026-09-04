#!/bin/bash
set -euo pipefail

# Define variables
INPUT_DIR="../stats/data/fsc"
OUTPUT_DIR="../stats/fsc/bs"

DATA_NAME="Yakushima"
FSC="/home/itots/programs/fsc28_linux64/fsc28"

BEST_MODEL="model2" # Best model determined from previous analyses. Adjust if needed.

# Run fastsimcoal2
process_bs_run () {
    local bs=$1
    local run=$2

    echo "Running non-parametric bootstrap fsc for bs ${bs} run ${run}"

    mkdir -p "${OUTPUT_DIR}/bs_${bs}/run${run}"
   
    cp "${INPUT_DIR}/tpl_est/${BEST_MODEL}.tpl" "${OUTPUT_DIR}/bs_${bs}/run${run}/${DATA_NAME}.tpl"
    cp "${INPUT_DIR}/tpl_est/${BEST_MODEL}.est" "${OUTPUT_DIR}/bs_${bs}/run${run}/${DATA_NAME}.est"
    cp "${INPUT_DIR}/bs/sfs/bs_${run}/fastsimcoal2/${DATA_NAME}_DAFpop0.obs" "${OUTPUT_DIR}/bs_${bs}/run${run}/${DATA_NAME}_DAFpop0.obs"
    cd "${OUTPUT_DIR}/bs_${bs}/run${run}"
    fsc28 \
        -t "${DATA_NAME}.tpl" \
        -e "${DATA_NAME}.est" \
        -n 200000 \
        -L 50 \
        -M \
        -d \
        -q \
        &> "fsc.log"
    cd ../../../../../shell

}

export -f process_bs_run
export INPUT_DIR OUTPUT_DIR DATA_NAME FSC BEST_MODEL
parallel -j62 process_bs_run ::: {1..100} ::: {1..100}
