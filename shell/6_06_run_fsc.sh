#!/bin/bash
set -euo pipefail

# Define variables
INPUT_DIR="../stats/data/fsc"
OUTPUT_DIR="../stats/fsc"
DATA_NAME="Yakushima"

# Run fastsimcoal2
process_model_run () {
    local model=$1
    local run=$2

    echo "Running fsc for ${model} run ${run}"

    mkdir -p "${OUTPUT_DIR}/model_comp/${model}/run${run}"
   
    cp "${INPUT_DIR}/tpl_est/${model}.tpl" "${OUTPUT_DIR}/model_comp/${model}/run${run}/${DATA_NAME}.tpl"
    cp "${INPUT_DIR}/tpl_est/${model}.est" "${OUTPUT_DIR}/model_comp/${model}/run${run}/${DATA_NAME}.est"
    cp "${INPUT_DIR}/sfs/easysfs/fastsimcoal2/${DATA_NAME}_DAFpop0.obs" "${OUTPUT_DIR}/model_comp/${model}/run${run}/${DATA_NAME}_DAFpop0.obs"
    cd "${OUTPUT_DIR}/model_comp/${model}/run${run}"
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

export -f process_model_run
export INPUT_DIR OUTPUT_DIR DATA_NAME
models=("model0" "model1" "model2" "model3")
parallel -j62 process_model_run ::: "${models[@]}" ::: {1..100}