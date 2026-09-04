#!/bin/bash
set -euo pipefail

# Define variables
INPUT_DIR="../stats/data/fsc"
OUTPUT_DIR="../stats/fsc"
DATA_NAME="Yakushima"
BEST_RUNS="81,49,31,44" # Best runs for model0, model1, model2, model3 respectively

# Run fastsimcoal2
process_model_run () {
    local model=$1
    local run=$2

    echo "Running fsc for ${model} run ${run}"

    mkdir -p "${OUTPUT_DIR}/best/${model}/run${run}"

    IFS=',' read -r -a best_runs <<< "${BEST_RUNS}"
    
    if [[ "${model}" == "model0" ]]; then
        best_run=${best_runs[0]}
    elif [[ "${model}" == "model1" ]]; then
        best_run=${best_runs[1]}
    elif [[ "${model}" == "model2" ]]; then
        best_run=${best_runs[2]}
    elif [[ "${model}" == "model3" ]]; then
        best_run=${best_runs[3]}
    fi
    cp "${OUTPUT_DIR}/model_comp/${model}/run${best_run}/${DATA_NAME}/${DATA_NAME}_maxL.par" "${OUTPUT_DIR}/best/${model}/run${run}/${DATA_NAME}_maxL.par"
    cp "${OUTPUT_DIR}/model_comp/${model}/run${best_run}/${DATA_NAME}_DAFpop0.obs" "${OUTPUT_DIR}/best/${model}/run${run}/${DATA_NAME}_maxL_DAFpop0.obs"

    cd "${OUTPUT_DIR}/best/${model}/run${run}"
    fsc28 \
        -i "${DATA_NAME}_maxL.par" \
        -n 10000000 \
        -d \
        -q \
        &> "fsc.log"
    cd ../../../../../shell

}

export -f process_model_run
export INPUT_DIR OUTPUT_DIR DATA_NAME BEST_RUNS
models=("model0" "model1" "model2" "model3")
parallel -j62 process_model_run ::: "${models[@]}" ::: {1..100}