#!/usr/bin/env bats

setup() {
    N_SAMPLES=20
    N_KMERS=5000
    KMER_LENGTH=21
    FRAC_ASSOCIATED=0.1

    export TEST_TEMP_DIR=`mktemp -u --tmpdir kintsugi-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}

    export DATA_PATH="${TEST_TEMP_DIR}/simulation"
    export TRAIN_WORKDIR="${TEST_TEMP_DIR}/train_workdir"
    export TEST_WORKDIR="${TEST_TEMP_DIR}/test_workdir"

    simulate_kmer_data \
	--n-samples ${N_SAMPLES} \
	--n-unique-kmers ${N_KMERS} \
	--kmer-length ${KMER_LENGTH} \
	--fraction-associated-kmers ${FRAC_ASSOCIATED} \
	--simulation-dir ${DATA_PATH}
}

@test "Workflow" {
      run kintsugi_cli \
      	  train \
	  --num-dimensions 8192 \
	  --sig-threshold 0.05 \
	  --n-partitions 4 \
	  --n-cores 1 \
	  --data-description-fl ${DATA_PATH}/data_labels_paths.tsv \
	  --workdir ${TRAIN_WORKDIR} \
	  --model-output-fl ${TEST_TEMP_DIR}/model.pkl

    [ "$status" -eq 0 ]
    [ -e ${TEST_TEMP_DIR}/model.pkl ]

    run kintsugi_cli \
    	inspect-model \
	--model-fl ${TEST_TEMP_DIR}/model.pkl

    [ "$status" -eq 0 ]

    run kintsugi_cli \
    	apply-model \
	--model-fl ${TEST_TEMP_DIR}/model.pkl \
	--workdir ${TEST_WORKDIR} \
	--data-paths-fl ${DATA_PATH}/data_paths.tsv \
	--predictions-output-fl ${TEST_TEMP_DIR}/predictions.tsv \
	--n-cores 1

    [ "$status" -eq 0 ]
    [ -e ${TEST_TEMP_DIR}/predictions.tsv ]
}

