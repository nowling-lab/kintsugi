#!/usr/bin/env bats

setup() {
    N_SAMPLES=20
    N_KMERS=5000
    KMER_LENGTH=21
    FRAC_ASSOCIATED=0.1

    export TEST_TEMP_DIR=`mktemp -u --tmpdir kintsugi-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}

    export DATA_PATH="${TEST_TEMP_DIR}/simulation"
    export KMER_PATH="${TEST_TEMP_DIR}/simulation/kmer_counts"
    export PHENO_PATH="${TEST_TEMP_DIR}/simulation/phenotypes.labels"

    simulate_kmer_data \
	--n-samples ${N_SAMPLES} \
	--n-unique-kmers ${N_KMERS} \
	--kmer-length ${KMER_LENGTH} \
	--fraction-associated-kmers ${FRAC_ASSOCIATED} \
	--simulation-dir ${DATA_PATH}
}

@test "Workflow" {
    start=0
    end=${N_SAMPLES}
    for ((i=start; i < end; i++))
    do
	run bash -c "zcat ${KMER_PATH}/sample_${i}_kmer_counts.tsv.gz | partition_kmer_counts --n-partitions 1 --sample-name sample_${i} --output-dir ${TEST_TEMP_DIR}/kmer_count_partitions"

	[ "$status" -eq 0 ]
	[ -e ${TEST_TEMP_DIR}/kmer_count_partitions/partition_0/sample_${i}.tsv.gz ]
    done

    run bash -c "merge_kmer_count_partitions --partition-dir ${TEST_TEMP_DIR}/kmer_count_partitions/partition_0 | gzip -c > ${TEST_TEMP_DIR}/kmer_count_partitions/partition_0_counts.tsv.gz"

    [ "$status" -eq 0 ]
    [ -e ${TEST_TEMP_DIR}/kmer_count_partitions/partition_0_counts.tsv.gz ]

    run extract_features \
	--kmer-count-fl ${TEST_TEMP_DIR}/kmer_count_partitions/partition_0_counts.tsv.gz \
	--num-dimensions 8192 \
	--sig-threshold 0.05 \
	--labels-fl ${PHENO_PATH} \
	--features-fl ${TEST_TEMP_DIR}/partition_features.pkl

    [ "$status" -eq 0 ]
    [ -e ${TEST_TEMP_DIR}/partition_features.pkl ]

    run merge_features \
	--num-dimensions 8192 \
	--labels-fl ${PHENO_PATH} \
	--merged-features-fl ${TEST_TEMP_DIR}/merged_features.pkl \
	--features-fl ${TEST_TEMP_DIR}/partition_features.pkl ${TEST_TEMP_DIR}/partition_features.pkl ${TEST_TEMP_DIR}/partition_features.pkl

    [ "$status" -eq 0 ]
    [ -e ${TEST_TEMP_DIR}/merged_features.pkl ]

    run train_genotype_classifier \
	--features-fl ${TEST_TEMP_DIR}/merged_features.pkl \
	--labels-fl ${PHENO_PATH} \
	--model-type lr \
	--model-fl ${TEST_TEMP_DIR}/model.pkl

    [ "$status" -eq 0 ]
    [ -e ${TEST_TEMP_DIR}/model.pkl ]

    run predict_inversion_genotypes \
	--model-fl ${TEST_TEMP_DIR}/model.pkl \
	--predictions-fl ${TEST_TEMP_DIR}/sample_0_predictions.tsv \
	--kmer-counts-fl ${TEST_TEMP_DIR}/simulation/kmer_counts/sample_0_kmer_counts.tsv.gz

    [ "$status" -eq 0 ]
    [ -e ${TEST_TEMP_DIR}/sample_0_predictions.tsv ]
}

