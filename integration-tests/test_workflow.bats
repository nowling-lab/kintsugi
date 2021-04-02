#!/usr/bin/env bats

setup() {
    N_INDIVIDUALS=20
    N_KMERS=5000
    KMER_LENGTH=21
    FRAC_ASSOCIATED=0.1

    export TEST_TEMP_DIR=`mktemp -u --tmpdir kintsugi-tests.XXXX`
    mkdir -p ${TEST_TEMP_DIR}

    export DATA_PATH="${TEST_TEMP_DIR}/simulation"
    export KMER_PATH="${TEST_TEMP_DIR}/simulation/kmer_counts"
    export PHENO_PATH="${TEST_TEMP_DIR}/phenotypes.labels"

    simulate_kmer_data \
	--n-samples ${N_SAMPLES} \
	--n-unique-kmers ${N_KMERS} \
	--kmer-length ${KMER_LENGTH} \
	--fraction-associated-kmers ${FRAC_ASSOCIATED} \
	--simulation-dir ${DATA_PATH}
}

@test "Workflow" {
    run merge_kmer_count_partitions \
	--partition-dir ${KMER_PATH} > \
	gzip -c > ${TEST_TEMP_DIR}/merged_kmer_counts.tsv.gz
    
    [ "$status" -eq 0 ]
    [ -e ${TEST_TEMP_DIR}/merged_kmer_counts.tsv.gz ]
}

