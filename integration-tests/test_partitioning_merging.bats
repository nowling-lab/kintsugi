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

@test "K-mer Count Partitioning and Merging" {
    for i in `seq 1 ${N_SAMPLES}`;
    do
	run bash -c "zcat ${KMER_PATH}/sample_${i}_kmer_counts.tsv.gz | partition_kmer_counts --n-partitions 4 --sample-name ${i} --output-dir ${TEST_TEMP_DIR}/kmer_count_partitions"

	[ "$status" -eq 0 ]
	[ -e ${TEST_TEMP_DIR}/kmer_count_partitions/partition_0/sample_${i}.tsv.gz ]
	[ -e ${TEST_TEMP_DIR}/kmer_count_partitions/partition_1/sample_${i}.tsv.gz ]
	[ -e ${TEST_TEMP_DIR}/kmer_count_partitions/partition_2/sample_${i}.tsv.gz ]
	[ -e ${TEST_TEMP_DIR}/kmer_count_partitions/partition_3/sample_${i}.tsv.gz ]
    done

    for i in `seq 1 4`;
    do
	run bash -c "merge_kmer_count_partitions --partition-dir ${TEST_TEMP_DIR}/kmer_count_partitions/partition_${i} > ${TEST_TEMP_DIR}/kmer_count_partitions/partition_${i}_counts.tsv"

	[ "$status" -eq 0 ]
	[ -e ${TEST_TEMP_DIR}/kmer_count_partitions/partition_${i}_counts.tsv ]
    done
}

