configfile: "config.yaml"

rule link_counts:
    input:
        lambda w: config["sample_data_paths"][w.sample_name]
    output:
        "data/sample_kmer_counts/{sample_name}.gz"
    threads:
        1
    shell:
        "ln -s {input} {output}"

rule partition_kmer_counts:
    input:
        "data/sample_kmer_counts/{sample_name}.gz"
    params:
        n_partitions = config["n_partitions"]
    output:
        ["data/partitioned_kmer_counts/partition_%s/{sample_name}.tsv.gz" % part_num for part_num in range(config["n_partitions"])]
    threads:
        1
    shell:
        "zcat {input} | partition_kmer_counts --n-partitions {params.n_partitions} --sample-name {wildcards.sample_name} --output-dir data/partitioned_kmer_counts"

rule merge_partitions:
    input:
        ["data/partitioned_kmer_counts/partition_{part_num}/%s.tsv.gz" % sample_name for sample_name in config["sample_data_paths"].keys()]
    output:
        "data/merged_kmer_counts/partition_{part_num}.gz"
    threads:
        1
    shell:
        "merge_kmer_count_partitions --partition-dir data/partitioned_kmer_counts/partition_{wildcards.part_num} | gzip -c > {output}"

rule extract_features:
    input:
        partition_fl="data/merged_kmer_counts/partition_{part_num}.gz",
        labels_fl=lambda w: config["labels_fl"]
    params:
        n_dimensions=config["n_features"],
        sig_threshold=config["sig_threshold"],
        output_fl=lambda w: "--sig-kmer-output-fl data/significant_kmers/sig_kmers_partition_{}.tsv.gz".format(w.part_num) if config["output_sig_kmers"] else ""
    output:
        "data/extracted_features/partition_{part_num}.pkl"
    threads:
        1
    shell:
        "extract_features --kmer-count-fl {input.partition_fl} --num-dimensions {params.n_dimensions} --sig-threshold {params.sig_threshold} --labels-fl {input.labels_fl} --features-fl {output} {params.output_fl}"

rule merge_features:
    input:
        feature_fls=["data/extracted_features/partition_{}.pkl".format(part_num) for part_num in range(config["n_partitions"])],
        labels_fl=lambda w: config["labels_fl"]
    params:
        n_dimensions=config["n_features"]
    output:
        "data/merged_features.pkl"
    threads:
        1
    shell:
        "merge_features --num-dimensions {params.n_dimensions} --labels-fl {input.labels_fl} --features-fl {input.feature_fls} --merged-features-fl {output}"

rule train_genotype_classifier:
    input:
        features="data/merged_features.pkl",
        labels_fl=lambda w: config["labels_fl"]
    output:
        "data/model.pkl"
    threads:
        1
    shell:
        "train_genotype_classifier --features-fl {input.features} --labels-fl {input.labels_fl} --model-type lr --model-fl {output}"
        
rule run_experiments:
    input:
	    model="data/model.pkl"
