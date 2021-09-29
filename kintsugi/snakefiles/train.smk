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
        n_partitions = config["n_partitions"],
        sampling_params = config["sampling_params"]
    output:
        ["data/partitioned_kmer_counts/partition_%s/{sample_name}.tsv.gz" % part_num for part_num in range(config["n_partitions"])]
    threads:
        1
    shell:
        "zcat {input} | partition_kmer_counts --n-partitions {params.n_partitions} --sample-name {wildcards.sample_name} --output-dir data/partitioned_kmer_counts {params.sampling_params}"

rule merge_partitions:
    input:
        ["data/partitioned_kmer_counts/partition_{part_num}/%s.tsv.gz" % sample_name for sample_name in config["sample_data_paths"].keys()]
    output:
        "data/merged_kmer_counts/partition_{part_num}.gz"
    threads:
        1
    shell:
        "merge_kmer_count_partitions --partition-dir data/partitioned_kmer_counts/partition_{wildcards.part_num} | gzip -c > {output}"

rule select_features:
    input:
        "data/merged_kmer_counts/partition_{part_num}.gz"
    params:
        labels_fl=config["labels_fl"],
        threshold=config["sig_threshold"],
        feature_scaling=lambda w: "--feature-scaling {}".format(config["feature_scaling"]) if config["feature_scaling"] else ""
    output:
        "data/significant_kmer_counts/sig_kmers_partition_{part_num}.gz"
    threads:
        1
    shell:
        "kmer_association_testing --kmer-count-fl {input} --sig-kmer-output-fl {output} --labels-fl {params.labels_fl} --sig-threshold {params.threshold} {params.feature_scaling}"

rule train_genotype_classifier:
    input:
        sig_kmers=["data/significant_kmer_counts/sig_kmers_partition_{}.gz".format(i) for i in range(config["n_partitions"])],
        labels_fl=lambda w: config["labels_fl"]
    params:
        feature_scaling=lambda w: "--scaling {}".format(config["feature_scaling"]) if config["feature_scaling"] else "",
        n_features=config["n_features"],
        feature_sample_prob=lambda w: "--feature-sample-prob {}".format(config["feature_sampling_prob"]) if config["feature_sampling_prob"] else ""
    output:
        "data/model.pkl"
    threads:
        1
    shell:
        "train_on_significant --kmer-count-fls {input.sig_kmers} --labels-fl {input.labels_fl} --num-dimensions {params.n_features} --model-output-fl {output} {params.feature_scaling} {params.feature_sample_prob}"
        
rule run_experiments:
    input:
	    model="data/model.pkl"
