configfile: "config.yaml"

import os

##
## Determine paths of k-mer count files and labels
## from base dir
##

label_fl_path = os.path.join(config["base_dir"],
                             config["labels"])

sample_data_paths = dict()
data_fl_path = os.path.join(config["base_dir"],
                            config["input_data"])
with open(data_fl_path) as fl:
    for ln in fl:
        cols = ln.strip().split()
        sample_data_paths[cols[0]] = os.path.join(config["base_dir"],
                                                  cols[1])      


rule link_counts:
    params:
        input_path = lambda w: sample_data_paths[w.sample_name]
    output:
        "data/sample_kmer_counts/{sample_name}.gz"
    threads:
        1
    shell:
        "ln -s ../../{params.input_path} {output}"

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
        ["data/partitioned_kmer_counts/partition_{part_num}/%s.tsv.gz" % sample_name for sample_name in sample_data_paths.keys()]
    output:
        "data/merged_kmer_counts/partition_{part_num}.gz"
    threads:
        1
    shell:
        "merge_kmer_count_partitions --partition-dir data/partitioned_kmer_counts/partition_{wildcards.part_num} | gzip -c > {output}"
        
rule run_experiments:
    input:
        merged_partitions=expand("data/merged_kmer_counts/partition_{part_num}.gz",
                                 part_num=range(config["n_partitions"]))
