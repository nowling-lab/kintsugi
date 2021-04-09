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
        
rule predict_genotypes:
    input:
        sample_counts="data/sample_kmer_counts/{sample_name}.gz",
        model_fl=config["model_fl"]
    output:
        "data/predictions/pred_{sample_name}.tsv"
    threads:
        1
    shell:
        "predict_inversion_genotypes --model-fl {input.model_fl} --kmer-counts-fl {input.sample_counts} --predictions-fl {output} --sample-name {wildcards.sample_name}"

rule merge_predictions:
    input:
        lambda w: expand("data/predictions/pred_{sample_name}.tsv",
                         sample_name=config["sample_data_paths"].keys())
    output:
        "data/all_predictions.tsv"
    threads:
        1
    shell:
        "cat {input} > {output}"

rule run_experiments:
    input:
        predictions="data/all_predictions.tsv"
