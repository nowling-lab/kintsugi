configfile: "config.yaml"

rule link_bam:
    input:
        src=lambda wildcards: "%s/{sample}.bam" % config["input_bam_dir"]
    output:
        "data/input_alignments/{sample}.bam"
    threads:
        1
    shell:
        "ln -s {input.src} {output}"

# filter bams
rule index_bam:
    input:
        "data/input_alignments/{sample}.bam"
    output:
        "data/input_alignments/{sample}.bai"
    threads:
        4
    shell:
        "samtools index -@{threads} {input} {output}"
        
rule filter_bam:
    input:
        bam="data/input_alignments/{sample}.bam",
        bai="data/input_alignments/{sample}.bai"
    params:
        min_mapping_qual=config["min_mapping_quality"],
        chrom=config["chrom"],
        kmer_size=config["kmer_size"]
    output:
        jf="data/kmer_counts/{sample}.jf"
    threads:
        4
    shell:
        "samtools view -b -@{threads} -F 0x0200 -F 0x0100 -F 0x004 -q {params.min_mapping_qual} {input.bam} {params.chrom} | samtools fasta - |  jellyfish count -t {threads} -m {params.kmer_size} -s 1000M -C -o {output.jf} /dev/fd/0"


# jellyfish dump is single threaded
# but specify number of threads to avoid
# overwhelming I/O
rule dump_kmers:
    input:
        jf="data/kmer_counts/{sample}.jf"
    output:
        counts="data/kmer_counts/{sample}.counts"
    threads:
        6
    shell:
        "jellyfish dump -c -L 2 {input.jf} > {output.counts}"

rule extract_features:
    input:
        kmer_counts="data/kmer_counts/{sample}.counts",
 
    params:
        n_features=config["n_features"],
        use_binary=config["use_binary_features"]
    output:
        feature_matrix="data/feature_extraction/{sample}.features"
    threads:
        2
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} {params.use_binary} --feature-matrix {output.feature_matrix}"

#log1p feature scaling Combinations w/PCA

rule log1p_binary_feature:
    input: 
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
        feature_matrix="data/log1p_features_before_after/{sample}.features"
    threads:
        4
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before log1p --feature-scaling-after binary"

rule log1p_count_feature:
    input: 
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
       feature_matrix="data/log1p_features_count/{sample}.features"
    threads:
        4
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before log1p --feature-scaling-after counts"

rule log1p_log1p_feature:
    input: 
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
        feature_matrix="data/log1p_features_log1p/{sample}.features"
    threads:
        4
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before log1p --feature-scaling-after log1p"


#Binarizer feature scaling combinations w/PCA

rule binarizer_binary_feature:
    input:
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
        feature_matrix="data/binarizer_features_binary/{sample}.features"
    threads:
        2
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before binary --feature-scaling-after binary"

rule binarizer_count_feature:
    input:
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
        feature_matrix="data/binarizer_features_count/{sample}.features"
    threads:
        2
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before binary --feature-scaling-after counts"
 
rule binarizer_log1p_feature:
    input:
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
        feature_matrix="data/binarizer_features_log1p/{sample}.features"
    threads:
        2
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before binary --feature-scaling-after log1p"

# Counts feature combination w/PCA

rule counts_count_feature:
    input:
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
        feature_matrix="data/count_features_counts/{sample}.features"
    threads:
        2
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before counts --feature-scaling-after counts"

rule counts_log1p_feature:
    input:
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
        feature_matrix="data/count_features_log1p/{sample}.features"
    threads:
        2
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before counts --feature-scaling-after log1p"


rule counts_binary_feature:
    input:
        kmer_counts="data/kmer_counts/{sample}.counts",
    params:
        n_features=config["n_features"]
    output:
        feature_matrix="data/count_features_binary/{sample}.features"
    threads:
        2
    shell:
        "scripts/feature_extractor.py --kmer-freq-fl {input.kmer_counts} --n-features {params.n_features} --feature-matrix {output.feature_matrix} --feature-scaling-before counts --feature-scaling-after binary"


rule pca_log1p_before:
    input:
        feature_matrices=expand("data/log1p_features_before_after/{sample}.features",
                                sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_before.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"

rule pca_log1p_count:
    input:
        feature_matrices=expand("data/log1p_features_count/{sample}.features",
                            sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_log1p_count.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"

rule pca_log1p_log1p:
    input:
       feature_matrices=expand("data/log1p_features_log1p/{sample}.features",
                            sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_log1p_log1p.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"

rule pca_binarizer_binary:
    input:
        feature_matrices=expand("data/binarizer_features_binary/{sample}.features",
                            sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_binary_binarizer.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"

rule pca_binarizer_log1p:
    input:
        feature_matrices=expand("data/binarizer_features_log1p/{sample}.features",
                            sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_binarizer_log1p.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"

rule pca_binarizer_count:
    input:
        feature_matrices=expand("data/binarizer_features_count/{sample}.features",
                            sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_binarizer_count.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"

rule pca_count_count:
    input:
        feature_matrices=expand("data/count_features_counts/{sample}.features",
                            sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_count_count.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"

rule pca_count_log1p:
    input:
        feature_matrices=expand("data/count_features_log1p/{sample}.features",
                            sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_count_log1p.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"

rule pca_count_binary:
    input:
        feature_matrices=expand("data/count_features_binary/{sample}.features",
                            sample=config["samples"])
    params:
        groups_fl=config["groups_fl"]
    output:
        plot="data/pca/pca_count_binarizer.png"
    threads:
        4
    shell:
        "scripts/pca.py --feature-matrices {input.feature_matrices} --groups-fl {params.groups_fl} --plot-fl {output.plot}"




# top-level rules
rule setup_inputs:
    input:
        filtered_bams=expand("data/filtered_alignments/{sample}.filtered.bam",
                             sample=config["samples"])

rule count_all_kmers:
    input:
        kmer_counts=expand("data/kmer_counts/{sample}.counts",
                           sample=config["samples"])

