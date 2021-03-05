configfile: "config.yaml"

# filter bams
rule index_bam:
    input:
        "data/input_alignments/{sample}.bam"
    output:
        "data/input_alignments/{sample}.bai"
    threads:
        6
    shell:
        "samtools index -@{threads} {input} {output}"
        
rule filter_bam:
    input:
        bam="data/input_alignments/{sample}.bam",
        bai="data/input_alignments/{sample}.bai"
    params:
        min_mapping_qual=config["min_mapping_quality"],
        kmer_size=config["kmer_size"]
    output:
        jf="data/kmer_counts_{chrom}/{sample}_{chrom}.jf"
    threads:
        6
    shell:
        "samtools view -b -@2 -F 0x0200 -F 0x0100 -F 0x004 -q {params.min_mapping_qual} {input.bam} {wildcards.chrom} | samtools fasta - | jellyfish count -t 4 -m {params.kmer_size} -s 1000M -C -o {output.jf} /dev/fd/0 {output}"

# jellyfish dump is single threaded
# but specify number of threads to avoid
# overwhelming I/O
rule dump_kmers:
    input:
        jf="data/kmer_counts_{chrom}/{sample}_{chrom}.jf"
    output:
        counts="data/kmer_counts_{chrom}/{sample}_{chrom}_L{min_count}.counts.gz"
    threads:
        6
    shell:
        "jellyfish dump -c -L {wildcards.min_count} {input.jf} | pigz -p {threads} -c > {output.counts}"

rule merge_counts:
    input:
        lambda w: ["data/kmer_counts_{}/{}_{}_L{}.counts.gz".format(w.chrom, sample, w.chrom, w.min_count) \
                   for sample in config["samples"]]
    output:
        "data/merged_counts/merged_counts_{chrom}_L{min_count}.gz"
    threads:
        6
    shell:
        "scripts/merge_kmer_counts --count-fls {input} | pigz -p {threads} -c > {output}"
        
rule count_kmers:
    input:
        merged_kmer_counts=expand("data/merged_counts/merged_counts_{chrom}_L{min_count}.gz",
                                  chrom=["2L", "2R"],
                                  min_count=[1, 2, 3, 4, 5])
