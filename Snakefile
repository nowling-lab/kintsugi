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

# -S speeds up sort by allow it to pre-allocate space
# sort is blocking -- consumes all info, sorts, and then dumps all info
# -- so it doesn't run at the same time as pigz
rule dump_kmers:
    input:
        jf="data/kmer_counts_{chrom}/{sample}_{chrom}.jf"
    output:
        counts="data/kmer_counts_{chrom}/L{min_count}/{sample}.counts.gz"
    threads:
        6
    shell:
        "jellyfish dump -c -L {wildcards.min_count} {input.jf} | sort -S 16G --parallel {threads} | pigz -p {threads} -c > {output.counts}"

rule merge_counts:
    input:
        lambda w: ["data/kmer_counts_{}/L{}/{}.counts.gz".format(w.chrom, w.min_count, sample) \
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
