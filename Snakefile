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
        chrom=config["chrom"],
        kmer_size=config["kmer_size"]
    output:
        jf="data/kmer_counts/{sample}.jf"
    threads:
        6
    shell:
        "samtools view -b -@{threads} -F 0x0200 -F 0x0100 -F 0x004 -q {params.min_mapping_qual} {input.bam} {params.chrom} | samtools fasta - | jellyfish count -t {threads} -m {params.kmer_size} -s 1000M -C -o {output.jf} /dev/fd/0 {output}"

# k-mer counting with jellyfish
#rule extract_bam_reads:
#    input:
#        bam="data/filtered_alignments/{sample}.filtered.bam"
#    params:
#        kmer_size=config["kmer_size"]
#    output:
#        jf="data/kmer_counts/{sample}.jf"
#    threads:
#        6
#    shell:
#        "samtools fasta {input.bam} | jellyfish count -t {threads} -m {params.kmer_size} -s 1000M -C -o {output.jf} /dev/fd/0"

# jellyfish dump is single threaded
# but specify number of threads to avoid
# overwhelming I/O
rule dump_kmers:
    input:
        jf="data/kmer_counts/{sample}.jf"
    output:
        counts="data/kmer_counts_L{min_count}/{sample}.counts.gz"
    threads:
        6
    shell:
        "jellyfish dump -c -L {wildcards.min_count} {input.jf} | pigz -p {threads} -c > {output.counts}"

rule merge_counts:
    input:
        lambda w: ["data/kmer_counts_L{}/{}.counts.gz".format(w.min_count, sample) \
                   for sample in config["samples"]]
    output:
        "data/merged_counts/kmer_counts_L{min_count}.gz"
    threads:
        6
    shell:
        "scripts/merge_kmer_counts --count-fls {input} | pigz -p {threads} -c > {output}"
        
rule count_kmers:
    input:
        merged_kmer_counts=expand("data/merged_counts/kmer_counts_L{min_count}.gz",
                           min_count=[1, 2, 3, 4, 5])
