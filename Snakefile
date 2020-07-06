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
        6
    shell:
        "samtools index -@{threads} {input} {output}"
        
rule filter_bam:
    input:
        bam="data/input_alignments/{sample}.bam",
        bai="data/input_alignments/{sample}.bai"
    params:
        min_mapping_qual=config["min_mapping_quality"],
        chrom=config["chrom"]
    output:
        "data/filtered_alignments/{sample}.filtered.bam"
    threads:
        6
    shell:
        "samtools view -b -@{threads} -F 0x0200 -F 0x0100 -F 0x004 -q {params.min_mapping_qual} {input.bam} {params.chrom} > {output}"

# get reads from bams
rule extract_bam_reads:
    input:
        bam="data/filtered_alignments/{sample}.filtered.bam"
    output:
        fasta="data/aligned_reads/{sample}.fasta"
    threads:
        6
    shell:
        "samtools fasta {input.bam} > {output.fasta}"

# k-mer counting with jellyfish
rule count_kmers:
    input:
        fasta="data/aligned_reads/{sample}.fasta"
    params:
        kmer_size=config["kmer_size"]
    output:
        jf="data/kmer_counts/{sample}.jf"
    threads:
        6
    shell:
        "jellyfish count -t {threads} -m {params.kmer_size} -s 1000M -C -o {output.jf} {input.fasta}"

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
        "jellyfish dump -c -o {output.counts} {input.jf}"

# find document frequencies
# note that pre-allocating a buffer (-S) makes
# a massive difference in run time
rule sort_kmers:
    input:
        counts="data/kmer_counts/{sample}.counts"
    output:
        sorted_kmers="data/doc_freq/{sample}.sorted.counts"
    threads:
        6
    shell:
        "cut -f 1 -d \" \" {input.counts} | sort -S 16G --parallel {threads} > {output.sorted_kmers}"

rule merge_kmer_counts:
    input:
        sorted_kmers=expand("data/doc_freq/{sample}.sorted.counts",
                            sample=config["samples"])
    output:
        kmers_doc_freq="data/doc_freq/all_kmers.doc_freq"
    params:
        batch_size=config["merge_batch_size"]
    threads:
        24
    shell:
        "sort --batch-size {params.batch_size} --parallel {threads} -S 96G -m {input.sorted_kmers} | uniq -c > {output.kmers_doc_freq}"
        
# top-level rules
rule setup_inputs:
    input:
        filtered_bams=expand("data/filtered_alignments/{sample}.filtered.bam",
                             sample=config["samples"])

rule count_all_kmers:
    input:
        kmer_counts=expand("data/kmer_counts/{sample}.counts",
                           sample=config["samples"])

