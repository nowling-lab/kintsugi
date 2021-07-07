# kintsugi
Pipeline for k-mer counting

## Requirements

* [jellyfish](https://github.com/gmarcais/Jellyfish) v2
* [pigz](https://zlib.net/pigz/)
* [samtools](http://www.htslib.org/)

## Running the Pipeline

1. Put BAM files under `data/input_alignments` (directory will need to be created)
1. Update sample list in `config.yaml` file
1. Run `snakemake --cores 16 count_kmers`

## Pipeline Output
The pipeline will output gzipped kmer count text files under `data/merged_counts`.  The first line starts with "kmer" followed by the sample names, separated by tabs.  Each line after that starts with a k-mer string folllowed by the k-mer counts for each sample.

## Running Smoke Tests
To run the smoke tests:

```bash
$ bats smoke-tests/*.bats
```
