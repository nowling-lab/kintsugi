#!/usr/bin/env python

import argparse

from pybloom import BloomFilter

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--kmer-doc-freqs",
                        type=str,
                        required=True)

    parser.add_argument("--output-bf",
                        type=str,
                        required=True)

    parser.add_argument("--min-df",
                        type=int,
                        default=2)

    parser.add_argument("--max-elements",
                        type=int,
                        default=2000000000)

    parser.add_argument("--max-error",
                        type=float,
                        default=0.001)
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    progress = 1
    passlist = BloomFilter(capacity = args.max_elements, error_rate = args.max_error)
    with open(args.kmer_doc_freqs) as fl:
        for i, ln in enumerate(fl):
            cols = ln.strip().split()
            if int(cols[0]) >= args.min_df:
                passlist.add(cols[1])

    with open(args.output_bf, "w") as fl:
        passlist.tofile(fl)
