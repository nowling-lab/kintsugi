#!/usr/bin/env python

import argparse

from joblib import dump
from joblib import load
from pybloom import BloomFilter
import scipy.sparse as sp
from sklearn.feature_extraction import FeatureHasher

def load_rand_proj(flname):
    return load(flname)

def load_bloomfilter(flname):
    with open(flname) as fl:
        bfilter = BloomFilter.fromfile(fl)
    return bfilter

def read_kmers(flname):
    processed = 1
    with open(flname) as infl:
        for i, ln in enumerate(infl):
            cols = ln.strip().split()
            yield (processed, cols[0], int(cols[1]))
            processed += 1

def filter_passlist(kmers, passlist):
    for t in kmers:
        processed, kmer, count = t 
        if kmer in passlist:
            yield t

def to_binary(kmers):
    for processed, kmer, count in kmers:
        yield processed, kmer, 1

def progress(kmers):
    next_output = 1
    for processed, kmer, count in kmers:
        if next_output <= processed:
            print "Processed", processed, "kmers"
            next_output *= 2
        yield processed, kmer, count

def get_counts(kmers):
    for processed, kmer, count in kmers:
        yield kmer, count
        
def hash_features(kmers, n_features):
    extractor = FeatureHasher(n_features = n_features,
                              input_type = "pair",
                              non_negative=True)
    features = extractor.transform([kmers])

    return features
        
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rand-proj-fl",
                        type=str)

    parser.add_argument("--feature-matrix",
                        type=str,
                        required=True)
    
    parser.add_argument("--passlist-bf",
                        type=str,
                        help="Optional passlist in the form a bloom filter")

    parser.add_argument("--kmer-freq-fl",
                        type=str,
                        required=True)

    parser.add_argument("--binary",
                        action="store_true",
                        help="Use binary features")

    parser.add_argument("--n-features",
                        type=int,
                        default=20,
                        help="Number of hashed features in base 2")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    print "Opening kmer file"
    kmers = read_kmers(args.kmer_freq_fl)

    if args.passlist_bf:
        print "Reading passlist bloomfilter"
        passlist = load_bloomfilter(args.passlist_bf)
        kmers = filter_passlist(kmers, passlist)

    if args.binary:
        print "Using binary features"
        kmers = to_binary(kmers)

    #kmers = progress(kmers)
    kmers = get_counts(kmers)

    n_features = 2 ** args.n_features
    features = hash_features(kmers, n_features=n_features)

    print "Done with feature extraction"
    
    if args.rand_proj_fl:
        print "Loading random projection"
        srp = load_rand_proj(args.rand_proj_fl)
        features = srp.transform(features)

    print "Saving feature matrix"
    dump(features, args.feature_matrix)
