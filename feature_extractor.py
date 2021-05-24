#!/usr/bin/env python

import argparse

from joblib import dump
from joblib import load
from pybloom import BloomFilter
import scipy.sparse as sp
from sklearn.feature_extraction import FeatureHasher
from sklearn.preprocessing import Binarizer
import numpy as np


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

def progress(kmers):
    next_output = 1
    for processed, kmer, count in kmers:
        if next_output <= processed:
            print "Processed", processed, "kmers"
            next_output *= 2
        yield processed, kmer, count

def get_counts(kmers, count_scaling):
    for processed, kmer, count in kmers:
        if count_scaling == 'counts':
            yield kmer, count
        elif count_scaling == 'binary':
            yield kmer, 1
        else: 
            yield kmer, np.log1p(count) 
        
def hash_features(kmers, n_features):
    extractor = FeatureHasher(n_features = n_features,
                              input_type = "pair",
                              non_negative=True)
    features = extractor.transform([kmers]).toarray()
    return features
        
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--rand-proj-fl",
                        type=str,
                        required=False)

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

    parser.add_argument('--feature-scaling-before', choices=['binary', 'counts', 'log1p'], required=True)
   
    parser.add_argument('--feature-scaling-after', choices=['binary', 'counts', 'log1p'], required=True)
  
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    print "Opening kmer file"
    kmers = read_kmers(args.kmer_freq_fl)
  
    kmers = progress(kmers)
    kmers = get_counts(kmers, args.feature_scaling_before)

    n_features = 2 ** args.n_features
    features = hash_features(kmers, n_features=n_features)
     
    if args.feature_scaling_after == "binary":
        print("Using binary features after")
        features = Binarizer().fit_transform(features)

    elif args.feature_scaling_after == "log1p":
        print("Using log1p features after")
        features = np.log1p(features)
    else:
        print("Usng counts features after")
    
    dump(features, args.feature_matrix)
