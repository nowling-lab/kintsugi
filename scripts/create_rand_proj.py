#!/usr/bin/env python

import argparse

from joblib import dump
from scipy.sparse import csr_matrix
from sklearn.random_projection import SparseRandomProjection
        
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--n-features",
                        type=int,
                        required=True,
                        help="base 2")

    parser.add_argument("--n-components",
                        type=int,
                        default="auto")

    parser.add_argument("--n-samples",
                        type=int,
                        required=True)

    parser.add_argument("--output-fl",
                        type=str,
                        required=True)

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    print "Creating sparse random projection"
    srp = SparseRandomProjection(n_components = args.n_components)
    n_features = 2 ** args.n_features
    empty_X = csr_matrix((args.n_samples, n_features))
    srp.fit(empty_X)

    print "Saving sparse random projection"
    dump(srp, args.output_fl)
        
    
