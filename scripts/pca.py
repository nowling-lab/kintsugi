#!/usr/bin/env python

import argparse

from joblib import load
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

def load_matrix(flname):
    return load(flname)

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--feature-matrices",
                        nargs="+",
                        type=str)

    parser.add_argument("--sample-names",
                        nargs="+",
                        type=str)

    parser.add_argument("--plot-fl",
                        type=str,
                        required=True)

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    
    matrices = []
    labels = []
    for flname in args.feature_matrices:
        matrix = load_matrix(flname)
        for name in args.sample_names:
            if name in flname:
                matrices.append(matrix)
                labels.append(name)

    feature_matrix = sp.vstack(matrices).toarray()

    pca = PCA(n_components = 6, whiten=True)
    proj = pca.fit_transform(feature_matrix)

    labels = np.array(labels)
    plt.scatter(proj[:, 0], proj[:, 1], label=labels)
    plt.legend()
    plt.xlabel("Component 1", fontsize=18)
    plt.ylabel("Component 2", fontsize=18)

    plt.savefig(args.plot_fl, DPI=300)
