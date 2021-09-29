"""
Copyright 2021 Ronald J. Nowling

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from setuptools import find_packages
from setuptools import setup

setup(name="kintsugi",
      version=0.1,
      description="Genotype classification from k-mers",
      author="Ronald J. Nowling",
      author_email="rnowling@gmail.com",
      license="Apache License, Version 2.0",
      zip_safe=False,
      python_requires=">=3.6",
      install_requires=["numpy>=0.19.1", "scipy>=0.19.1", "sklearn", "mmh3", "trashcompactor", "snakemake", "pyyaml", "lz4"],
      packages=find_packages(include=["kintsugi", "kintsugi.*"]),
      package_data={
          "kintsugi" : ["snakefiles/*.smk"]
      },
      scripts=[
          "bin/evaluate_predictions",
          "bin/merge_kmer_count_partitions",
          "bin/partition_kmer_counts",
          "bin/predict_inversion_genotypes",
          "bin/simulate_kmer_data",
          "bin/split_data_set",
          "bin/kintsugi_cli",
          "bin/train_on_significant",
          "bin/kmer_association_testing"
      ]
)
