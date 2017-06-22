#!/bin/bash
set -euo pipefail
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/anaconda
export PATH=~/anaconda/bin:$PATH
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels bioconda
~/anaconda/bin/python setup.py install
