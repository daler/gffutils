#!/bin/bash

# Tests the installation of requirements via bioconda channel and via pip for
# both Python 2 and 3.
#
# Runs main tests and doctests for Python 2 and 3.

set -e
set -x

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for PYTHON in 2 3; do

    # clone into a separate directory just for this python version
    src=/tmp/gffutils_py$PYTHON
    mkdir $src
    git clone $HERE/.. $src
    cd $src

    # extract version
    VERSION=$(python -c 'exec(open("gffutils/version.py").read());print(version)')

    # ------------------------------------------------------------------------
    # Installs everything from bioconda channel for this python version.
    # Don't run any tests, just make sure things can be installed using conda.
    conda create -n bioconda_py${PYTHON} python=$PYTHON
    set +x; source activate bioconda_py${PYTHON}; set -x
    conda install -c bioconda --file requirements.txt
    python setup.py develop
    (cd / && python -c 'import gffutils')
    set +x; source deactivate; set -x


    # ------------------------------------------------------------------------
    # Installs everything from sdist (simulates from PyPI)
    conda create -n py$PYTHON python=$PYTHON
    set +x; source activate py$PYTHON; set -x
    python setup.py clean sdist
    pip install dist/gffutils-${VERSION}.tar.gz
    (cd / && python -c 'import gffutils')


    # Prepare for testing by installing just nose for main tests, then run 'em
    conda install nose
    nosetests -x --with-doctest

    # Install tools and run doctests
    pip install -r docs-requirements.txt
    (cd doc && make clean && make doctest)
    set +x; source deactivate; set -x

done
