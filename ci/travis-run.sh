#!/bin/bash

set -eo pipefail
set -x

# Full test suite after conda-installing deps
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
nosetests -x --with-doctest -a '!slow'
conda install -y --file $HERE/../docs-requirements.txt
(cd doc && make clean && make doctest)

# Fresh environment, pip-installed from just-created sdist tarball
conda create -y -n new python=$TRAVIS_PYTHON_VERSION
source activate new
python setup.py sdist
pip install dist/gffutils-*.tar.gz
python -c 'import gffutils; print(gffutils.__version__)'
source deactivate
