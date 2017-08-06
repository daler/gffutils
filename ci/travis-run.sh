#!/bin/bash

set -eo pipefail
set -x

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
nosetests -x --with-doctest -a '!slow'
conda install -y --file $HERE/../docs-requirements.txt
(cd doc && make clean && make doctest)


conda create -n new python=$TRAVIS_PYTHON_VERSION
source activate new
python setup.py sdist
pip install dist/gffutils-*.tar.gz
python -c 'import gffutils; print(gffutils.__version__)'
source deactivate
