#!/bin/bash

set -eo pipefail
set -x

# Full test suite after conda-installing deps
source activate tmp$PY
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
nosetests -x --with-doctest -a '!slow'
conda install -y --file $HERE/../docs-requirements.txt
(cd doc && make clean && make doctest)

# Fresh environment, pip-installed from just-created sdist tarball
conda create -y -n new python=$PY
source activate new
python setup.py clean sdist
pip install dist/gffutils-*.tar.gz
python -c 'import gffutils; print(gffutils.__version__)'
source deactivate
