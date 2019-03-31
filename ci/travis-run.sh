#!/bin/bash

set -eo pipefail
set -x

# Full test suite after conda-installing deps
source activate tmp$PY
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
nosetests -x --with-doctest -a '!slow'

# explicitly test two versions of bedtools for the pybedtools_integration
# module

conda install -y "bedtools<2.27"
GFFUTILS_USES_BEDTOOLS_227_OR_LATER="false" nosetests --with-doctest gffutils/pybedtools_integration.py
conda install -y "bedtools>=2.27"
GFFUTILS_USES_BEDTOOLS_227_OR_LATER="true" nosetests --with-doctest gffutils/pybedtools_integration.py

conda install -y --file $HERE/../docs-requirements.txt
(cd doc && make clean && make doctest)

# Fresh environment, pip-installed from just-created sdist tarball
if [[ $PY == 3 ]]; then
  pyver="python>=3.7"
else
  pyver="python=2"
fi

conda create -y -n new python=$PY
source activate new
python setup.py clean sdist
pip install dist/gffutils-*.tar.gz
python -c 'import gffutils; print(gffutils.__version__)'
source deactivate
