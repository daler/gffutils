#!/bin/bash

set -eo pipefail
set -x

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
nosetests -x --with-doctest -a '!slow'
conda install -y --file $HERE/../docs-requirements.txt
(cd doc && make clean && make doctest)
