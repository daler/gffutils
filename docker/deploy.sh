#!/bin/bash

HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
git clone $HERE/.. /tmp/gffutils
cd /tmp/gffutils
python setup.py sdist upload
