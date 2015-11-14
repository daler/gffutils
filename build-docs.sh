#!/bin/bash

# Build docs here, then copy them over to a fresh, temporary checkout of the
# gh-pages branch from github. Then upload 'em. After a few minutes, you'll see
# the newly-generated docs at daler.github.io/gffutils.

# Ideas from:
# http://executableopinions.readthedocs.org/en/latest/labs/gh-pages/gh-pages.html
set -e
set -x

(cd doc && make html)
HERE=$(pwd)
MSG="Adding gh-pages docs for $(git log --abbrev-commit | head -n1)"
DOCSOURCE=$HERE/doc/build/html
TMPREPO=/tmp/docs
rm -rf $TMPREPO
mkdir -p -m 0755 $TMPREPO
git clone git@github.com:daler/gffutils.git $TMPREPO
cd $TMPREPO
git checkout gh-pages
cp -r $DOCSOURCE/* $TMPREPO
touch $TMPREPO/.nojekyll
git add -A
git commit -m "$MSG"
git push origin gh-pages
cd $HERE
