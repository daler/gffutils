#!/bin/bash

HERE=`pwd`
echo "-------------------------------------------------"
echo " Running nosetests"
echo "-------------------------------------------------"
nosetests -x --with-doctest

echo "-------------------------------------------------"
echo " Running doctests in Sphinx docs and source"
echo "-------------------------------------------------"
cd doc && make doctest && cd $HERE
