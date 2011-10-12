gffutils
========
.. contents::

.. note:

    An older, somewhat slower version can still be found at
    https://github.com/daler/GFFutils_old, but that version is no longer
    actively developed.

Overview and motivation
-----------------------
This module is used for doing things with GFF files that are too
complicated for a simple ``awk`` or ``grep`` command line call.

For example, to get a BED file of genes from a GFF file, you can use something
simple like::

    grep "gene" chr2.gff | awk '{print $1,$4,$5}' > out.bed

But how would you use commandline tools to get:

* all the exons of a gene
* exons that are found in all isoforms of a gene
* a BED file of 3' exons from genes longer than 5 kb
* average number of isoforms for genes on the plus strand

These more complex questions are actually quite easy to answer using
``gffutils``.
