Installation
============

Using `conda <http://conda.pydata.org/docs/index.html>`_::

    conda install --channel bioconda gffutils

Or `pip`::

    pip install gffutils


`gffutils` is `tested <https://travis-ci.org/daler/gffutils>`_ with Python 2.7
and Python 3.3.

Optional requirements
---------------------

* BioPython (for creating SeqFeatures and SeqRecords)
* pybedtools (for integration with pybedtools)
* BEDTools (used by pybedtools)

Install them all with `conda`::

    conda install --channel bioconda pybedtools bedtools biopython
