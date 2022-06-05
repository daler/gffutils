Installation
============

Using `conda <http://conda.pydata.org/docs/index.html>`_::

    conda install --channel conda-forge --channel bioconda gffutils

Or `pip`::

    pip install gffutils


`gffutils` is `tested <https://github.com/daler/gffutils/actions>`_ with Python
3.6, 3.7, 3.8, 3.9.

Optional requirements
---------------------

* BioPython (for creating SeqFeatures and SeqRecords)
* pybedtools (for integration with pybedtools)
* BEDTools (used by pybedtools)

Install them all with `conda`::

    conda install --channel conda-forge --channel bioconda pybedtools bedtools biopython
