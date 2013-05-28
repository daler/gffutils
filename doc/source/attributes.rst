Attributes
==========
The last field of a GFF or GTF file contains attributes.  As described in
:ref:`dialects`, these can be inconsistently formatted, but we try to the best
we can.  Once the attributes have been parsed, they can be accessed via
:attr:`Feature.attributes`.

:attr:`Feature.attributes` behaves like a combination of Python's
:class:`collections.defaultdict` and :class:`collections.OrderedDict`.  The
keys in the dictionary are ordered as they were in the attributes field from
the original line, and *all* values are lists.  Even single-value attributes.

As an example, the attribute string::

    ID="gene1"; biotype="protein_coding", Aliases="g1, examplegene";

becomes the dictionary::

    'ID': ['gene1'],
    'biotype': ['protein_coding']
    'Aliases': ['g1', 'examplegene']
