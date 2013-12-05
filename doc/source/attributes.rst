Attributes
==========
The last field of a GFF or GTF file contains attributes.  As described in
:ref:`dialects`, these can be inconsistently formatted, but we try to the best
we can.  Once the attributes have been parsed, they can be accessed via
:attr:`Feature.attributes` or using getitem syntax on the :class:`Feature`
itself.

*All* attribute values are lists of unicode objects.
As an example, the attribute string::

    ID="gene1"; biotype="protein_coding", Aliases="g1, examplegene";

becomes the dictionary::

    'ID': [u'gene1'],
    'biotype': [u'protein_coding']
    'Aliases': [u'g1', u'examplegene']

We could access the biotype of the feature by::

    feature['biotype'][0]



