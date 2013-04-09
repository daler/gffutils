.. _parser-dictionaries:

Parser dictionaries
===================
When reading a GTF or GFF file, :class:`gffutils.Parser` instance yields
a dictionary for each line.

The keys in this dictionary correspond to the fields in the GFF3 spec, along
with an additional key, `extra`:

* `seqid`
* `source`
* `featuretype` (rather than "type", to avoid confusion with Python's built-in
  :func:`type` function)
* `start`
* `end`
* `score`
* `strand`
* `phase`
* `attributes`
* `extra`

All but `attributes` and `extra` are always strings -- including `start` and
`end`.

`attributes` is an :class:`OrderedDict` instance, which is populated in the
same order as the key/val pairs in the attributes string.

`extra` is a list of extra fields that extend past the standard 9 GFF fields,
again with everything as a string.


.. seealso::

    The `id_spec` and `transform` kwargs for :func:`gffutils.create_db` can be
    functions (or other callable objects) that accept a dictionary as returned
    by a :class:`gffutils.Parser` instance.

>>> import gffutils
>>> p = gffutils.Parser('chr1 . gene 1 100 . + . ID="gene1"', from_string=True)
>>> d = iter(p).next()
>>> for key, val in sorted(d.items()):
...     print key, ':', repr(val)
attributes : DefaultOrderedDict([('ID', ['gene1'])])
end : '100'
extra : []
featuretype : 'gene'
frame : '.'
score : '.'
seqid : 'chr1'
source : '.'
start : '1'
strand : '+'
