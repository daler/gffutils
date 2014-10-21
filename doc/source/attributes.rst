Attributes
==========
The last field of a GFF or GTF file contains attributes.  As described in
:ref:`dialects`, these can be inconsistently formatted, but we try to the best
we can.  Once the attributes have been parsed, they can be accessed via
:attr:`Feature.attributes` or using getitem syntax on the :class:`Feature`
itself.

An :class:`attributes.Attributes` object behaves much like a dictionary, except
that all of its values are stored internally as a list. By default, all
attribute values are returned as lists, even 1-item lists.  However, this can
be changed using the `constants.always_return_list` setting.

Let's get an example :class:`Attributes` object to work with, by parsing a GFF
line:

>>> f = gffutils.feature.feature_from_line(
... 'chr2L\tFlyBase\texon\t8193\t8589\t.\t+\t.\tID=exon1; Parent=FBtr0300689,FBtr0300690')

The :class:`attributes.Attributes` object is accessed like this:

>>> f.attributes
<gffutils.attributes.Attributes object at ...>

It behaves like a dictionary of lists:

>>> for i in sorted(f.attributes.items()):
...     print('{i[0]}: {i[1]}'.format(i=i))
ID: ['exon1']
Parent: ['FBtr0300689', 'FBtr0300690']

>>> f.attributes['ID']
['exon1']

Usually it's more convenient to access the attributes directly from the
feature, like this:

>>> f['ID']
['exon1']


We can add attributes, again directly from the feature:

>>> f['parent_type'] = 'mRNA'

.. _always_return_list:

By default, a list is always returned, even for 1-item lists:

>>> f['parent_type'] == f.attributes['parent_type'] == ['mRNA']
True

However, we can change this behavior like so:

>>> gffutils.constants.always_return_list = False

Now the single values are returned as strings rather than 1-item lists:

>>> for i in sorted(f.attributes.items()):
...     print('{i[0]}: {i[1]}'.format(i=i))
ID: exon1
Parent: ['FBtr0300689', 'FBtr0300690']
parent_type: mRNA

Reset back to the original behavior:

>>> gffutils.constants.always_return_list = True

