.. _examples:

Examples
--------


Often, GTF/GFF files have non-standard formatting.  While :mod:`gffutils` is
set up to handle standard formatting by default, it also allows great
flexibility for working with non-standard files.  These examples collect
various real-world, difficult-to-work-with files and show how to make them work
with :mod:`gffutils`.  If you have a troublesome file, the best thing to do is
to look at the first few lines for these examples and see if you can find one
with matching formatting.

These example files ship with :mod:`gffutils`.  For interactive use, their
absolute path can be obtained with the :func:`gffutils.example_filename`
function. Each example shows the contents of the file, along with text
explaining the troublesome parts and which settings to use in order to
successfully import into a database and retrieve features.


.. note::

    Be mindful of the difference between the ID attribute of a feature and the
    primary key of that feature in the database.

    :attr:`Feature.id` will always return the primary key.  `Feature["ID"]` will
    return the ID attribute of the feature, if it has one.  In many (but not
    all) cases, `Feature.id == Feature["ID"]`

.. note::

    Also be mindful of the difference between an attribute of a Python object
    vs the 9th field of a GFF line, which is called the "attributes" field
    according to the spec. Confusing?  Yes!


This documentation undergoes automated testing, and in order to support both
Python 2 and 3 in the same test suite, we need to do this:

>>> from __future__ import print_function

.. _mouse_extra_comma.gff3:

.. rst-class:: html-toggle

mouse_extra_comma.gff3
~~~~~~~~~~~~~~~~~~~~~~
This file has gene, mRNA, protein, UTR, CDS, and exon annotations. However, not
all of them have a useful ID in their attributes.


File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/mouse_extra_comma.gff3


Import
``````
The defaults for :mod:`gffutils` will load this example into a database.
However, the CDSs do not have ID attributes, and the default for GFF files is
`id_spec="ID"`.  As described at :ref:`id_spec`, they will be accessed by the
names `"CDS_1"`, `"CDS_2"`, and so on.

Note that using `id_spec=["ID", "Name"]` won't work, because the CDSs share the
same Name attribute, `CDS:NC_000083.5:LOC100040603`. We could use
`id_spec=["ID", "Name"]` in combination with `merge_strategy="create_unique"`,
which would give ids like `CDS:NC_000083.5:LOC100040603_1` (with an underscore
an integer at the end for each consecutive occurence of the same ID).  But
we'll go with the simpler default strategy for this file:

>>> import gffutils
>>> fn = gffutils.example_filename('mouse_extra_comma.gff3')
>>> db = gffutils.create_db(fn, ':memory:')

Example access
``````````````
Since the mRNA on line 2 has an ID, we can simply access it by ID:

>>> db['XM_001475631.1']
<Feature mRNA (chr17:6797760-6818159[+]) at 0x...>

It's not very convenient to access the CDSs by their newly created ids
`"CDS_1"`, `"CDS_2"`, etc.  (Which one was "1"?  Is that from the same gene as
CDS "5"?) But we can access those CDSs as children of the mRNA:

>>> [f.id for f in db.children('XM_001475631.1', featuretype='CDS')]
['CDS_1', 'CDS_2', 'CDS_3', 'CDS_4', 'CDS_5']

Or as "grandchildren" of the gene:

>>> [f.id for f in db.children('NC_000083.5:LOC100040603', featuretype='CDS', level=2)]
['CDS_1', 'CDS_2', 'CDS_3', 'CDS_4', 'CDS_5']

Note that the database id values are all unique for the CDSs, but their
`"Name"` attributes are still all the same, as expected:

>>> list(set([f['Name'][0] for f in db.children('XM_001475631.1', featuretype='CDS')]))
['CDS:NC_000083.5:LOC100040603']

.. _c_elegans_WS199_ann_gff.txt:

c_elegans_WS199_ann_gff.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a minimal file, with only one feature and that one feature has no
genomic coordinates.  It illustrates how things behave when there is very
little information for a feature.


File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/c_elegans_WS199_ann_gff.txt

Import
``````
Import is straightforward:

>>> fn = gffutils.example_filename('c_elegans_WS199_ann_gff.txt')
>>> db = gffutils.create_db(fn, ':memory:')


Access
``````
Confirm there's only a single feature imported:

>>> list(db.all_features())
[<Feature experimental_result_region (I:.-.[+]) at 0x...]

With no obvious ID fields in the attributes, we can access the single feature
by name using the feature type:

>>> db['experimental_result_region_1']
<Feature experimental_result_region (I:.-.[+]) at 0x...>

This shows that the internal representation of null values is None:

>>> f = db['experimental_result_region_1']
>>> assert f.start is None
>>> assert f.stop is None

But the string representation shows "`.`" placeholders:

>>> print(f) #doctest: +NORMALIZE_WHITESPACE
I	Expr_profile	experimental_result_region	.	.	.	+	.	expr_profile=B0019.1


.. _c_elegans_WS199_shortened_gff.txt:

c_elegans_WS199_shortened_gff.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This file contains many different kinds of features; much like
:ref:`mouse_extra_comma.gff3` above, some features have IDs, some do not, and
there are many that share the same feature ID.  This example shows a common
way of subsetting features (here, getting all the non-historical CDSs for
a gene).

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/c_elegans_WS199_shortened_gff.txt


Import
``````
If `merge_strategy="merge"`, is used here, exceptions will be raised because
all those CDSs cannot be logically merged: they have different
coordinates. So `merge_strategy="create_unique"` is used instead. This kwarg
will add an underscore and an integer to each duplicate value used as a primary
key:

>>> fn = gffutils.example_filename('c_elegans_WS199_shortened_gff.txt')
>>> db = gffutils.create_db(fn, ":memory:", merge_strategy='create_unique', keep_order=True)

Access
``````
Getting a gene:

>>> g = db['Gene:WBGene00017003']
>>> g
<Feature gene (I:4580693-4583815[-]) at 0x...>

Get all "non-historical" CDSs for this gene.  This illustrates that:

#. you can use a :class:`Feature` object instead of a string for the database primary key
#. printing a :class:`Feature` reconstructs the line, and
#. the string representation of a line does not include a `\\n`.

>>> for i in db.children(g, featuretype='CDS', order_by='start'):
...     if i.source != 'history':
...         print(i)  # doctest:+NORMALIZE_WHITESPACE
I	Coding_transcript	CDS	4580993	4581241	.	-	0	ID=CDS:D1007.5a;Parent=Transcript:D1007.5a;status=Confirmed;wormpep=CE:CE29034
I	Coding_transcript	CDS	4581214	4581237	.	-	0	ID=CDS:D1007.5b;Parent=Transcript:D1007.5b.2,Transcript:D1007.5b.1;status=Confirmed;wormpep=WP:CE33577
I	Coding_transcript	CDS	4581664	4582026	.	-	0	ID=CDS:D1007.5a;Parent=Transcript:D1007.5a;status=Confirmed;wormpep=CE:CE29034
I	Coding_transcript	CDS	4581664	4582026	.	-	0	ID=CDS:D1007.5b;Parent=Transcript:D1007.5b.2,Transcript:D1007.5b.1;status=Confirmed;wormpep=WP:CE33577
I	Coding_transcript	CDS	4582412	4582718	.	-	1	ID=CDS:D1007.5a;Parent=Transcript:D1007.5a;status=Confirmed;wormpep=CE:CE29034
I	Coding_transcript	CDS	4582412	4582718	.	-	1	ID=CDS:D1007.5b;Parent=Transcript:D1007.5b.2,Transcript:D1007.5b.1;status=Confirmed;wormpep=WP:CE33577
I	Coding_transcript	CDS	4583190	4583374	.	-	0	ID=CDS:D1007.5a;Parent=Transcript:D1007.5a;status=Confirmed;wormpep=CE:CE29034
I	Coding_transcript	CDS	4583190	4583374	.	-	0	ID=CDS:D1007.5b;Parent=Transcript:D1007.5b.2,Transcript:D1007.5b.1;status=Confirmed;wormpep=WP:CE33577
I	Coding_transcript	CDS	4583426	4583509	.	-	0	ID=CDS:D1007.5a;Parent=Transcript:D1007.5a;status=Confirmed;wormpep=CE:CE29034
I	Coding_transcript	CDS	4583426	4583509	.	-	0	ID=CDS:D1007.5b;Parent=Transcript:D1007.5b.2,Transcript:D1007.5b.1;status=Confirmed;wormpep=WP:CE33577
I	Coding_transcript	CDS	4583560	4583805	.	-	0	ID=CDS:D1007.5a;Parent=Transcript:D1007.5a;status=Confirmed;wormpep=CE:CE29034
I	Coding_transcript	CDS	4583560	4583805	.	-	0	ID=CDS:D1007.5b;Parent=Transcript:D1007.5b.2,Transcript:D1007.5b.1;status=Confirmed;wormpep=WP:CE33577


.. _ensembl_gtf.txt:


ensembl_gtf.txt
~~~~~~~~~~~~~~~

In this GTF file, gene_id and transcript_id are identical, creating problems
for unique IDs in the database.


File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/ensembl_gtf.txt

Import
``````
As part of the import process, :mod:`gffutils` creates new "derived" features
for genes and transcripts (see :ref:`gtf`).  These are not explicitly defined in a GTF file, but
can be inferred.  Even though gene and transcript features do not exist in the
file, they will still be created.  This is why `id_spec={'gene': 'gene_id'}`
works: the primary key of "gene" features -- which don't exist in the GTF file
but are inferred by :mod:`gffutils` -- will be the "gene_id" attribute.

Note that in this example, gene_id and transcript_id are identical.  This
causes all sorts of problems.  To get around this we can write a transform
function (see :ref:`transform` for details) that applies an arbitrary
transformation to a dictionary returned by the parser.  This will modify the
data upon import into the database:

>>> def transform_func(x):
...     # adds some text to the end of transcript IDs
...     if 'transcript_id' in x.attributes:
...         x.attributes['transcript_id'][0] += '_transcript'
...     return x


Now we can supply this tranform function to :func:`create_db`:

>>> fn = gffutils.example_filename('ensembl_gtf.txt')
>>> db = gffutils.create_db(fn, ":memory:",
... id_spec={'gene': 'gene_id', 'transcript': "transcript_id"},
... merge_strategy="create_unique",
... transform=transform_func,
... keep_order=True)

Access
``````
We can access the derived "gene" feature by its ID (recall it doesn't actually
exist in the GTF file):

>>> g = db["B0019.1"]
>>> g
<Feature gene (I:12759579-12764949[-]) at 0x...>

>>> print(g)  # doctest:+NORMALIZE_WHITESPACE
I	gffutils_derived	gene	12759579	12764949	.	-	.	gene_id "B0019.1";

And the derived "transcript" by its ID -- which now has "_transcript" on the
end because of that transform function:

>>> t = db["B0019.1_transcript"]
>>> t
<Feature transcript (I:12759579-12764949[-]) at 0x...>

>>> print(t)  #doctest:+NORMALIZE_WHITESPACE
I	gffutils_derived	transcript	12759579	12764949	.	-	.	gene_id "B0019.1"; transcript_id "B0019.1_transcript";



.. _F3-unique-3.v2.gff:

F3-unique-3.v2.gff
~~~~~~~~~~~~~~~~~~
This file does not contain hierarchical information, and has many directives.
Since much of the benefit of :mod:`gffutils` comes from being able to navigate
the hierarchical relationships, it might be overkill to use :mod:`gffutils` on
this file.  An exeception would be if you want to be able to look up the read
names by their ID.  Here, we use yet a differnt way of overriding the way
:mod:`gffutils` decides on a primary key for the database.


File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/F3-unique-3.v2.gff

Import
``````
By default, the `id_spec="some_string"` parameter tells :mod:`gffutils` to
create a primary key based on the value of the "some_string" attribute.  For
this example, however, we would like to make the first field the primary key.
Luckily, :mod:`gffutils` will accept the name of a GFF field, surrounded by
":".  The names of the fields are in `gffutils.constants._gffkeys`:

>>> gffutils.constants._gffkeys
['seqid', 'source', 'featuretype', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

We want to use the first field, `seqid`, as the primary key for the
database, so we use `id_spec=":seqid:"`.

>>> fn = gffutils.example_filename('F3-unique-3.v2.gff')
>>> db = gffutils.create_db(fn, ':memory:', id_spec=[':seqid:'])

>>> db.directives[0]
'solid-gff-version 0.2'

Now we can access the features by their sequence name:

>>> db['3_8_425_F3']
<Feature read (3_8_425_F3:1119124-1119143[-]) at 0x...>


.. _glimmer_nokeyval.gff3:

.. rst-class:: html-toggle

glimmer_nokeyval.gff3
~~~~~~~~~~~~~~~~~~~~~
This is a problematic file in which the attributes are not all in format
"key=val".  This example demonstrates how :mod:`gffutils` handles cases like
this: by implicitly treating them as keys and setting the value to an empty
list.

Like the the :ref:`ensembl_gtf.txt` file, gene and transcript names are not
different.  Based on IDs alone, it's not clear whether the CDS's parent is the
gene or the mRNA.

However, the way to fix this is different from `ensembl_gtf.txt`.  Recall in
that case, we constructed a function to modify transcript IDs, and when
:mod:`gffutils` inferred the new transcripts from the component exons, it used
that new ID.

Here, we need to do several things in the transform function:  if the feature
type contains "RNA", then set the ID to `ID + "_transript"` (this will handle
things like ncRNA, tRNA, etc and is more general than just "mRNA").  Next, we
need to edit the Parent attribute of CDS, exon, etc.  How this is handled will
depend completely on the file you are using.


.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/glimmer_nokeyval.gff3

Import
``````
>>> def transform(f):
...     if f.featuretype == 'gene':
...         return f
...     elif "RNA" in f.featuretype:
...         f.attributes['ID'][0] += '_transcript'
...     else:
...         # assume that anything else is a child of a transcript, so we need
...         # to edit the "Parent" attribute
...         if 'Parent' in f.attributes:
...             f.attributes['Parent'] \
...                = [i + '_transcript' for i in f.attributes['Parent']]
...     return f
...

>>> fn = gffutils.example_filename('glimmer_nokeyval.gff3')
>>> db = gffutils.create_db(fn, ":memory:", id_spec="ID", transform=transform)

Access
``````

This shows that keys with missing values are assigned an empty list.  It's up
to the calling code to decide how to handle this (say, by checking for the
presence of the key "Complete".

>>> for k, v in sorted(db['GL0000007'].attributes.items()):
...     print(k, '=', v)
Complete = []
ID = ['GL0000007']
Name = ['GL0000007']

This illustrates that the CDS finds the proper transcript parent.


>>> for f in db.parents(db['CDS_1'], level=1):
...     print(f.id)
GL0000006_transcript

.. _hybrid1.gff3:

.. rst-class:: html-toggle

hybrid1.gff3
~~~~~~~~~~~~
This file contains a FASTA sequence at the end.  Currently, :mod:`gffutils`
assumes that there are no further annotations and stops parsing the file at
this point (specifically, when it sees the string "##FASTA" at the beginning of
a line).

Another issue with this file is that the "Note" attributes field contains
url-encoded characters. These are kept verbatim.


File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/hybrid1.gff3


Import
``````
Straightforward:

>>> fn = gffutils.example_filename('hybrid1.gff3')
>>> db = gffutils.create_db(fn, ":memory:")

Access
``````
Print out the attributes:

>>> for k, v in sorted(db['A00469'].attributes.items()):
...     print(k, '=', v)
Alias = ['GH1']
Dbxref = ['AFFX-U133:205840_x_at', 'Locuslink:2688', 'Genbank-mRNA:A00469', 'Swissprot:P01241', 'PFAM:PF00103', 'AFFX-U95:1332_f_at', 'Swissprot:SOMA_HUMAN']
ID = ['A00469']
Note = ['growth hormone 1']

.. _jgi_gff2.txt:


jgi_gff2.txt
~~~~~~~~~~~~
This file is GFF2 format -- the attribute format is like GTF, but without the
required "gene_id" and "transcript_id" fields that are used to infer genes and
transcripts. To get around this, we can supply the `gtf_transcript_key` kwarg,
which will override the default `"transcript_id"` (see :ref:`gtf` for more
background on this).

There is information in the file on which exon goes to which transcript, but no
such information about which CDS goes to which transcript.  There is, however,
information about the parent protein.

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/jgi_gff2.txt

Import
``````
Inspecting the file, we see that the "name" attribute is repeated on each line.
We can use that as the gene ID for inferring the gene extents, like "gene_id"
in a normal GTF file.  Similarly, the "transcriptId" attribute here can be used
like "transcript_id" in a normal GTF file for inferring transcript extents.

We also need to make a choice about how we're going to use this database.  Do
we want to have the CDS features be considered children of transcripts?  In
that case, we'll need the transform function, which creates a "transcriptId"
attribute out of an existing "proteinId" attribute:

>>> def transform(d):
...     try:
...         d['transcriptId'] = d['proteinId']
...     except KeyError:
...         pass
...     return d


>>> fn = gffutils.example_filename('jgi_gff2.txt')
>>> db = gffutils.create_db(fn, ":memory:",
... id_spec={'transcript': 'transcriptId', 'gene': 'name'},
... gtf_transcript_key='transcriptId', gtf_gene_key='name',
... transform=transform)

Access
``````
Since we used that transform function, the exons and CDSs are all children of
the "873" transcript:

>>> sorted(list(db.children('873', level=1)), key=lambda x: (x.featuretype, x.start))
[<Feature CDS (chr_1:37061-37174[-]) at 0x...>, <Feature CDS (chr_1:37315-37620[-]) at 0x...>, <Feature CDS (chr_1:37752-38216[-]) at 0x...>, <Feature exon (chr_1:37061-37174[-]) at 0x...>, <Feature exon (chr_1:37315-37620[-]) at 0x...>, <Feature exon (chr_1:37752-38216[-]) at 0x...>]

Here we can see that all children of the gene are accounted for:

>>> for f in db.children("fgenesh1_pg.C_chr_1000007", order_by='featuretype'):
...     print('{0.featuretype:>10}: {0.id}'.format(f))
       CDS: CDS_1
       CDS: CDS_2
       CDS: CDS_3
      exon: exon_1
      exon: exon_2
      exon: exon_3
transcript: 873


.. _ncbi_gff3.txt:

ncbi_gff3.txt
~~~~~~~~~~~~~
This file is problematic because all genes have the same ID attribute.  We can
get around this by using the "db_xref" field as the key for genes.


File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/ncbi_gff3.txt

Import
``````
For "gene" featuretypes, we want to use the "db_xref" field as the primary key;
all other featuretypes will be incremented versions of the featuretype
("CDS_1", "CDS_2", etc).

>>> fn = gffutils.example_filename('ncbi_gff3.txt')
>>> db = gffutils.create_db(fn, ":memory:", id_spec={'gene': 'db_xref'})

Access
``````
Genes can now be accessed by the db_xref:

>>> db['GeneID:4537201']
<Feature gene (NC_008596.1:12272-13301[+]) at 0x...>



.. _wormbase_gff2_alt.txt:


wormbase_gff2_alt.txt
~~~~~~~~~~~~~~~~~~~~~
This file has no gene_id or transcript_id fields; it appears to use "CDS" as
a gene-level kind of object.  So we can use a transform function to add
"gene_id" and "transcript_id" fields to all non-CDS features so that the file
conforms to a GTF standard and gene extents can be inferred.  Furthermore, by
default :mod:`gffutils` uses `'exon'` as the default feature type to merge into
genes.  Here, we need to specify `gtf_subfeature='coding_exon'`.


File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/wormbase_gff2_alt.txt

Import
``````
The transform function to manipulate the attribute dictionary: add "gene_id"
and "transcript_id" attributes to intron and coding_exon features:

>>> def transform(f):
...     if f.featuretype in ['coding_exon', 'intron']:
...         _id = f.attributes['CDS'][0]
...         f.attributes['gene_id'] = [_id]
...         f.attributes['transcript_id'] = [_id + '_transcript']
...     return f

>>> fn = gffutils.example_filename('wormbase_gff2_alt.txt')
>>> db = gffutils.create_db(fn, ":memory:", id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'}, transform=transform, gtf_subfeature='coding_exon')

Access
``````
Note that the inferred genes have a source of "gffutils_derived": 

>>> print(db["cr01.sctg102.wum.2.1"])  #doctest:+NORMALIZE_WHITESPACE
Contig102	gffutils_derived	gene	1629	3377	.	-	.	gene_id "cr01.sctg102.wum.2.1";


Get a report of the childrent of the gene:

>>> for f in db.children("cr01.sctg102.wum.2.1", order_by='featuretype'):
...     print('{0.featuretype:>12}: {0.id}'.format(f))
 coding_exon: coding_exon_1
 coding_exon: coding_exon_2
 coding_exon: coding_exon_3
 coding_exon: coding_exon_4
      intron: intron_1
      intron: intron_2
      intron: intron_3
  transcript: cr01.sctg102.wum.2.1_transcript


.. _wormbase_gff2.txt:


wormbase_gff2.txt
~~~~~~~~~~~~~~~~~
This file is problematic because of inconsistent formatting: the `SAGE_tag`
features have a different attributes format than the rest of the features.  So
we need to do some work with the dialect used for parsing attributes into
dictionaries.

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/wormbase_gff2.txt

Import
``````
We will need `force_dialect_check=True` in this case if we want to be able to have
the parents of the `SAGE_tag` features be recognized, because the attributes
for `SAGE_tag` features are formatted differently.  That is, they don't have
quoted values, and the separating semicolon does not have a space on either
side as is the case for the other features.

The default of `False` means that only the first few lines are checked to
determine the dialect to be used for all lines.  But that won't work here
because of the inconsistent formatting.  So we have to take the more
time-consuming approach of checking every dialect in order to figure out how to
parse them into a dictionary.


Transform function to create appropriate "Parent" attributes:

>>> def transform(f):
...     if f.featuretype == 'Transcript':
...         f.attributes['Parent'] = f.attributes['Gene']
...     else:
...         try:
...             f.attributes['Parent'] = f.attributes['Transcript']
...         except KeyError:
...             pass
...     return f

>>> fn = gffutils.example_filename('wormbase_gff2.txt')


>>> db = gffutils.create_db(fn,
... ":memory:",
... transform=transform,
... id_spec={'Transcript': "Transcript"},
... force_gff=True,
... force_dialect_check=True,
... keep_order=True)

The dialect for the database will be None:

>>> assert db.dialect is None, db.dialect

For cases like this, we should probably construct our own dialect to force all
attributes to have the same format.  To help with this, we can use the
:func:`helpers.infer_dialect` function by providing attributes:

>>> from gffutils import helpers
>>> dialect = helpers.infer_dialect([
... 'Transcript "B0019.1" ; WormPep "WP:CE40797" ; Note "amx-2" ; Prediction_status "Partially_confirmed" ; Gene "WBGene00000138" ; CDS "B0019.1" ; WormPep "WP:CE40797" ; Note "amx-2" ; Prediction_status "Partially_confirmed" ; Gene "WBGene00000138"',
... ])

>>> db.dialect = dialect


Access
``````
>>> t = db["B0019.1"]

Since we've set the dialect for the database, any features returned from the
database should follow that dialect:

>>> print(t)  #doctest:+NORMALIZE_WHITESPACE
I	Coding_transcript	Transcript	12759582	12764949	.	-	.	Transcript "B0019.1" ; WormPep "WP:CE40797" ; WormPep "WP:CE40797" ; Note "amx-2" ; Note "amx-2" ; Prediction_status "Partially_confirmed" ; Prediction_status "Partially_confirmed" ; Gene "WBGene00000138" ; Gene "WBGene00000138" ; CDS "B0019.1" ; Parent "WBGene00000138" ; Parent "WBGene00000138"

>>> len(list(db.children(t, featuretype='exon')))
15

>>> db['SAGE_tag_1'].attributes['Transcript']
['B0019.1']


>>> len(list(db.children(t, featuretype='SAGE_tag')))
4

.. _gencode-v190.gtf:

gencode-v19.gtf
~~~~~~~~~~~~~~~
This example file contains the first gene from the Gencode v19 annotations.
It's in GTF format, but in contrast to an on-spec GTF file, which only has
exon, CDS, UTRs, and start/stop codon features, this file already has genes and
transcripts on separate lines.  This means we can avoid the gene and transcript
inference, which saves time.

.. literalinclude:: ../../gffutils/test/data/gencode-v19.gtf

Import
``````
We can use `disable_infer_genes=True` and `disable_infer_transcripts=True` to
disable the gene and transcript inference.  This will dramatically speed up
import, but will result in identical results:

>>> fn = gffutils.example_filename('gencode-v19.gtf')
>>> db = gffutils.create_db(fn,
... ":memory:",
... keep_order=True,
... disable_infer_genes=True, disable_infer_transcripts=True)



Access
``````
Access is the same as other files:

>>> db['ENSG00000223972.4']
<Feature gene (chr1:11869-14412[+]) at 0x...>


Ensure that inferring gene extent results in the same database -- it just takes
longer to import:

>>> db2 = gffutils.create_db(fn,
... ":memory:",
... keep_order=True,
... disable_infer_genes=True, disable_infer_transcripts=True)

>>> for i, j in zip(db.all_features(), db2.all_features()):
...     assert i == j
