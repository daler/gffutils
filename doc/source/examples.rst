.. _examples:

Examples
--------
These example files ship with :mod:`gffutils`.  For interactive use, their
absolute path can be obtained with the :func:`gffutils.example_filename`
function.

Each section has a **"File contents"** subsection that contains the contents of the
example file but is hidden by default.

.. _mouse_extra_comma.gff3:

.. rst-class:: html-toggle

`mouse_extra_comma.gff3`
~~~~~~~~~~~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/mouse_extra_comma.gff3
    :lines: 1

This file has gene, mRNA, protein, UTR, CDS, and exon annotations. However, not
all of them have a useful ID in their attributes.




.. rst-class:: html-toggle

File contents
`````````````


.. literalinclude:: ../../gffutils/test/data/mouse_extra_comma.gff3

Import
``````

Imagine if we were to use `id_spec=["ID", "Name"]` and
`merge_strategy="merge"`.  :mod:`gffutils` would attempt to merge those CDS
features -- which all share the same name -- but would raise an exception
because their non-attribute fields (e.g., genomic coords) don't match.

Instead, we could use `merge_strategy="create_unique"`, which would create
unique IDs for those CDSs, appending a `_1`, `_2`, `_3`, and `_4` for each.

Another alternative would be to use `id_spec="ID"`, in which case only the gene
and mRNA features would be easily accessible.

>>> import gffutils
>>> fn = gffutils.example_filename('mouse_extra_comma.gff3')
>>> db = gffutils.create_db(fn, ':memory:', id_spec=['ID', 'Name'], merge_strategy='create_unique')

Example access
``````````````

Since the mRNA on line 2 has an ID, we can access it by its name:

>>> db['XM_001475631.1']
<Feature mRNA (chr17:6797760-6818159[+]) at 0x...>

However, the protein on line 3 does not have an ID.  This is a limitation of
how this particular GFF file is formatted.  :mod:`gffutils` falls back to using
an auto-increment of the feature type, "protein", so we can access it like
this:

>>> db['protein_1']
<Feature protein (chr17:6806527-6812289[+]) at 0x...>

For a large file, it would be difficult to know whether the protein you want is
`protein_1` or `protein_972` . . . luckily, it's still possible to get this
protein another way -- as a child of an mRNA with a known ID:

>>> db.children('XM_001475631.1', featuretype='protein').next()
<Feature protein (chr17:6806527-6812289[+]) at 0x...>

Or, as a level-2 child of the gene:
>>> db.children('NC_000083.5:LOC100040603', level=2, featuretype='protein')
<Feature protein (chr17:6806527-6812289[+]) at 0x...>

Access those CDSs either as children of the mRNA:

>>> [f.id for f in db.children('XM_001475631.1', featuretype='CDS')]
['CDS:NC_000083.5:LOC100040603', 'CDS:NC_000083.5:LOC100040603_1', 'CDS:NC_000083.5:LOC100040603_2', 'CDS:NC_000083.5:LOC100040603_3', 'CDS:NC_000083.5:LOC100040603_4']

Or get them individually.  This demonstrates that 1) the database ID is
different from the "Name" field in the attributes -- which doesn't change, and
2) the attributes are reconstructed.  Note the trailing comma, which was in the
original line, is represented by an extra item in the "Parent" attribute
list.

>>> cds4 = db['CDS:NC_000083.5:LOC100040603_4']
>>> for attr in cds4.attributes.items():
...     print attr
('Name', ['CDS:NC_000083.5:LOC100040603'])
('Parent', ['XM_001475631.1', ''])

>>> print cds4 # doctest:+NORMALIZE_WHITESPACE
chr17	RefSeq	CDS	6812219	6812289	.	+	2	Name=CDS:NC_000083.5:LOC100040603;Parent=XM_001475631.1,




.. _c_elegans__WS199_ann_gff.txt:

.. rst-class:: html-toggle

`c_elegans_WS199_ann_gff.txt`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/c_elegans_WS199_ann_gff.txt
    :lines: 2


This is a minimal file, with only one feature and that one feature has no
genomic coordinates.

.. rst-class:: html-toggle

File contents
`````````````
.. literalinclude:: ../../gffutils/test/data/c_elegans_WS199_ann_gff.txt

Import
``````
>>> fn = gffutils.example_filename('c_elegans_WS199_ann_gff.txt')
>>> db = gffutils.create_db(fn, ':memory:')


Access
``````
>>> list(db.all_features())
[<Feature experimental_result_region (I:.-.[+]) at 0x...]

With no obvious ID fields in the attributes, we can access it by name using the
feature type:

>>> db['experimental_result_region_1']
<Feature experimental_result_region (I:.-.[+]) at 0x...>


.. _c_elegans_WS199_shortened_gff.txt:

.. rst-class:: html-toggle

`c_elegans_WS199_shortened_gff.txt`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/c_elegans_WS199_shortened_gff.txt
    :lines: 7


This file contains many different kinds of features; much like
:ref:`mouse_extra_comma.gff3` above, some features have IDs, some do not, and
there are many that share the same feature ID.  This example shows a common
way of subsetting features (here, getting all the non-historical CDSs for
a gene).

.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/c_elegans_WS199_shortened_gff.txt

Import
``````
If `merge_strategy="merge"`, is used here, exceptions will be raised because
all those CDSs cannot be logically merged since they have different
coordinates. So `merge_strategy="create_unique"` is used instead.

>>> fn = gffutils.example_filename('c_elegans_WS199_shortened_gff.txt')
>>> db = gffutils.create_db(fn, ":memory:", id_spec="ID", merge_strategy='create_unique')

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

>>> for i in db.children(g, featuretype='CDS'):
...     if i.source != 'history':
...         print i  # doctest:+NORMALIZE_WHITESPACE
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

.. rst-class:: html-toggle

`ensembl_gtf.txt`
~~~~~~~~~~~~~~~~~

Example line:

.. literalinclude:: ../../gffutils/test/data/ensembl_gtf.txt
    :lines: 2

In this GTF file, gene_id and transcript_id are identical, creating problems
for unique IDs in the database.

.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/ensembl_gtf.txt

Import
``````
As part of the import process, :mod:`gffutils` creates new "derived" features
for genes and transcripts (see :ref:`gtf`).  These are not explicitly defined in a GTF file, but
can be inferred.  Even though gene and transcript features do not exist in the
file, they will still be created.  This is why `id_spec={'gene': 'gene_id'}`
works.

Note that in this example, gene_id and transcript_id are identical.  This
causes all sorts of problems; but to get around this we can write a transform
function that applies an arbitrary transformation to a dictionary returned by
the parser.  This will modify the data upon import into the database:

>>> def transform_func(x):
...     # adds some text to the end of transcript IDs
...     if 'transcript_id' in x['attributes']:
...         x['attributes']['transcript_id'][0] += '_transcript'
...     return x


Now we can supply this tranform function to :func:`create_db`:

>>> fn = gffutils.example_filename('ensembl_gtf.txt')
>>> db = gffutils.create_db(fn, ":memory:", id_spec={'gene': 'gene_id', 'transcript': "transcript_id"}, merge_strategy="create_unique", transform=transform_func)

Access
``````
We can access the derived "gene" feature by its ID:

>>> g = db["B0019.1"]
>>> g
<Feature gene (I:12759579-12764949[-]) at 0x...>

>>> print g  # doctest:+NORMALIZE_WHITESPACE
I	gffutils_derived	gene	12759579	12764949	.	-	.	gene_id "B0019.1";

And the derived "transcript" by its ID -- which now has "_transcript" on the
end because of that transform function:

>>> t = db["B0019.1_transcript"]
>>> t
<Feature transcript (I:12759579-12764949[-]) at 0x...>

>>> print t  #doctest:+NORMALIZE_WHITESPACE
I	gffutils_derived	transcript	12759579	12764949	.	-	.	transcript_id "B0019.1_transcript"; gene_id "B0019.1";

How many annotated transcripts are there for this gene?

>>> list(db.children(g, level=1))
[<Feature transcript (I:12759579-12764949[-]) at 0x...>]


.. _F3-unique-3.v2.gff:

.. rst-class:: html-toggle

`F3-unique-3.v2.gff`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/F3-unique-3.v2.gff
    :lines: 16

This file does not contain hierarchical information, and has many directives.
Since much of the benefit of :mod:`gffutils` comes from being able to navigate
the hierarchical relationships, it might be overkill to use :mod:`gffutils` on
this file.  An exeception would be if you want to be able to look up the read
names by their ID.


.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/F3-unique-3.v2.gff

Import
``````
By default, `id_spec` looks into the attributes field of a line to determine
what the ID of a feature should be.  For this example, however, we would like
to make the first field, "seqid", the primary key.  This is accomplished by
surrounding the key with `:`, so here, `id_spec=":seqid:"`.

>>> fn = gffutils.example_filename('F3-unique-3.v2.gff')
>>> db = gffutils.create_db(fn, ':memory:', id_spec=[':seqid:'])

>>> db.directives[0]
'solid-gff-version 0.2'

Now we can access the features by their sequence name:
>>> db['3_8_425_F3']
<Feature read (3_8_425_F3:1119124-1119143[-]) at 0x...>


.. _glimmer_nokeyval.gff3:

.. rst-class:: html-toggle

`glimmer_nokeyval.gff3`
~~~~~~~~~~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/glimmer_nokeyval.gff3
    :lines: 3

This is a problematic file in which the attributes are not all in format
"key=val".  This demonstrates how :mod:`gffutils` handles cases like this: by
implicitly treating them as keys and setting the value to an empty list.

Similar to the `ensembl_gtf.txt` file, gene and transcript names are not
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
...     if f['featuretype'] == 'gene':
...         return f
...     elif "RNA" in f['featuretype']:
...         f['attributes']['ID'][0] += '_transcript'
...     else:
...         if 'Parent' in f['attributes']:
...             f['attributes']['Parent'] \
...                = [i + '_transcript' for i in f['attributes']['Parent']]
...     return f
...

>>> fn = gffutils.example_filename('glimmer_nokeyval.gff3')
>>> db = gffutils.create_db(fn, ":memory:", id_spec="ID", transform=transform)

Access
``````

This shows that keys with missing values are assigned an empty list.  It's up
to the calling code to decide how to handle this (say, by checking for the
presence of the key "Complete".


>>> for k, v in db['GL0000007'].attributes.items():
...     print k, '=', v
ID = ['GL0000007']
Name = ['GL0000007']
Complete = []

This illustrates that the CDS finds the proper transcript parent.

>>> for f in db.parents(db['CDS_1'], level=1):
...     print f.id
GL0000006_transcript

.. _hybrid1.gff3:

.. rst-class:: html-toggle

`hybrid1.gff3`
~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/hybrid1.gff3
    :lines: 7

This file contains a FASTA sequence at the end.  Currently, :mod:`gffutils`
assumes that there are no further annotatins and stops parsing the file at this
point.  Also note that the "Note" attributes field contains url-encoded
characters. These are kept verbatim.

.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/hybrid1.gff3


Import
``````

>>> fn = gffutils.example_filename('hybrid1.gff3')
>>> db = gffutils.create_db(fn, ":memory:")

Access
``````

>>> for k, v in db['A00469'].attributes.items():
...     print k, '=', v
ID = ['A00469']
Dbxref = ['AFFX-U133:205840_x_at', 'Locuslink:2688', 'Genbank-mRNA:A00469', 'Swissprot:P01241', 'PFAM:PF00103', 'AFFX-U95:1332_f_at', 'Swissprot:SOMA_HUMAN']
Note = ['growth%20hormone%201']
Alias = ['GH1']

.. _jgi_gff2.txt:

.. rst-class:: html-toggle

`jgi_gff2.txt`
~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/jgi_gff2.txt
    :lines: 1

This file is GFF2 format -- the attribute format is like GTF, but without the
required "gene_id" and "transcript_id" fields.  To get around this, we can
supply the `gtf_transcript_key` kwarg, which will override the default
`"transcript_id"`.

.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/jgi_gff2.txt

Import
``````
Inspecting the file, we see that the "name" attribute is repeated on each line.
We can use that as the gene ID for inferring the gene extents.

>>> fn = gffutils.example_filename('jgi_gff2.txt')
>>> db = gffutils.create_db(fn, ":memory:",
... id_spec={'transcript': 'transcriptId', 'gene': 'name'},
... gtf_transcript_key='transcriptId', gtf_gene_key='name')

Access
``````
>>> list(db.children('873', level=1))
[<Feature exon (chr_1:37061-37174[-]) at 0x...>, <Feature exon (chr_1:37315-37620[-]) at 0x...>, <Feature exon (chr_1:37752-38216[-]) at 0x...>]

>>> for f in db.children("fgenesh1_pg.C_chr_1000007"):
...     print '{0.featuretype:>10}: {0.id}'.format(f)
transcript: 873
       CDS: CDS_1
      exon: exon_1
       CDS: CDS_2
      exon: exon_2
       CDS: CDS_3
      exon: exon_3


.. _ncbi_gff3.txt:

.. rst-class:: html-toggle

`ncbi_gff3.txt`
~~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/ncbi_gff3.txt
    :lines: 6

This file is problematic because all genes have the same ID attribute.  We can
get around this by using the "db_xref" field as the key for genes.

.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/ncbi_gff3.txt

Import
``````

>>> fn = gffutils.example_filename('ncbi_gff3.txt')
>>> db = gffutils.create_db(fn, ":memory:", id_spec={'gene': 'db_xref'})

Access
``````
>>> db['GeneID:4537201']
<Feature gene (NC_008596.1:12272-13301[+]) at 0x...>



.. _wormbase_gff2_alt.txt:

.. rst-class:: html-toggle

`wormbase_gff2_alt.txt`
~~~~~~~~~~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/wormbase_gff2_alt.txt
    :lines: 2

This file has no gene_id or transcript_id fields; it appears to use "CDS" as
a gene-level kind of object.  So we can use a transform function to add
"gene_id" and "transcript_id" fields to all non-CDS features so that the file
conforms to a GTF standard and gene extents can be inferred.  Furthermore, by
default :mod:`gffutils` uses `'exon'` as the default feature type to merge into
genes.  Here, we need to specify `gtf_subfeature='coding_exon'`.

.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/wormbase_gff2_alt.txt

Import
``````

>>> def transform(f):
...     if f['featuretype'] in ['coding_exon', 'intron']:
...         _id = f['attributes']['CDS'][0]
...         f['attributes']['gene_id'] = [_id]
...         f['attributes']['transcript_id'] = [_id + '_transcript']
...     return f

>>> fn = gffutils.example_filename('wormbase_gff2_alt.txt')
>>> db = gffutils.create_db(fn, ":memory:", id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'}, transform=transform, gtf_subfeature='coding_exon')

Access
``````

>>> print db["cr01.sctg102.wum.2.1"]  #doctest:+NORMALIZE_WHITESPACE
Contig102	gffutils_derived	gene	1629	3377	.	-	.	gene_id "cr01.sctg102.wum.2.1";

>>> for f in db.children("cr01.sctg102.wum.2.1"):
...     print '{0.featuretype:>12}: {0.id}'.format(f)
 coding_exon: coding_exon_4
  transcript: cr01.sctg102.wum.2.1_transcript
      intron: intron_3
 coding_exon: coding_exon_3
      intron: intron_2
 coding_exon: coding_exon_2
      intron: intron_1
 coding_exon: coding_exon_1


.. _wormbase_gff2.txt:

.. rst-class:: html-toggle

`wormbase_gff2.txt`
~~~~~~~~~~~~~~~~~~~
Example line:

.. literalinclude:: ../../gffutils/test/data/wormbase_gff2.txt
    :lines: 3

.. rst-class:: html-toggle

File contents
`````````````

.. literalinclude:: ../../gffutils/test/data/wormbase_gff2.txt

Import
``````
>>> def transform(f):
...     if f['featuretype'] == 'Transcript':
...         f['attributes']['Parent'] = f['attributes']['Gene']
...     else:
...         f['attributes']['Parent'] = f['attributes']['Transcript']
...     return f

>>> fn = gffutils.example_filename('wormbase_gff2.txt')
>>> db = gffutils.create_db(fn, ":memory:", transform=transform, id_spec={'Transcript': "Transcript"}, force_gff=True)

Access
``````
>>> t = db["B0019.1"]
>>> print t  #doctest:+NORMALIZE_WHITESPACE
I	Coding_transcript	Transcript	12759582	12764949	.	-	.	Transcript "B0019.1" ; WormPep "WP:CE40797" ; WormPep "WP:CE40797" ; Note "amx-2" ; Note "amx-2" ; Prediction_status "Partially_confirmed" ; Prediction_status "Partially_confirmed" ; Gene "WBGene00000138" ; Gene "WBGene00000138" ; CDS "B0019.1" ; Parent "WBGene00000138" ; Parent "WBGene00000138"

>>> len(list(db.children(t, featuretype='exon')))
15

>>> len(list(db.children(t, featuretype='SAGE_tag')))
4

