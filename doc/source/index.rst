.. currentmodule:: gffutils


.. _introduction:

Introduction
============
:mod:`gffutils` is a Python package for working with `GFF
<http://www.sanger.ac.uk/resources/software/gff/spec.htm>`_ and `GTF
<http://mblab.wustl.edu/GTF22.html>`_ formats in a hierarchical manner.

Below is a short demonstration of :mod:`gffutils`. For the full documentation,
see :doc:`contents`.

Example file
------------
Consider a brief example that comes with :mod:`gffutils`.  This GFF3 file
describes a gene with two isoforms, each annotated with UTRs, introns, exons,
and CDSs.

>>> import gffutils
>>> fn = gffutils.example_filename('intro_docs_example.gff')
>>> print open(fn).read() #doctest:+NORMALIZE_WHITESPACE
chr2L	FlyBase	gene	7529	9484	.	+	.	ID=FBgn0031208;Name=CG11023;
chr2L	FlyBase	mRNA	7529	9484	.	+	.	ID=FBtr0300689;Name=CG11023-RB;Parent=FBgn0031208;
chr2L	FlyBase	mRNA	7529	9484	.	+	.	ID=FBtr0300690;Name=CG11023-RC;Parent=FBgn0031208;
chr2L	FlyBase	five_prime_UTR	7529	7679	.	+	.	ID=five_prime_UTR_FBgn0031208:1_737;Name=CG11023-u5;Parent=FBtr0300689,FBtr0300690;
chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690;
chr2L	FlyBase	CDS	7680	8116	.	+	0	ID=CDS_FBgn0031208:1_737;Name=CG11023-cds;Parent=FBtr0300689,FBtr0300690;
chr2L	FlyBase	intron	8117	8192	.	+	.	ID=intron_FBgn0031208:1_FBgn0031208:2;Name=CG11023-in;Parent=FBtr0300690;
chr2L	FlyBase	intron	8117	8192	.	+	.	ID=intron_FBgn0031208:1_FBgn0031208:3;Name=CG11023-in;Parent=FBtr0300689;
chr2L	FlyBase	exon	8193	8589	.	+	.	Parent=FBtr0300690;
chr2L	FlyBase	exon	8193	9484	.	+	.	ID=FBgn0031208:3;Name=CG11023:3;Parent=FBtr0300689;
chr2L	FlyBase	CDS	8193	8610	.	+	2	ID=CDS_FBgn0031208:3_737;Name=CG11023-cds;Parent=FBtr0300689;
chr2L	FlyBase	CDS	8193	8589	.	+	2	ID=CDS_FBgn0031208:2_737;Name=CG11023-cds;Parent=FBtr0300690;
chr2L	FlyBase	exon	8668	9484	.	+	.	ID=FBgn0031208:4;Name=CG11023:4;Parent=FBtr0300690;
chr2L	FlyBase	CDS	8668	9276	.	+	0	ID=CDS_FBgn0031208:4_737;Name=CG11023-cds;Parent=FBtr0300690;
chr2L	FlyBase	intron	8590	8667	.	+	.	ID=intron_FBgn0031208:2_FBgn0031208:4;Name=CG11023-in;Parent=FBtr0300690;
chr2L	FlyBase	three_prime_UTR	8611	9484	.	+	.	ID=three_prime_UTR_FBgn0031208:3_737;Name=CG11023-u3;Parent=FBtr0300689;
chr2L	FlyBase	three_prime_UTR	9277	9484	.	+	.	ID=three_prime_UTR_FBgn0031208:4_737;Name=CG11023-u3;Parent=FBtr0300690;
<BLANKLINE>

Create the database
-------------------
:mod:`gffutils` must pre-process input files by importing into a local sqlite3
file-based database.  `force=True` will overwrite any existing database:

>>> db = gffutils.create_db(fn, dbfn='test.db', force=True)

Since real-world files don't always completely adhere to the format
specifications, :mod:`gffutils` provides lots of configuration options to
tailor the pre-processing to your particular file -- see :doc:`database-import`
for more details.

Use the database
----------------

Once the database has been created, downstream code can connect to it simply
by:

>>> import gffutils
>>> db = gffutils.FeatureDB('test.db')

Access features by name, which returns a :class:`Feature` instance:

>>> gene = db['FBgn0031208']
>>> gene
<Feature gene (chr2L:7529-9484[+]) at 0x...>

Inspect fields from the GFF line:

>>> gene.start
7529

>>> gene.end
9484

Access items within the attributes field (values are always in a list, even if
only one item):

>>> gene.attributes['Name']
['CG11023']

Printing a :class:`Feature` returns the original GFF line:

>>> print gene  #doctest:+NORMALIZE_WHITESPACE
chr2L	FlyBase	gene	7529	9484	.	+	.	ID=FBgn0031208;Name=CG11023;


Traversing the hierarchy returns a generator that can be iterated over.  Here
are the first-level children of the gene:

>>> for i in db.children(gene, level=1, order_by='start'):
...     print i #doctest:+NORMALIZE_WHITESPACE
chr2L	FlyBase	mRNA	7529	9484	.	+	.	ID=FBtr0300689;Name=CG11023-RB;Parent=FBgn0031208;
chr2L	FlyBase	mRNA	7529	9484	.	+	.	ID=FBtr0300690;Name=CG11023-RC;Parent=FBgn0031208;

Here are the second-level children, constrained to return only the exons:

>>> for i in db.children(gene, level=2, featuretype='exon', order_by='start'):
...     print i #doctest:+NORMALIZE_WHITESPACE
chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690;
chr2L	FlyBase	exon	8193	9484	.	+	.	ID=FBgn0031208:3;Name=CG11023:3;Parent=FBtr0300689;
chr2L	FlyBase	exon	8193	8589	.	+	.	Parent=FBtr0300690;
chr2L	FlyBase	exon	8668	9484	.	+	.	ID=FBgn0031208:4;Name=CG11023:4;Parent=FBtr0300690;

Given an exon, we can go up the hierarchy as well:

>>> for exon in db.features_of_type('exon'):
...     for g in db.parents(exon, featuretype='gene'):
...         assert g == gene

Retrieve entries by genomic coordinates, which uses the `UCSC binning strategy
<http://genome.cshlp.org/content/12/6/996/F7.expansion.html>`_:

>>> list(db.region(region=('chr2L', 9277, 10000), completely_within=True))
[<Feature three_prime_UTR (chr2L:9277-9484[+]) at 0x...>]

Or that overlap a feature:

>>> for UTR in db.region(gene, featuretype=['three_prime_UTR', 'five_prime_UTR']):
...     print UTR #doctest:+NORMALIZE_WHITESPACE
chr2L	FlyBase	five_prime_UTR	7529	7679	.	+	.	ID=five_prime_UTR_FBgn0031208:1_737;Name=CG11023-u5;Parent=FBtr0300689,FBtr0300690;
chr2L	FlyBase	three_prime_UTR	8611	9484	.	+	.	ID=three_prime_UTR_FBgn0031208:3_737;Name=CG11023-u3;Parent=FBtr0300689;
chr2L	FlyBase	three_prime_UTR	9277	9484	.	+	.	ID=three_prime_UTR_FBgn0031208:4_737;Name=CG11023-u3;Parent=FBtr0300690;


