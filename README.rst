`gffutils`
==========


Overview and motivation
-----------------------
`gffutils` is a Python package used for doing things with GFF files that are
too complicated for a simple ``awk`` or ``grep`` command line call.  For example,
to get a BED file of genes from a GFF file, you can use something simple like::

    grep "gene" chr2.gff | awk '{print $1,$4,$5}' > out.bed

But how would you use command line tools to get:

* all the exons of a gene
* exons that are found in all isoforms of a gene
* a BED file of 3' exons from genes longer than 5 kb
* average number of isoforms for genes on the plus strand

These more complex questions are actually quite easy to answer using
`gffutils`.  This is made possible by importing a GFF or GTF file into
a sqlite3 database -- a one-time process -- and then searching and manipulating
the contents of the database to do useful things.

.. note::

    An older, somewhat slower version can still be found at
    https://github.com/daler/GFFutils_old, but that version is no longer
    actively developed.

Installation
------------

::

    pip install gffutils

Or, after cloning the git repository::

    python setup.py develop

Creating a database
-------------------
A GFF or GTF file needs to be imported into a local, single-file database
befor using with `gffutils`.  This database is a standard sqlite3 database, so
it can also be used independently of `gffutils` in other downstream tasks.  It
can take several minutes to create (depending on your input file) but it's
a one-time procedure.


.. .. doctest::
    :hide:
::

    >>> import os
    >>> import gffutils
    >>> gff_fn = gffutils.example_filename('dm3-5-genes.gff3')
    >>> db_fn = os.path.basename(gff_fn) + '.db'
    >>> if os.path.exists(db_fn):
    ...    os.unlink(db_fn)


Here we get an example file that ships with `gffutils`, and create a database
filename named after it but in the current working directory: In practice you
would provide your own GFF or GTF file.

.. .. doctest::
::

    # example filenames
    >>> import os
    >>> import gffutils
    >>> gff_fn = gffutils.example_filename('dm3-5-genes.gff3')
    >>> db_fn = os.path.basename(gff_fn) + '.db'

Then, create the database:

.. .. doctest::
::

    >>> import gffutils
    >>> gffutils.create_db(gff_fn, db_fn)

You will see some output, something like this::

    Populating features table and first-order relations: 192 (99%)...done in 0.3s
    Updating 2nd-order relations...done in 0.8s
    Creating indexes...done in 2.8s

Now you're ready to use `gffutils`...




Using the database interactively
--------------------------------
To connect to the database, simply create a `FeatureDB` class by passing the
database filename:

.. doctest::

    >>> import gffutils
    >>> db_fn = 'dm3-5-genes.gff3.db'
    >>> db = gffutils.FeatureDB(db_fn)

For reference, the schema is:

.. doctest::

    >>> db.print_schema()
    CREATE TABLE features (
                                    id text,
                                    chrom text,
                                    start int,
                                    stop int,
                                    strand text,
                                    featuretype text,
                                    score float,
                                    source text,
                                    frame text,
                                    attributes text,
                                    primary key (id)
                                  )
    CREATE TABLE relations (
                parent text,
                child text,
                level int,
                primary key(parent,child,level) )
    CREATE TABLE meta (
                filetype text
                )


Note that most `FeatureDB` methods return iterators for performance.  Show the
featuretypes that were in the GFF file:

.. doctest::

    >>> # Which kinds of featuretypes are in the database?
    >>> print list(db.featuretypes())
    ['CDS', 'exon', 'five_prime_UTR', 'gene', 'intron', 'mRNA', 'three_prime_UTR']


Return an iterator of all the features of one type:

.. doctest::

    >>> genes = db.features_of_type('gene')
    >>> gene = genes.next()
    >>> type(gene)
    <type 'gffutils.gfffeature.Feature'>

`Feature` objects have attributes like:

.. doctest::

    >>> gene.chrom
    '2L'

    >>> gene.start
    114726

    >>> gene.stop
    156030

    >>> gene.featuretype
    'gene'

    >>> len(gene)
    41305

`Feature.attributes` is a dictionary-like object:

.. doctest::

    >>> gene.attributes.keys()
    ['Ontology_term', 'gbunit', 'derived_computed_cyto', 'Alias', 'Dbxref', 'ID', 'Name']


    >>> gene.attributes['Name']
    'CG11455'

The primary key in the database for a feature is the ``ID`` field.  So you can
access features by their ID directly if you know it:

.. doctest::

    >>> ID = gene.attributes['ID']

    >>> print ID
    FBgn0031228

    >>> assert db[ID] == gene

Instead of a string ID, you can also use the `Feature` object itself:

.. doctest::

    >>> assert db[gene] == gene

Printing a `Feature` prints the full GFF line:

.. doctest::
    :options: -ELLIPSIS, +NORMALIZE_WHITESPACE

    >>> print gene
    2L	FlyBase	gene	114726	156030	.	+	.	ID=FBgn0031228;Name=CG11455;Alias=NADH ubiquinone oxidoreductase 15 kDa,NADH:ubiquinone oxidoreductase 15 kDa subunit;Ontology_term=SO:0000010,SO:0000087,GO:0006120,GO:0003954,GO:0005747;Dbxref=FlyBase:FBan0011455,FlyBase_Annotation_IDs:CG11455,GB_protein:AAF51538,GB_protein:ACZ94135,GB_protein:ACZ94134,GB_protein:AAN10510,GB_protein:ACZ94133,GB:AI404167,GB:AY069186,GB_protein:AAL39331,GB:CZ476154,MITODROME:MTDROME11455,UniProt/TrEMBL:Q7K1C0,INTERPRO:IPR019342,OrthoDB5_Drosophila:EOG5GHZJ8,OrthoDB5_Diptera:EOG5BRW72,OrthoDB5_Insecta:EOG5WPZSP,OrthoDB5_Arthropoda:EOG5N5TDD,OrthoDB5_Metazoa:EOG5PCDK4,EntrezGene:33179,InterologFinder:33179,BIOGRID:59439,FlyAtlas:CG11455-RA,GenomeRNAi:33179;gbunit=AE014134;derived_computed_cyto=21B3-21B3

The major advantage of `gffutils` is the ability to navigate the hierarchy of
relationships.  The `FeatureDB.children()` and `FeatureDB.parents()` methods
are the workhorses for this.

By default, all child (and grandchild, etc) features will be returned using the
`FeatureDB.children()` method.

.. doctest::

    >>> len(list(db.children(gene)))
    73

Looks like a pretty complex gene:

.. doctest::

    >>> from collections import Counter
    >>> Counter(i.featuretype for i in db.children(gene))
    Counter({'intron': 17, 'five_prime_UTR': 16, 'exon': 13, 'CDS': 13, 'mRNA': 11, 'three_prime_UTR': 3})

We can restrict the children to only a selected featuretype:

.. doctest::

    >>> len(list(db.children(gene, featuretype='mRNA')))
    11


Are any of these exons constitutive (present in all isoforms)?

.. doctest::

    >>> # All isoforms for this gene
    >>> isoforms = set(i.id for i in db.children(gene, featuretype='mRNA'))

    >>> # Iterate through the child exons; if the exon's parent mRNAs are the
    >>> # same as all the isoforms for the gene, then it's consitutive.

    >>> constitutive = []
    >>> for exon in db.children(gene, featuretype='exon'):
    ...     parent_isoforms = set(i.id for i in db.parents(exon, featuretype='mRNA'))
    ...
    ...     if isoforms == parent_isoforms:
    ...         constitutive.append(exon.id)

    >>> constitutive
    ['FBgn0031228:13']

Inspect that exon:

.. doctest::
    :options: -ELLIPSIS, +NORMALIZE_WHITESPACE

    >>> exon = db['FBgn0031228:13']
    >>> print exon
    2L	FlyBase	exon	155858	156030	.	+	.	ID=FBgn0031228:13;Name=CG11455:13;Parent=FBtr0078117,FBtr0078118,FBtr0301886,FBtr0301887,FBtr0301888,FBtr0306542,FBtr0330638,FBtr0330639,FBtr0330640,FBtr0330641,FBtr0330642;parent_type=mRNA

    >>> len(exon.attributes['Parent'])
    11

Exonic bp of the gene:

.. doctest::
    :options: -ELLIPSIS, +NORMALIZE_WHITESPACE

    >>> # These exons overlap quite a bit; summing the length of all exons
    >>> # wouldn't make sense if we wanted to calculate RPKM or something
    >>> exons = list(db.children(gene, featuretype='exon'))
    >>> for exon in exons:
    ...     print exon.start, exon.stop
    114726 114991
    155089 155178
    155089 155784
    155250 155784
    155333 155410
    155333 155429
    155333 155784
    155466 155784
    155494 155784
    155546 155784
    155567 155784
    155638 155784
    155858 156030

    >>> # So we can merge them to get the total exonic bp for this gene:
    >>> merged_exons = list(db.merge_features(db.children(gene, featuretype='exon')))
    >>> for i in merged_exons:
    ...     print i
    2L	.	merged_exon	114726	114991	.	+	.	
    2L	.	merged_exon	155089	155784	.	+	.	
    2L	.	merged_exon	155858	156030	.	+	.

    >>> sum(len(i) for i in merged_exons)
    1135

Longest protein for this gene:

.. doctest::

    >>> lengths = {}
    >>> for isoform in db.children(gene, featuretype='mRNA'):
    ...     lengths[isoform.id] = sum(len(i) for i in db.children(isoform, featuretype='CDS'))
    >>> sorted(lengths.items(), key=lambda x: x[1], reverse=True)[0]
    ('FBtr0301887', 306)


Longest transcript for this gene:

.. doctest::

    >>> lengths = {}
    >>> for isoform in db.children(gene, featuretype='mRNA'):
    ...     lengths[isoform.id] = sum(len(i) for i in db.children(isoform, featuretype='exon'))
    >>> sorted(lengths.items(), key=lambda x: x[1], reverse=True)[0]
    ('FBtr0330641', 869)


Gene in the database with the most exons:

.. doctest::

    >>> gene_with_most, exon_count = None, 0
    >>> for g in db.features_of_type('gene'):
    ...     this_count = sum(1 for _ in db.children(g, featuretype='exon'))
    ...     if this_count > exon_count:
    ...         gene_with_most = g
    ...         exon_count = this_count
    >>> gene_with_most.id, exon_count
    ('FBgn0031220', 18)
