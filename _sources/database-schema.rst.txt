Database schema
===============
.. _schema:

Schema
------
The following is the schema used for :mod:`gffutils`. Feel free to skip this if
you're not familiar with SQL. An explanation of each table can be found below.

>>> print(gffutils.constants.SCHEMA)  #doctest:+NORMALIZE_WHITESPACE
<BLANKLINE>
<BLANKLINE>
CREATE TABLE features (
    id text,
    seqid text,
    source text,
    featuretype text,
    start int,
    end int,
    score text,
    strand text,
    frame text,
    attributes text,
    extra text,
    bin int,
    primary key (id)
    );
<BLANKLINE>
CREATE TABLE relations (
    parent text,
    child text,
    level int,
    primary key (parent, child, level)
    );
<BLANKLINE>
CREATE TABLE meta (
    dialect text,
    version text
    );
<BLANKLINE>
CREATE TABLE directives (
    directive text
    );
<BLANKLINE>
CREATE TABLE autoincrements (
    base text,
    n int,
    primary key (base)
    );
<BLANKLINE>
CREATE TABLE duplicates (
    idspecid text,
    newid text,
    primary key (newid)
    );
<BLANKLINE>
<BLANKLINE>
<BLANKLINE>




.. _featurestable:

`features` table
----------------
The `features` table stores the primary information from each line in the
original file (and any additional information added by the user).

:id:

    Primary key for features.  The content of this field is determined by the
    user and the file format at creation time

    .. seealso::

        :doc:`database-ids` has more information about how the contents of this
        field are determined.

:seqid, source, feature, start, end, score, strand, frame:
    These fields correspond exactly to the fields in the GFF/GTF lines

:attributes:
    A JSON-serialized dictionary of attributes.  Note that the string
    representation of attributes is not stored; rather, it is reconstructed as
    needed using the dialect

    .. seealso::

        See :doc:`dialect` for what dialects are and how they are constructed

:extra:
    A JSON-serialized list of non-standard extra fields.  These are sometimes
    added by analysis tools (e.g., BEDTools).  For standard GFF/GTF, this
    field will be empty.

:bin:
    The genomic bin, according to the `UCSC binning
    <http://genome.cshlp.org/content/12/6/996/F7.expansion.html>`_ strategy.

`relations` table
-----------------
The `relations` table stores the heierarchical information.  It's sort of
a simple directed acyclic graph that seems to work well for GFF/GTF files with
[relatively] simple graph structure.

:parent:
    Foreign key to `features.id` -- a gene, for example.

:child:
    Foreign key to `feature.id` -- an mRNA or an exon, for example.

:level:
    In graph terms, the number of edges between `child` and `parent`.  In
    biological terms, if parent=gene and child=mRNA, then level=1.  If
    parent=gene and child=exon, then level=2.

`meta` table
------------
This table stores extra information about the database in general.

:dialect:
    A JSON-serialized version of the dialect empirically determined when
    parsing the original file.

    .. seealso::

        :doc:`dialect`

:version:
    The :mod:`gffutils` version used to create the database.

`directives` table
------------------
A table that acts as a simple list of directives (lines starting with `##`) in
the original GFF file.

:directive:
    String directive, without the leading `##`.

`autoincrements` table
----------------------
When items have conflicting primary keys based on the user-provided criteria
then :mod:`gffutils` can autoincrement in order to get a unique -- yet
reasonably meaningful -- primary key. For example, if the user specified that
the "ID" attributes field for a GFF3 file should be used for primary keys, but
two lines have the same `ID="GENE_A"` field, then the second line's ID will be
autoincremented to `ID="GENE_A_1"`.

After database creation, this table stores the autoincrementing information so
that when features are added later, autoincrementing can start at the correct
integer (rather than 0).

.. seealso::

    :doc:`database-ids`

:base:
    By default the feature type (`gene`, `exon`, etc) but can also be the value
    of any GFF field or attribute (e.g., the seqid or "GENE_1" (in the case
    of multiple features with ID="GENE_1").

:n:
    Current extent of autoincrementing -- add 1 to this when autoincrementing
    next time.

