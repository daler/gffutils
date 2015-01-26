
.. _database-ids:

Database IDs
============
A primary key is a unique identifier used in a database.  When
importing a GFF or GTF file with :mod:`gffutils` into a database, each unique
feature in the file should have its own primary key.

Primary keys are important because they are used to retrieve information from
the database using dictionary syntax.  For example in the :ref:`introduction`,
the example file has a line like this::

    chr2L   FlyBase gene    7529    9484    .       +       .       ID=FBgn0031208;Name=CG11023;

By default, the primary key GFF3 features (like this one) is the "ID" field of
the attributes.  So the unique key used for the database for this feature is
`FBgn0031208`.  This means we can access the gene from the database like this::

    gene = db['FBgn0031208']

Now, imagine we wanted to get the 5' UTRs for the gene.  Looking at the
example file in the :ref:`introduction`, we actually do
have an ID.  So we could access it like this::

    utr = db['five_prime_UTR_FBgn0031208:1_737']

This is quite awkward to type.  Plus, using this method we would have to type
all the unique IDs for each of the UTRs we wanted!

In a sense, we only need one good "hook" into the database by meaningful IDs,
and then we can access the other features based on parents and children.  That
is, we could get *all* the 5'UTRs for the gene without knowing their individual
IDs like this::

    utrs = db.children('FBgn0031208', featuretype='five_prime_UTR')

This works because 1) we have a unique ID for the gene, 2) we have unique IDs
for each 5'UTR, and 3) the relationships in the database are constructed using
these unique IDs.

If your input GFF or GTF file is formatted in the canonical way, the default
settings should work fine.  The rest of this section details strategies for
instructing :mod:`gffutils` to use the most meaningful primary key for your
particular input file.


.. _id_spec:

`id_spec`
---------


.. seealso::

    Examples that show the use of `id_spec`:

    * :ref:`F3-unique-3.v2.gff` (uses `id_spec=:seqid:`)
    * :ref:`jgi_gff2.txt`
    * :ref:`ncbi_gff3.txt`
    * :ref:`wormbase_gff2.txt`

The `id_spec` (ID specification) kwarg determines how to extract information
from each line in order to construct a primary key for the feature.  It can
have several different forms -- None, string, list, dictionary, or callable.

:None:
    The primary key for each feature will be an auto-incremented version of the
    feature type (e.g., "gene_1", "gene_2", etc).

:string:
    Use the attribute value.

    For example, `id_spec="ID"`. The primary key for each feature will be the
    value of the "ID" attribute.  If this is not found, then an
    auto-incremented version of the feature type is used

:list:

    Use the first available attribute value from the list.

    For example, `id_spec=["ID", "Name"]`: the primary key for each feature
    will be the value of the "ID" attribute.  If no "ID" attribute is found,
    then use the value of the "Name" field.  If this is not found, then an
    auto-incremented version of the feature type is used.

:dict:

    Use different strategies according to the featuretype.

    For example, `id_spec={"gene": "Name", "mRNA": ["ID", "transcript_id"]}`:

    * For "gene" features, the primary key will be the value of the "Name" attribute.
      If this is not found, then use an auto-incremented version of "gene".
    * For "mRNA" features, if the "ID" attribute exists, then use it as the
      primary key will be the value of the "ID" attribute.  If not, then use
      the "transcript_id" value.  If this is not found, then use an
      auto-incremented version of "mRNA".
    * For any other feature, the primary key will be the an auto-incremented version fo
      the feature type.

:special string:

    Use a GFF field value (from the first 8 columns) rather than an attributes
    value.  Must be surrounded by `:`.  The options to use can be found in the
    list :obj:`gffutils.constants._gffkeys` [:-1].

    For example, `id_spec=":seqid:"`:  use the "seqid" field as the primary
    key.


:function:

    Apply a custom function (or other callable object), and use its return
    value as the primary key.

    The function must accept a single :class:`gffutils.Feature` object.  It can
    return one of the following:

    * None, in which case the behavior is the same as `id_spec=None`.
    * A special string starting with `autoincrement:X`, which will
      auto-increment based on the value of `X`.  That is, if a function returns
      `autoincrement:chr21`, then the primary key of the first feature will be
      `chr21_1`, the second will be `chr21_2`, and so on.
    * A string to be used as the primary key.


The default for GFF3 files is `id_spec="ID"`.  If a feature has an "ID"
attribute, it will be used for the primary key.  If not, then an
auto-incremented key, based on the featuretype, will be used.

The default for GTF files is `id_spec={"gene": "gene_id", "transcript":
"transcript_id"}`.  Even though "gene" and "transcript" features do not exist
in the original file, :mod:`gffutils` infers the gene and transcript boundaries
(as described in :ref:`gtf`, and will use this `id_spec` for those inferred
regions.


.. _transform:

`transform`
-----------

.. seealso::

    Examples that show the use of `transform`:

    * :ref:`ensembl_gtf.txt`
    * :ref:`glimmer_nokeyval.gff3`
    * :ref:`wormbase_gff2_alt.txt`
    * :ref:`wormbase_gff2.txt`

The `transform` kwarg is a function that accepts single
:class:`gffutils.Feature` object and that returns a (possibly modified)
:class:`gffutils.Feature` object.  It is used to modify, on-the-fly, items as
they are being imported into the database.  It is generally used for files that
don't fit the standard GFF3 or GTF specs.

One example use-case is that FlyBase GFF3 files do have have a leading "chr"
for the seqid GFF field.  If we wanted to add this to each feature as it is
imported into the database, then we could use the following function::

    def add_chr(d):
        d['seqid'] = "chr" + d['seqid']
        return d



`merge_strategy`
----------------

.. seealso::

    Examples that show the use of `merge_strategy`:

    * :ref:`c_elegans_WS199_shortened_gff.txt`
    * :ref:`mouse_extra_comma.gff3`

This parameter specifies the behavior when two items have an identical
primary key.

For example, consider the following attribute strings for two
consecutive lines.  Assume that `id_spec="ID"`, in which case these two
lines have the same primary key::

    ID="exon_1"; Parent="transcript_1";
    ID="exon_1"; Parent="transcript_2";


Using `merge_strategy="merge"`, then there will be a single entry in
the database for `"exon_1"`, but the attributes will be merged and only
unique values will be retained.  The new, edited feature will end up
looking like this::

   ID="exon_1"; Parent="transcript_1,transcript_2";  # db key: "exon_1"

Using `merge_strategy="create_unique"`, then the second entry will have
a unique, autoincremented primary key assigned to it, and both lines
will be in the database, accessible by two different keys::

    ID="exon_1"; Parent="transcript_1";  # database key: "exon_1"
    ID="exon_1"; Parent="transcript_2";  # database key: "exon_1_1"


Using `merge_strategy="error"`, a :class:`gffutils.DuplicateIDError`
exception will be raised.  This means you will have to edit the file
yourself to fix the duplicated IDs.

Using `merge_strategy="warning"`, a warning will be printed to the
logger, and the feature will be skipped.

