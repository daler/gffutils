Developer's docs
================
This section serves as an entry point for learning about the internals of
:mod:`gffutils`.

Package modules
---------------

* `create.py` for creating new databases
* `parser.py` for parsing GFF/GTF files and determining the dialect
* `feature.py` contains the class definition for individual Feature objects
* `bins.py` implements the UCSC genomic binning strategy
* `constants.py` stores things like SELECT queries, GFF field names, database
  schema, and the default dialect
* `interface.py` provides the :class:`FeatureDB` class for interfacing with an
  existing db

General workflow
----------------
::

    Creation:

        Parser -> dictionary for each line -> imported into database using appropriate _DBCreator subclass

    Access:

        Interface with db using FeatureDB(database) -> Feature instances

        str(Feature) -> GFF/GTF line


Parsing
~~~~~~~
The :class:`parser.Parser` class handles the conversion of GFF/GTF lines into
dictionaries.  The :func:`parser._split_keyvals` function handles the
conversion of attribute strings -- field 9 of a GFF/GTF -- into dictionaries.  The
dictionaries are actually :class:`helpers.DefaultOrderedDict` instances, which
combine features of :class:`DefaultDict` and :class:`OrderedDict`.  This way,
the order of attribute key/val pairs is preserved upon reconstruction.  This
same :func:`parser._split.keyvals` function also inspects the format of the
provided attribute string to identify the "dialect" (see :ref:`dialects` for
more on this).

Import
~~~~~~
While the format of each line in a GFF and GTF file are *syntactically* similar
(same number of fields, roughly the same attributes string formatting), in the
context of the file as a whole they can be very *semantically* different.

For example, GTF files do not have "gene" features defined.  The genomic
coordinates of a gene must be inferred from the various "exon" features
comprising a particular gene.  For a GTF file, it's easy to figure out which
gene an exon belongs to by the "gene_id" attribute.

In contrast, GFF files typically have a "Parent" attribute.  For an exon, the
parent is the transcript; in order to figure out which gene an exon belongs to
requires looking at the parent of that transcript.  The transcript may be
defined many lines away in the GFF file, making it difficult to work with using
a line-by-line parsing approach.

The point of :mod:`gffutils` is to make access to the underlying data uniform
across both formats and to allow inter-conversion for use by downstream tools.
It does this by creating a relational database of features and parent-child
relationships.

Since the formats are so different, they require different methods of creation.
The :class:`create._DBCreator` class abstracts out common creation tasks.  The
:class:`create._GFFDBCreator` and :class:`create._GTFDBCreator` classes take
care of the format-specific routines.

:class:`_DBCreator` takes care of:
    * setting up the parser
    * logic for autoincrementing and handling primary keys
    * initializing the database
    * finalizing the db after format-specific tasks are complete -- things like
      writing version info, dialect, autoincrent info, etc.

:class:`_GFFDBCreator` and :class:`_GTFDBCreator` subclass :class:`_DBCreator`
and override the :meth:`_populate_from_lines` and :meth:`_update_relations`
methods.  Details are best left to the source code itself and the comments in
those methods.

The :func:`create.create_db` function delegates out to the appropriate class,
and all the docs for the kwargs are in this function.

Access
~~~~~~
Since the db creation imported the data into a uniform format, access requires
only a single class, :class:`interface.FeatureDB`.  Most methods on this class
simply perform queries on the database and return iterators of
:class:`feature.Feature` instances.
