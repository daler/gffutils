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
Three kinds of input data can be provided: a filename, an iterator of
:class:`Feature` objects, or a string containing GFF/GTF lines.

These are all provided to :class:`DataIterator`, which delegates out to the
right iterator class (:class:`FileIterator`, :class:`FeatureIterator`,
:class:`StringIterator`).  The iterator class must define
a :class:`_custom_iter` method, which is responsible for taking input in
whatever form (filename, Feature objects, string of lines, respectively) and
always yielding :class:`Feature` objects.

These iterator classes are subclasses of :class:`BaseIterator`.

Upon initialization, a :class:`BaseIterator` figures out what the dialect is.
This means consuming `checklines` :class:`Feature` objects from
:meth:`_custom_iter`.  It uses the :func:`iterators.peek` function to do this,
which returns a list of the first `checklines` objects, along with a new,
re-wound iterator.  This new iterator is stored as the
:attr:`BaseIterator._iter` attribute.

At this point, the :class:`BaseIterator` has a dialect.  Its :meth:`__iter__`
method iterates over the :meth:`_custom_iter`, Upon iterating, it
will add this dialect to every :class:`Feature`.  This means that, no matter
what format `data` is, the following will print the lines correctly::

    >>> for feature in DataIterator(data):
    ...     print feature

A dialect can be optionally provided, which will disable the automatic dialect
inference.  This makes it straightforward to sanitize input, or convert to
a new dialect.  For example, to convert from GTF to GFF dialects::

    >>> for feature in DataIterator(GTF_data, dialect=GFF_dialect):
    ...     print feature

If `dialect` is not None, then that dialect will be used; otherwise, it will be
auto-detected.

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

The :class:`Feature` instances yielded from these iterators inherit the
database's dialect so that they print correctly.
