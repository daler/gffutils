.. currentmodule:: gffutils

Change log
==========

Changes in v0.10.1
------------------

- Fix issue with new merge routine (`#152
  <https://github.com/daler/gffutils/issues/152>`_)

Changes in v0.10
----------------

- Support very large chromosomes (fixed issues `#94
  <https://github.com/daler/gffutils/issues/94>`_ and `#112
  <https://github.com/daler/gffutils/issues/112>`_)

- Expand `~` to user's home directory for filenames
  (issue `#105 <https://github.com/daler/gffutils/issues/105>`_).

- When merging, make merging attributes optional (issue `#107
  <https://github.com/daler/gffutils/issues/107>`_)

- Use a proper context manager for open files, fixes `issue #110
  <https://github.com/daler/gffutils/issues/110>`_.

- Update code to reflect changes in later Python versions (`#121
  <https://github.com/daler/gffutils/pull/121>`_ and  `#123
  <https://github.com/daler/gffutils/pull/123>`_ thanks @abhishekkumaresan)

- Dramatically improved merging routine -- many thanks to Nolan Wood
  @innovate-invent (`#130 <https://github.com/daler/gffutils/pull/130>`_).

- Previously, when merging the second feature's attributes were not
  deep-copied, resulting in unintended changes to the underlying dict (`#133
  <https://github.com/daler/gffutils/pull/133>`_, thanks Nolan Wood
  @innovate-invent)

- Fixed an issue that when imputing intron features, attributes were being
  pulled from the first (or last) exon (`#139
  <https://github.com/daler/gffutils/pull/139>`_, thanks @stekaz).

- Support creating Feature objects using empty values for attributes (`#144
  <https://github.com/daler/gffutils/pull/144>`_).

- Ensure that tests work post-installation (`#145
  <https://github.com/daler/gffutils/pull/145`_, thanks Michael Crusoe @mr-c)

- Removed redundant ``inspection.py`` module (`#147
  <https://github.com/daler/gffutils/pull/147>`_).

- Improvements to ``FeatureDB.update``, especially with respect to handling
  autoincrementing feature IDs. Previously, upon updating a db with another,
  autoincrement integers would restart at 1. Thanks Nolan Wood
  (@innovate-invent) and @abhishekkumaresan (`#149
  <https://github.com/daler/gffutils/pull/149>`_)


Changes in v0.9
---------------
Long-overdue release with performance improvements and better handling of
corner-case GFF and GTF files.

- performance tests (thanks Andrew Lando)
- performance improvements by building additional indexes (thanks Andrew Lando)
- performance improvments by running ``analyze features`` on created table
  (thanks Andrew Lando). Existing databases that have not had this run will
  trigger a warning suggesting that this should be run to speed up queries
  dramatically.
- add test for corner-case GTFs (`issue #79 <https://github.com/daler/gffutils/issues/79>`_)
- add fix for corner-case GFFs where `"="` is both a separator between fields
  as well as part of a value inside a field even when not quoted (`issue #82
  <https://github.com/daler/gffutils/issues/82>`_)
- fix handling of strandedness in :meth:`gffutils.feature.Feature.sequence`
  (`issue #87 <https://github.com/daler/gffutils/issues/87>`_)
- fix handling of corner-case GFFs that are completely missing a start or end
  position (`issue #85 <https://github.com/daler/gffutils/issues/85>`_)
- improvements to test framework
- All percent-encoded characters are decoded upon parsing (regardless of if the
  GFF3 spec says they should have been encoded in the first place), and then
  re-encoded when converting the Feature to a string (`issue #98
  <https://github.com/daler/gffutils/issues/98>`_). Only characters specified
  in the GFF3 spec are re-encoded. Note that some GFF files have spaces encoded
  as `%20`, but spaces should not be encoded according to the GFF3 specs. In
  this case, they will be decoded into spaces upon parsing, but not re-encoded
  when converting to string. Set
  `gffutils.constants.ignore_url_escape_characters=True` to disable any
  encoding/decoding behavior.
- improved testing framework

Changes in v0.8.7.1
-------------------
Fixes bug in `gffutils.pybedtools_integration.tsses` where iterating over large
databases and using the `as_bed6=True` argument could cause a deadlock.

Changes in v0.8.7
-----------------
New module, :mod:`gffutils.pybedtools_integration`. In particular, the
:func:`gffutils.pybedtools_integration.tsses` function provides many options
for creating a GTF, GFF, or BED file of transcription start sites (TSSes) from
an annotation.

Changes in v0.8.6.1
-------------------
Only a warning -- and not an ImportError -- is raised if BioPython is not installed.

Lots of updates in the testing framework to use docker containers on travis-ci.org.

Changes in v0.8.4
-----------------
This version addresses issues `#48
<https://github.com/daler/gffutils/issues/48>`_ and `#20
<https://github.com/daler/gffutils/issues/20>`_. It only affects database
creation using certain GTF files.

To summarize, there are some publicly available GTF files that don't match the
GTF specification and have transcripts and genes already added. By default,
`gffutils` assumes a GTF matches spec and that there are no transcript or gene
features. It infers transcript and gene extents from exons alone. So for these
off-spec GTF files, `gffutils` would do a lot of extra work inferring the
transcript and gene extents, and then it would try to the inferred features
back into the database. Since they were already there, it triggered `gffutils`'
feature-merging machinery.

The point is, if you didn't specifically tell `gffutils` to skip this step, all
of this extra merging work would cause database creation to take far longer
than it should have (possibly 10-100x longer).

With v0.8.4, if you create a database out of a GTF file and there are
transcript or gene features in it, `gffutils` will emit a warning and
a recommendation to disable inferring transcripts and/or genes to speed things
up dramatically.

The new keyword arguments for controlling this in :func:`gffutils.create_db`
are `disable_infer_transcripts` and `disable_infer_genes`.  These are both set
to *False* by default.

The previous, soon-to-be-deprecated way of doing this was to use
`infer_gene_extent=False`.  The new equivalent is to use
`disable_infer_transcripts=True` and `disable_infer_genes=True`. If you use the
old method, it will be automatically converted to the new method and a warning
will be emitted.

This new behavior is more flexible since it gives us the ability to infer
transcripts if genes exist, or infer genes if transcripts exist (rather than
the previous all-or-nothing approach).


Changes in v0.8.3.1
-------------------
Thanks to Sven-Eric Schelhorn (@schellhorn on github), this version fixes a bug
where, if multiple gffutils processes try to create databases from GTF files
simultaneously, the resulting databases would be incomplete and incorrect.


Changes in v0.8.3
-----------------
New features
~~~~~~~~~~~~
- New :func:`inspect.inspect` function for examining the contents of
  a GFF or GTF file.

- New :meth:`Feature.sequence` method to extract the sequence for a feature (uses
  `pyfaidx <https://github.com/mdshw5/pyfaidx>`_).

- Expose `ignore_strand` kwarg in :meth:`FeatureDB.children_bp` method

- When creating or updating a database, the provided `transform` function can
  return a value evaluating to False which will cause that feature to be
  skipped.

- :func:`create_db()` can use remote gzipped files as input
- New :meth:`FeatureDB.delete` method to delete features from
  a database

- Initial support for BioPython SeqFeature objects

- `limit` kwarg can now be used for :meth:`FeatureDB.parents` and
  :meth:`FeatureDB.children` to restrict returned features to
  a genomic range

- :meth:`FeatureDB.interfeatures` can now update attributes

- Much more flexible :meth:`FeatureDB.region` that allows slice-like
  operations.

- Improve :meth:`FeatureDB.update` so that entire features (rather
  than just attributes) can be replaced or updated (thanks Rintze Zelle for
  ideas and testing)

Bug fixes
~~~~~~~~~
- fix a bug when using a function as an `id_spec` for `create_db()` function
  (thanks @moritzbuck on github)
