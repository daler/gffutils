.. currentmodule:: gffutils

Change log
==========
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
