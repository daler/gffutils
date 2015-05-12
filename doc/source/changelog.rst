.. currentmodule:: gffutils

Change log
==========
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
