Importing data into a database
==============================
The core of :mod:`gffutils` is a `sqlite3 <http://www.sqlite.org/>`_ database.
If you're familiar with SQL, you can check out the :ref:`schema`.

If you have a well-formatted GFF or GTF file, simply create a database with the
:func:`gffutils.create_db` function::

    import gffutils
    db = gffutils.create_db(filename, database_filename)

This is a one-time operation. This will parse the file, infer the relationships
among the features in the file, and store the features and relationships in the
file `database_filename`.  Once it is complete, from now on you just have to
attach to the existing `database_filename` like this::

    db = gffutils.FeatureDB(database_filename)

You can now use the tools in :mod:`gffutils` to work with the data.

In practice, however, GFF and GTF files do not always exactly match
the official specification.  These difficult-to-work-with files are the *raison
d'être* for :mod:`gffutils`, which provides many different parameters for
customization and for handling special cases.

The goal of these different parameters is to ultimately assign unique IDs to
each feature in the file (which is used as the primary key for the
:ref:`featurestable`.  The section :ref:`database-ids` details these settings,
and the :ref:`examples` show all sorts of tricks for getting
improperly-formatted files to work with :mod:`gffutils`.
