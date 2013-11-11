Importing data into a database
==============================
The core of :mod:`gffutils` is a `sqlite3 <http://www.sqlite.org/>`_ database.
If you're familiar with SQL, you can check out the :ref:`schema`.


If you have a well-formatted GFF or GTF file, create a database with the
:func:`gffutils.create_db` function::

    import gffutils
    db = gffutils.create_db(filename, database_filename)

In practice, however, GFF and GTF files do not always exactly match
specification.  In this case, there are different parameters to adjust when
creating the database. What parameters you choose completely depends on your
particular file.

In the end, we need to assign unique IDs to each feature in the file (which is
used as the primary key for the :ref:`featurestable`, which can get tricky
depending on your GFF or GTF file.  The section :ref:`database-ids` details
these settings, and the :ref:`examples` show all sorts of tricks for getting
improperly-formatted files to work with :mod:`gffutils`.
