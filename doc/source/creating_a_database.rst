Creating a database
===================

A database is created using the :func:`gffutils.create_db` function.  This
function automatically detects the file format (GFF or GTF).  The only required
input argument is the filename.  By default, the resulting database will be the
input file with a ``.db`` appended; alternatively you can supply your own
filename.

If `verbose=True`, then progress will be reported to the terminal.

Creating a database then is simply::

    >>> import gffutils
    >>> gffutils.create_db(gfffn, dbfn, verbose=True)

This may take several minutes to complete, depending on the size and complexity
of the input file.

Preparing your file
-------------------
GTF or GFF files directly downloaded from data providers may not be formatted
as you would like.  For example, FlyBase GFF files come with the entire genome
as a FASTA entry at the end of the GFF file.  Additionally, FlyBase chromosome
names do not match those used by the UCSC genome browser (they lack the leading
"chr").  Finally, there may be feature types that are not of interest, which
would serve to increase the file size and lookup time of the resulting
database.

The :func:`clean_gff` function can be used to:

* filter out featuretypes that are not of interest (for example, remove
  "transposon_insertion_sites" from file to be used strictly for gene models)
* add "chr" to the beginning of chromosomes
* remove FASTA sequence from the end of the file
* ensure start/stop coordinates are > 0 and that start is always <= stop

After running :func:`clean_gff`, the resulting file can be used to create the
database.
