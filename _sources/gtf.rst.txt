.. _gtf:

GTF files
=========
From the `GTF definition <http://mblab.wustl.edu/GTF22.html>`_:

    *The following feature types are required: 'CDS', 'start_codon',
    'stop_codon'. The features '5UTR', '3UTR', 'inter', 'inter_CNS',
    'intron_CNS' and 'exon' are optional. All other features will be ignored.*

So genes and transcript are not explicitly defined.  The transcript extent, is,
after all, implied by all exons with a single `transcript_id`, and the extent
of a gene is implied by all exons with a single `gene_id`.  However, this can
be tedious to calculate by hand.

:mod:`gffutils` **infers the gene and transcript extents when the file is
imported into a database**, and adds new "derived" features for each gene and
transcript.  That way, a gene can be easily accessed by its ID, just like for
GFF files.

However, not all files meet the official GTF specifications where each feature
has `transcript_id` to indicate its parent feature and `gene_id` to indicate
its "grandparent" feature.  To accommodate this, :mod:`gffutils` provides some
extra options for the :func:`gffutils.create_db` function:

`transcript_key` and `gene_key`
-------------------------------

.. seealso::

    Example that shows the use of `transcript_key`:

    * :ref:`jgi_gff2.txt`

These kwargs are used to extract the parent and grandparent feature
respectively.

By default, `transcript_key="transcript_id"` and `gene_key="gene_id"`.  But if
your particular data file does not conform to this, then they can be changed.

`subfeature`
------------
.. seealso::

    Examples that show the use of `subfeature`:

    * :ref:`wormbase_gff2_alt.txt`

Genes and transcripts are inferred from their component "exon" features.  If
your particular data file does not conform to the GTF standard, you can use the
`subfeature` kwarg to change this.  By default, `subfeature="exon"`, but see
the example :ref:`wormbase_gff2_alt.txt` for an instance of where this needs to
be changed.

