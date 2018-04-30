"""
Module for integration with pybedtools
"""

import os
import pybedtools
from pybedtools import featurefuncs
from gffutils import helpers
import six


def to_bedtool(iterator):
    """
    Convert any iterator into a pybedtools.BedTool object.

    Note that the supplied iterator is not consumed by this function. To save
    to a temp file or to a known location, use the `.saveas()` method of the
    returned BedTool object.
    """
    def gen():
        for i in iterator:
            yield helpers.asinterval(i)
    return pybedtools.BedTool(gen())


def tsses(db, merge_overlapping=False, attrs=None, attrs_sep=":",
          merge_kwargs=None, as_bed6=False, bedtools_227_or_later=True):
    """
    Create 1-bp transcription start sites for all transcripts in the database
    and return as a sorted pybedtools.BedTool object pointing to a temporary
    file.

    To save the file to a known location, use the `.moveto()` method on the
    resulting `pybedtools.BedTool` object.

    To extend regions upstream/downstream, see the `.slop()` method on the
    resulting `pybedtools.BedTool object`.

    Requires pybedtools.

    Parameters
    ----------
    db : gffutils.FeatureDB
        The database to use

    as_bed6 : bool
        If True, output file is in BED6 format; otherwise it remains in the
        GFF/GTF format and dialect of the file used to create the database.

        Note that the merge options below necessarily force `as_bed6=True`.

    merge_overlapping : bool
        If True, output will be in BED format. Overlapping TSSes will be merged
        into a single feature, and their names will be collapsed using
        `merge_sep` and placed in the new name field.

    merge_kwargs : dict
        If `merge_overlapping=True`, these keyword arguments are passed to
        pybedtools.BedTool.merge(), which are in turn sent to `bedtools merge`.
        The merge operates on a BED6 file which will have had the name field
        constructed as specified by other arguments here.  See the available
        options for your installed version of BEDTools; the defaults used here
        are `merge_kwargs=dict(o='distinct', c=4, s=True)`.

        Any provided `merge_kwargs` are used to *update* the default. It is
        recommended to not override `c=4` and `s=True`, otherwise the
        post-merge fixing may not work correctly. Good candidates for tweaking
        are `d` (merge distance), `o` (operation), `delim` (delimiter to use
        for collapse operations).

    attrs : str or list
        Only has an effect when `as_bed6=True` or `merge_overlapping=True`.

        Determines what goes in the name field of an output BED file. By
        default, "gene_id" for GTF databases and "ID" for GFF. If a list of
        attributes is supplied, e.g. ["gene_id", "transcript_id"], then these
        will be joined by `attr_join_sep` and then placed in the name field.

    attrs_sep: str
        If `as_bed6=True` or `merge_overlapping=True`, then use this character
        to separate attributes in the name field of the output BED. If also
        using `merge_overlapping=True`, you'll probably want this to be
        different than `merge_sep` in order to parse things out later.

    bedtools_227_or_later : bool
        In version 2.27, BEDTools changed the output for merge. By default,
        this function expects BEDTools version 2.27 or later, but set this to
        False to assume the older behavior.

        For testing purposes, the environment variable
        GFFUTILS_USES_BEDTOOLS_227_OR_LATER is set to either "true" or "false"
        and is used to override this argument.

    Examples
    --------

    >>> import gffutils
    >>> db = gffutils.create_db(
    ...    gffutils.example_filename('FBgn0031208.gtf'),
    ...    ":memory:",
    ...    keep_order=True,
    ...    verbose=False)

    Default settings -- no merging, and report a separate TSS on each line even
    if they overlap (as in the first two):


    >>> print(tsses(db))                        # doctest: +NORMALIZE_WHITESPACE
    chr2L	gffutils_derived	transcript_TSS	7529	7529	.	+	.	gene_id "FBgn0031208"; transcript_id "FBtr0300689";
    chr2L	gffutils_derived	transcript_TSS	7529	7529	.	+	.	gene_id "FBgn0031208"; transcript_id "FBtr0300690";
    chr2L	gffutils_derived	transcript_TSS	11000	11000	.	-	.	gene_id "Fk_gene_1"; transcript_id "transcript_Fk_gene_1";
    chr2L	gffutils_derived	transcript_TSS	12500	12500	.	-	.	gene_id "Fk_gene_2"; transcript_id "transcript_Fk_gene_2";
    <BLANKLINE>


    Default merging, showing the first two TSSes merged and reported as
    a single unique TSS for the gene. Note the conversion to BED:

    >>> x = tsses(db, merge_overlapping=True)
    >>> print(x)  # doctest: +NORMALIZE_WHITESPACE
    chr2L	7528	7529	FBgn0031208	.	+
    chr2L	10999	11000	Fk_gene_1	.	-
    chr2L	12499	12500	Fk_gene_2	.	-
    <BLANKLINE>

    Report both gene ID and transcript ID in the name. In some cases this can
    be easier to parse than the original GTF or GFF file. With no merging
    specified, we must add `as_bed6=True` to see the names in BED format.

    >>> x = tsses(db, attrs=['gene_id', 'transcript_id'], as_bed6=True)
    >>> print(x)  # doctest: +NORMALIZE_WHITESPACE
    chr2L	7528	7529	FBgn0031208:FBtr0300689	.	+
    chr2L	7528	7529	FBgn0031208:FBtr0300690	.	+
    chr2L	10999	11000	Fk_gene_1:transcript_Fk_gene_1	.	-
    chr2L	12499	12500	Fk_gene_2:transcript_Fk_gene_2	.	-
    <BLANKLINE>

    Use a 3kb merge distance so the last 2 features are merged together:

    >>> x = tsses(db, merge_overlapping=True, merge_kwargs=dict(d=3000))
    >>> print(x)  # doctest: +NORMALIZE_WHITESPACE
    chr2L	7528	7529	FBgn0031208	.	+
    chr2L	10999	12500	Fk_gene_1,Fk_gene_2	.	-
    <BLANKLINE>


    The set of unique TSSes for each gene, +1kb upstream and 500bp downstream:

    >>> x = tsses(db, merge_overlapping=True)
    >>> x = x.slop(l=1000, r=500, s=True, genome='dm3')
    >>> print(x)  # doctest: +NORMALIZE_WHITESPACE
    chr2L	6528	8029	FBgn0031208	.	+
    chr2L	10499	12000	Fk_gene_1	.	-
    chr2L	11999	13500	Fk_gene_2	.	-
    <BLANKLINE>


    """
    _override = os.environ.get('GFFUTILS_USES_BEDTOOLS_227_OR_LATER', None)
    if _override is not None:
        if _override == 'true':
            bedtools_227_or_later = True
        elif _override == 'false':
            bedtools_227_or_later = False
        else:
            raise ValueError(
                "Unknown value for GFFUTILS_USES_BEDTOOLS_227_OR_LATER "
                "environment variable: {0}".format(_override))

    if bedtools_227_or_later:
        _merge_kwargs = dict(o='distinct', s=True, c='4,5,6')
    else:
        _merge_kwargs = dict(o='distinct', s=True, c='4')

    if merge_kwargs is not None:
        _merge_kwargs.update(merge_kwargs)

    def gen():
        """
        Generator of pybedtools.Intervals representing TSSes.
        """
        for gene in db.features_of_type('gene'):
            for transcript in db.children(gene, level=1):
                if transcript.strand == '-':
                    transcript.start = transcript.stop
                else:
                    transcript.stop = transcript.start
                transcript.featuretype = transcript.featuretype + '_TSS'
                yield helpers.asinterval(transcript)

    # GFF/GTF format
    x = pybedtools.BedTool(gen()).sort()

    # Figure out default attrs to use, depending on the original format.
    if attrs is None:
        if db.dialect['fmt'] == 'gtf':
            attrs = 'gene_id'
        else:
            attrs = 'ID'

    if merge_overlapping or as_bed6:

        if isinstance(attrs, six.string_types):
            attrs = [attrs]

        def to_bed(f):
            """
            Given a pybedtools.Interval, return a new Interval with the name
            set according to the kwargs provided above.
            """
            name = attrs_sep.join([f.attrs[i] for i in attrs])
            return pybedtools.Interval(
                f.chrom,
                f.start,
                f.stop,
                name,
                str(f.score),
                f.strand)

        x = x.each(to_bed).saveas()

    if merge_overlapping:
        if bedtools_227_or_later:
            x = x.merge(**_merge_kwargs)
        else:
            def fix_merge(f):
                f = featurefuncs.extend_fields(f, 6)
                return pybedtools.Interval(
                    f.chrom,
                    f.start,
                    f.stop,
                    f[4],
                    '.',
                    f[3])
            x = x.merge(**_merge_kwargs).saveas().each(fix_merge).saveas()

    return x
