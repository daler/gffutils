import gffutils
from textwrap import dedent

TEST_FILENAMES = [
    gffutils.example_filename(i)
    for i in [
        "c_elegans_WS199_ann_gff.txt",
        "ensembl_gtf.txt",
        "hybrid1.gff3",
        "ncbi_gff3.txt",
        "c_elegans_WS199_dna_shortened.fa",
        "F3-unique-3.v2.gff",
        "jgi_gff2.txt",
        "wormbase_gff2_alt.txt",
        "c_elegans_WS199_shortened_gff.txt",
        "glimmer_nokeyval.gff3",
        "mouse_extra_comma.gff3",
        "wormbase_gff2.txt",
    ]
]


def parser_smoke_test():
    """
    Just confirm we can iterate completely through the test files....
    """
    # Don't show the warnings for tests
    import logging

    gffutils.parser.logger.setLevel(logging.CRITICAL)
    for filename in TEST_FILENAMES:
        p = gffutils.iterators._FileIterator(filename)
        for i in p:
            continue


def test_file_iterator():
    db = gffutils.create_db(gffutils.example_filename("hybrid1.gff3"), ":memory:")
    n = len(list(db.all_features()))
    assert n == 6, n


def test_feature_iterator():
    it = gffutils.DataIterator(gffutils.example_filename("hybrid1.gff3"))

    # A DataIterator of a _FileIterator should be a _FeatureIterator
    it = gffutils.DataIterator((i for i in it))
    assert isinstance(it, gffutils.iterators._FeatureIterator), it

    it = gffutils.DataIterator(it)

    db = gffutils.create_db(it, ":memory:")
    n = len(list(db.all_features()))
    assert n == 6, n


def test_update():
    db = gffutils.create_db(
        gffutils.example_filename("FBgn0031208.gtf"),
        ":memory:",
        disable_infer_genes=True,
    )
    f = gffutils.feature.feature_from_line(
        'chr2L	gffutils_derived	gene	11500	12500	.	-	.	gene_id "fake";', strict=False
    )
    db.update([f], disable_infer_genes=True)
    for i in db.all_features():
        print(i.id)
    assert f == db["fake"]


def test_string_iterator():
    gtfdata = dedent(
        """
    chr1	a	testing	1	10	.	+	.	gene_id "fake"; n "2";
    chr1	b	testing	1	10	.	+	.	gene_id "fake"; n "1";
    """
    )
    data = gffutils.iterators.DataIterator(gtfdata, from_string=True)
    n = len(list(data))
    assert n == 2, n

    db = gffutils.create_db(gtfdata, ":memory:", from_string=True)
    n = len(list(db.all_features()))
    assert n == 2, n
