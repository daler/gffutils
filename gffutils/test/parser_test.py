import tempfile
from gffutils import parser, create, feature, iterators, constants, helpers, exceptions
from gffutils import example_filename, create_db
from . import attr_test_cases
from textwrap import dedent

import pytest

TEST_FILENAMES = [
    example_filename(i)
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


def test_directives():
    data = dedent(
        """
    ##directive1 example
    .	.	.	.	.	.	.	.
    .	.	.	.	.	.	.	.
    .	.	.	.	.	.	.	.
    .	.	.	.	.	.	.	.
    """
    )

    it = iterators.DataIterator(data, from_string=True)
    assert it.directives == ["directive1 example"]

    db = create_db(data, dbfn=":memory:", from_string=True, verbose=False)
    assert db.directives == ["directive1 example"], db.directives


@pytest.mark.parametrize("item", attr_test_cases.attrs)
def test_attrs_OK(item):
    """
    Given an attribute string and a dictionary of what you expect, test the
    attribute splitting and reconstruction (invariant roundtrip).

    There are some corner cases for the roundtrip invariance that don't work
    (see attr_test_cases.py for details); `acceptable_reconstruction` handles
    those.
    """
    attr_str, attr_dict, acceptable_reconstruction = item
    result, dialect = parser._split_keyvals(attr_str)
    result = dict(result)
    assert result == attr_dict, result

    reconstructed = parser._reconstruct(result, dialect, keep_order=True)
    if acceptable_reconstruction:
        assert reconstructed == acceptable_reconstruction, reconstructed
    else:
        assert reconstructed == attr_str, reconstructed


def parser_smoke_test():
    """
    Just confirm we can iterate completely through the test files....
    """
    # Don't show the warnings for tests
    import logging

    parser.logger.setLevel(logging.CRITICAL)
    for filename in TEST_FILENAMES:
        p = iterators._FileIterator(filename)
        for i in p:
            continue


def test_empty_recontruct():
    """
    reconstructing attributes with incomplete information returns empty string
    """
    assert parser._reconstruct(None, constants.dialect) == ""
    with pytest.raises(exceptions.AttributeStringError):
        parser._reconstruct(dict(ID="asdf"), None)
    with pytest.raises(exceptions.AttributeStringError):
        parser._reconstruct(None, None)


def test_empty_split_keyvals():
    attrs, dialect = parser._split_keyvals(keyval_str=None)
    assert attrs == feature.dict_class()
    assert dialect == constants.dialect


def test_repeated_keys_conflict():
    """
    if dialect says repeated keys, but len(vals) > 1, then the keys are not
    actually repeated....
    """
    #
    # This is now only checked when infer_dialect is True -- and is disabled
    # when a dialect is provided
    #
    # dialect = constants.dialect.copy()
    # dialect['repeated keys'] = True
    # assert_raises(exceptions.AttributeStringError, parser._split_keyvals, "Parent=1,2,3", dialect)


def test_parser_from_string():
    # DEPRECATED
    #
    # _StringIterator has been removed and is instead handled by DataIterator
    # creating a temp file and returning a _FileIterator.
    return True


def test_valid_line_count():
    p = iterators._FileIterator(example_filename("ncbi_gff3.txt"))
    assert len(list(p)) == 17

    p = iterators._FileIterator(example_filename("hybrid1.gff3"))
    assert len(list(p)) == 6

    p = iterators._FileIterator(example_filename("FBgn0031208.gff"))
    assert len(list(p)) == 27


def test_inconsistent_dialect():
    """
    The second feature does not have a trailing semicolon (wormbase_gff2_alt is
    like this).  But since the first feature does, that's what the db's dialect
    is set to, which can cause errors when parsing attributes.
    """
    db = create.create_db(
        """
    chr1	.	gene	1	100	.	+	.	gene_id "gene1";
    chr1	.	mRNA	1	100	.	+	.	transcript_id "mRNA1"
    """,
        ":memory:",
        from_string=True,
    )
    items = list(db.all_features())
    print(items[0])
    # before, was ['"mRNA1'] -- note extra "
    assert items[1].attributes["transcript_id"] == ["mRNA1"], items[1].attributes[
        "transcript_id"
    ]


def test_attributes():
    s = "chr2L	FlyBase	mRNA	7529	9484	.	+	.	ID=FBtr0300690;Name=CG11023-RC;Parent=FBgn0031208;"
    f = feature.feature_from_line(s)
    f.keep_order = True
    assert str(f) == s, str(f)
