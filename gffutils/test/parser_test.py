from .. import parser
from ..__init__ import example_filename
import attr_test_cases

TEST_FILENAMES = [example_filename(i) for i in [
    'c_elegans_WS199_ann_gff.txt',
    'ensembl_gtf.txt',
    'hybrid1.gff3',
    'ncbi_gff3.txt',
    'c_elegans_WS199_dna_shortened.fa',
    'F3-unique-3.v2.gff',
    'jgi_gff2.txt',
    'wormbase_gff2_alt.txt',
    'c_elegans_WS199_shortened_gff.txt',
    'glimmer_nokeyval.gff3',
    'mouse_extra_comma.gff3',
    'wormbase_gff2.txt']]


def test_split_attrs():
    # nosetests generator for all the test cases in attr_test_cases.  (note no
    # docstring for this test function so that nosetests -v will print the test
    # cases)
    for (attr_str, attr_dict, acceptable_reconstruction) \
            in attr_test_cases.attrs:
        yield attrs_OK, attr_str, attr_dict, acceptable_reconstruction


def attrs_OK(attr_str, attr_dict, acceptable_reconstruction=None):
    """
    Given an attribute string and a dictionary of what you expect, test the
    attribute splitting and reconstruction (invariant roundtrip).

    There are some corner cases for the roundtrip invariance that don't work
    (see attr_test_cases.py for details); `acceptable_reconstruction` handles
    those.
    """
    result, dialect = parser._split_keyvals(attr_str)
    assert result == attr_dict, result

    reconstructed = parser._reconstruct(result, dialect)
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
        p = parser.Parser(filename)
        for i in p:
            continue
