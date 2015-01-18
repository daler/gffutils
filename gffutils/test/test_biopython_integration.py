from gffutils import example_filename, create, parser, feature
import gffutils
import gffutils.biopython_integration as bp

def test_roundtrip():
    """
    Feature -> SeqFeature -> Feature should be invariant.
    """
    db_fname = gffutils.example_filename("gff_example1.gff3")
    db = gffutils.create_db(db_fname, ':memory:')
    feature = db['ENSMUSG00000033845']
    feature.keep_order = True
    dialect = feature.dialect
    s = bp.to_seqfeature(feature)
    f = bp.from_seqfeature(s, dialect=dialect, keep_order=True)
    assert feature == f
