from gffutils import parser, feature, helpers, constants

def test_feature_from_line():
    # spaces and tabs should give identical results
    line1 = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    line2 = "chr2L FlyBase exon 7529 8116 . + . Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    assert feature.feature_from_line(line1, strict=False, keep_order=True) == \
            feature.feature_from_line(line2, strict=False, keep_order=True)


def test_default_feature():
    # Default Feature is 8 tab-delimited ".", with a trailing tab
    assert str(feature.Feature()) == \
    ".	.	.	.	.	.	.	.	"


def test_attributes_representations():
    # These different ways of supplying attributes should yield identical
    # results:
    s = ".	.	.	.	.	.	.	.	ID=asdf"
    for item in (
        '{"ID": ["asdf"]}',
        dict(ID=["asdf"]),
        "ID=asdf"
    ):
        result = str(feature.Feature(attributes=item))
        assert result == s, result


def test_default_start_stop():
    # Whether start or end is "." or None, attribute should always be None and
    # printing should show "."
    c = ['.', None]
    for i1 in c:
        for i2 in c:
            f = feature.Feature(start=i1, end=i2)
            assert f.start is None
            assert f.end is None
            assert f.stop is None
            assert str(f) == ".	.	.	.	.	.	.	.	", str(f)

    # Make sure zero works (protects against sloppy "if start:")
    f = feature.Feature(start=0, end=0)
    assert f.start == f.end == f.stop == 0
    assert str(f) == ".	.	.	0	0	.	.	.	", str(f)


def test_aliases():
    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    f = feature.feature_from_line(line, keep_order=True)
    assert f.chrom == 'chr2L' == f.seqid
    assert f.end == 8116 == f.stop

    f.chrom = 'fake'
    f.stop = 1
    assert f.chrom == 'fake' == f.seqid
    assert f.stop == 1 == f.end


def test_string_representation():
    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    f = feature.feature_from_line(line, keep_order=True)
    assert line == str(f), str(f)

    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690	some	more	stuff"
    f = feature.feature_from_line(line, keep_order=True)
    assert line == str(f)


def test_pbt_interval_conversion():
    try:
        import pybedtools
    except ImportError:
        return
    line = "chr2L FlyBase exon 7529 8116 . + . Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    f = feature.feature_from_line(line, strict=False, keep_order=True)
    pbt = helpers.asinterval(f)
    assert pbt.chrom == f.chrom == f.seqid
    assert pbt.start == f.start -1
    assert pbt.stop == f.stop == f.end
    pn = pbt.name
    fn = f.attributes['Name'][0]
    assert pn == fn, '%s, %s' % (pn, fn)


def test_hash():
    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690	some	more	stuff"
    f = feature.feature_from_line(line, keep_order=True)
    assert hash(f) == hash(line)


def test_repr():
    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690	some	more	stuff"
    f = feature.feature_from_line(line, keep_order=True)
    print(repr(f))
    print(hex(id(f)))
    assert repr(f) == ("<Feature exon (chr2L:7529-8116[+]) at %s>" % hex(id(f)))

def test_attribute_order():

    # default order is gene_id, transcript_id.  But feature_from_line -- if
    # dialect not provided -- will infer its own dialect.  In this case,
    # transcript_id comes first.
    attributes = 'transcript_id "mRNA1"; gene_id "gene1";'
    a = feature.feature_from_line(
        """
        chr1	.	mRNA	1	100	.	+	.	%s
        """ % attributes, strict=False, keep_order=True)
    a.strict = True
    a.keep_order = True
    assert str(a) == 'chr1	.	mRNA	1	100	.	+	.	transcript_id "mRNA1"; gene_id "gene1";', str(a)

    # ensure that using the default dialect uses the default order (and
    # indidentally converts to GFF3 format)
    orig_dialect = a.dialect
    a.dialect = constants.dialect
    a.keep_order = True
    assert str(a) == 'chr1	.	mRNA	1	100	.	+	.	gene_id=gene1;transcript_id=mRNA1', str(a)

    # adding an attribute shoud always result in that attribute coming last (as
    # long as that attribute is not in the dialect order)
    a['dummy'] = ['asdf']
    a.strict = True
    assert str(a) == 'chr1	.	mRNA	1	100	.	+	.	gene_id=gene1;transcript_id=mRNA1;dummy=asdf', str(a)


def test_unjsonify():
    attributes, dialect = parser._split_keyvals('transcript_id "mRNA1"')
    assert attributes == {'transcript_id': ['mRNA1']}, attributes

    s = helpers._jsonify(attributes)
    assert s == '{"transcript_id":["mRNA1"]}', s

    d = helpers._unjsonify(s, isattributes=True)
    assert d == attributes

class IsolatedTestCase(object):
    """
    Isolated test case for checking that the module-level
    constants.always_return_list works.

    This was needed because having this test as a function caused other tests
    to fail even though constants.always_return_list was put back to its
    original setting.  Apparently nose runs tests concurrently in the same
    namespace or something?  Anyway, these setup/teardowns do the trick.
    """
    def setup(self):
        constants.always_return_list = False

    def teardown(self):
        constants.always_return_list = True

    def test_feature_single_item(self):
        line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690	some	more	stuff"
        f = feature.feature_from_line(line, keep_order=True)
        assert f['Name'] == ['CG11023:1']
