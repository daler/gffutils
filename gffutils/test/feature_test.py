from .. import feature
from .. import helpers


def test_feature_from_line():
    # spaces and tabs should give identical results
    line1 = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    line2 = "chr2L FlyBase exon 7529 8116 . + . Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    assert feature.feature_from_line(line1) == \
            feature.feature_from_line(line2)


def test_default_feature():
    # Default Feature is 8 tab-delimited ".", with a trailing tab
    assert str(feature.Feature()) == \
    ".	.	.	.	.	.	.	.	"


def test_attributes_representations():
    # These different ways of supplying attributes should yield identical
    # results:
    s = ".	.	.	.	.	.	.	.	ID=asdf"
    assert str(feature.Feature(attributes='{"ID": ["asdf"]}')) == s
    assert str(feature.Feature(attributes=dict(ID=["asdf"]))) == s
    assert str(feature.Feature(attributes="ID=asdf")) == s


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
    f = feature.feature_from_line(line)
    assert f.chrom == 'chr2L' == f.seqid
    assert f.end == 8116 == f.stop

    f.chrom = 'fake'
    f.stop = 1
    assert f.chrom == 'fake' == f.seqid
    assert f.stop == 1 == f.end


def test_string_representation():
    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    f = feature.feature_from_line(line)
    assert line == str(f)

    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690	some	more	stuff"
    f = feature.feature_from_line(line)
    assert line == str(f)


def test_pbt_interval_conversion():
    line = "chr2L FlyBase exon 7529 8116 . + . Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    f = feature.feature_from_line(line)
    pbt = helpers.asinterval(f)
    assert pbt.chrom == f.chrom == f.seqid
    assert pbt.start == f.start -1
    assert pbt.stop == f.stop == f.end
    pn = pbt.name
    fn = f.attributes['Name'][0]
    assert pn == fn, '%s, %s' % (pn, fn)


def test_hash():
    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690	some	more	stuff"
    f = feature.feature_from_line(line)
    assert hash(f) == hash(line)


def test_repr():
    line = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690	some	more	stuff"
    f = feature.feature_from_line(line)
    print repr(f)
    print hex(id(f))
    assert repr(f) == ("<Feature exon (chr2L:7529-8116[+]) at %s>" % hex(id(f)))
