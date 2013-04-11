from .. import feature
from .. import helpers


def test_feature_from_line():
    # spaces and tabs should give identical results
    line1 = "chr2L	FlyBase	exon	7529	8116	.	+	.	Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    line2 = "chr2L FlyBase exon 7529 8116 . + . Name=CG11023:1;Parent=FBtr0300689,FBtr0300690"
    assert feature.feature_from_line(line1) == \
            feature.feature_from_line(line2)


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
