from gffutils import Feature, GFFFile
import os
import tempfile

gfffn = os.path.join(
        os.path.abspath(
            os.path.dirname(__file__)),
            'data/hg19.gff')
line = open(gfffn).readline()


def test_feature_str():
    fields = line.split('\t')
    f = Feature(*fields)
    print line.strip()
    print str(f)
    assert line.strip() == str(f)


def test_id():
    fields = line.split('\t')
    f = Feature(*fields)
    print f.id
    assert f.id[0] == 'NR_024540'

def test_len():
    fields = line.split('\t')
    f = Feature(*fields)
    assert len(f) == 15009
    f.start = 1
    f.stop = 10
    assert len(f) == 10


def test_gfffile():
    gff = GFFFile(gfffn)
    f1 = Feature(*line.split('\t'))
    f2 = gff.next()
    assert f1.start == f2.start


def test_gfffile_filetype():
    gff = GFFFile(gfffn)
    assert gff.filetype == 'gff'
