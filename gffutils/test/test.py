import warnings
from textwrap import dedent
from . import expected
from gffutils import example_filename, create, parser, feature
import gffutils
import gffutils.helpers as helpers
import gffutils.gffwriter as gffwriter
import gffutils.inspect as inspect
import gffutils.iterators as iterators
import sys
import os
import six
import shutil
import threading
import tempfile
from textwrap import dedent
from nose.tools import assert_raises
from six.moves import SimpleHTTPServer
if sys.version_info.major == 3:
    import socketserver as SocketServer
else:
    import SocketServer

import multiprocessing
import json
import tempfile
import shutil
import glob
import difflib

testdbfn_gtf = ':memory:'
testdbfn_gff = ':memory:'



fn = gffutils.example_filename('FBgn0031208.gtf')
def make_db(i):
    """
    Module-level function that can be pickled across processes for
    multiprocessing testing.
    """
    gffutils.create_db(fn, ':memory:', _keep_tempfiles='.%s' % i)
    return i


def test_update():
    # check both in-memory and file-based dbs
    db = create.create_db(
        example_filename('FBgn0031208.gff'), ':memory:', verbose=False,
        keep_order=True,
        force=True)

    orig_num_features = len(list(db.all_features()))

    f = feature.feature_from_line(
        'chr2L . testing 1 10 . + . ID=testing_feature;n=1',
        dialect=db.dialect, strict=False)

    # no merge strategy required because we're adding a new feature
    db.update([f])
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    x = x[0]
    x.keep_order = True
    assert str(x) == "chr2L	.	testing	1	10	.	+	.	ID=testing_feature;n=1", str(x)

    # ought to be one more now . . .
    num_features = len(list(db.all_features()))
    assert num_features == orig_num_features + 1, num_features

    # Now try updating with the same feature, but using merge_strategy="merge",
    # which appends items to attributes ( n=1 --> n=1,2 )
    f = feature.feature_from_line(
        'chr2L . testing 1 10 . + . ID=testing_feature;n=1',
        dialect=db.dialect, strict=False)
    f.keep_order = True
    f.attributes['n'] = ['2']
    db.update([f], merge_strategy='merge')
    x = list(db.features_of_type('testing'))
    assert len(x) == 1

    # Merging does a list(set()) operation, so the order is not guaranteed.
    # Fix it here for testing...
    x = x[0]
    x.attributes['n'].sort()

    assert str(x) == "chr2L	.	testing	1	10	.	+	.	ID=testing_feature;n=1,2", str(x)

    # still should have the same number of features as before (still 2)
    num_features = len(list(db.all_features()))
    assert num_features == orig_num_features + 1, num_features


    # Merging while iterating.  e.g., if you're updating children with gene
    # IDs.
    db = create.create_db(example_filename('FBgn0031208.gff'), ':memory:',
                          verbose=False, force=True, keep_order=True)
    def gen():
        for gene in db.features_of_type('gene'):
            for child in list(db.children(gene)):
                # important: the FBgn0031208.gff file was designed to have some
                # funky features: there are two exons without ID attributes.  These
                # are assigned to ids "exon_1" and "exon_2".  Upon update, with
                # still no ID, we then have two new features "exon_3" and "exon_4".
                # To prevent this issue, we ensure that the ID attribute exists...
                child.attributes['gene_id'] = [gene.id]
                if 'ID' not in child.attributes:
                    child.attributes['ID'] = [child.id]
                yield child

    db.update(gen(), merge_strategy='replace')

    print("\n\nafter\n\n")
    for child in db.children('FBgn0031208'):
        print(child.id)
        assert child.attributes['gene_id'] == ['FBgn0031208'], (child, child.attributes)

    num_entries = 0
    for gene_recs in list(db.iter_by_parent_childs()):
        # Add attribute to each gene record
        rec = gene_recs[0]
        rec.attributes["new"] = ["new_value"]
        db.update([rec], merge_strategy='replace')
        num_entries += 1
    print(list(db.all_features()))


    assert (num_entries > 1), "Only %d left after update" % (num_entries)


    # Replace
    f = feature.feature_from_line(
        'chr2L . testing 1 10 . + . ID=testing_feature;n=1',
        dialect=db.dialect, strict=False)

    f.keep_order = True

    f.attributes['n'] = ['3']
    db.update([f], merge_strategy='replace')
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    assert str(x[0]) == "chr2L	.	testing	1	10	.	+	.	ID=testing_feature;n=3", str(x[0])
    # still should have the same number of features as before (still 2)
    num_features = len(list(db.all_features()))
    assert num_features == orig_num_features + 1, num_features


    # Same thing, but GTF instead of GFF.
    db = create.create_db(
        example_filename('FBgn0031208.gtf'), ':memory:', verbose=False,
        force=True, keep_order=True)
    f = feature.feature_from_line('chr2L . testing 1 10 . + . gene_id "fake"; n "1"', strict=False)
    f.keep_order = True
    db.update([f], merge_strategy='merge')
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    x = x[0]
    x.keep_order = True

    # note the trailing semicolon.  That's because the db's dialect has
    # ['trailing semicolon'] = True.
    assert str(x) == 'chr2L	.	testing	1	10	.	+	.	gene_id "fake"; n "1";', str(x)



class BaseDB(object):
    """
    Generic test class.  Run different versions by subclassing and overriding orig_fn.
    """
    orig_fn = None
    def setup(self):

        def gff_id_func(f):
            if 'ID' in f.attributes:
                return f.attributes['ID'][0]
            elif 'Name' in f.attributes:
                return f.attributes['Name'][0]
            else:
                return '{0.featuretype}:{0.seqid}:{0.start}-{0.end}:{0.strand}'.format(f)

        def gtf_id_func(f):
            if f.featuretype == 'gene':
                if 'gene_id' in f.attributes:
                    return f.attributes['gene_id'][0]
            elif f.featuretype == 'transcript':
                if 'transcript_id' in f.attributes:
                    return f.attributes['transcript_id'][0]
            else:
                return '{0.featuretype}:{0.seqid}:{0.start}-{0.end}:{0.strand}'.format(f)

        if self.orig_fn.endswith('.gtf'): id_func = gtf_id_func
        if self.orig_fn.endswith('.gff'): id_func = gff_id_func
        self.db = create.create_db(
            self.orig_fn,
            ':memory:',
            id_spec=id_func,
            merge_strategy='create_unique',
            verbose=False,
            keep_order=True
        )
        self.c = self.db.conn.cursor()
        self.dialect = self.db.dialect

    def table_test(self):
        expected_tables = ['features', 'relations', 'meta', 'directives', 'autoincrements', 'duplicates', 'sqlite_stat1']
        self.c.execute('select name from sqlite_master where type="table"')
        observed_tables = [i[0] for i in self.c.execute('select name from sqlite_master where type="table"')]
        assert set(expected_tables) == set(observed_tables), observed_tables

    def _count1(self,featuretype):
        """Count using SQL"""
        self.c.execute('select count() from features where featuretype = ?',(featuretype,))
        results = self.c.fetchone()[0]
        print('count1("%s") says: %s' % (featuretype,results))
        return results

    def _count2(self,featuretype):
        """Count GFF lines"""
        cnt = 0
        for line in open(self.orig_fn):
            if line.startswith('#'):
                continue
            L = line.split()

            if len(L) < 3:
                continue

            if L[2] == featuretype:
                cnt += 1
        print('count2("%s") says: %s' % (featuretype, cnt))
        return cnt

    def _count3(self,featuretype):
        """Count with the count_features_of_type method"""
        results = self.db.count_features_of_type(featuretype)
        print('count3("%s") says: %s' % (featuretype, results))
        return results

    def _count4(self,featuretype):
        """Count by iterating over all features of this type"""
        cnt = 0
        for i in self.db.features_of_type(featuretype):
            cnt += 1
        print('count4("%s") says: %s' % (featuretype,cnt))
        return cnt

    def featurecount_test(self):
        #  Right number of each featuretype, using multiple different ways of
        #  counting?
        print('format:', self.dialect['fmt'])
        expected_feature_counts = expected.expected_feature_counts[self.dialect['fmt']]
        for featuretype, expected_count in expected_feature_counts.items():
            rawsql_cnt = self._count1(featuretype)
            fileparsed_cnt = self._count2(featuretype)
            count_feature_of_type_cnt = self._count3(featuretype)
            iterator_cnt = self._count4(featuretype)
            print("expected count:", expected_count)
            assert rawsql_cnt == count_feature_of_type_cnt == iterator_cnt == fileparsed_cnt == expected_count

    def _expected_parents(self):
        if self.dialect['fmt'] == 'gff3':
            parents1 = expected.GFF_parent_check_level_1
            parents2 = expected.GFF_parent_check_level_2
        if self.dialect['fmt'] == 'gtf':
            parents1 = expected.GTF_parent_check_level_1
            parents2 = expected.GTF_parent_check_level_2
        return parents1, parents2

    def test_parents_level_1(self):
        parents1, parents2 = self._expected_parents()
        for child, expected_parents in parents1.items():
            observed_parents = [i.id for i in self.db.parents(child, level=1)]
            print('observed parents for %s:' % child, set(observed_parents))
            print('expected parents for %s:' % child, set(expected_parents))
            assert set(observed_parents) == set(expected_parents)


    def test_parents_level_2(self):
        parents1, parents2 = self._expected_parents()
        for child, expected_parents in parents2.items():
            observed_parents = [i.id for i in self.db.parents(child, level=2)]
            print(self.db[child])
            print('observed parents for %s:' % child, set(observed_parents))
            print('expected parents for %s:' % child, set(expected_parents))
            assert set(observed_parents) == set(expected_parents)


    def test_bed12(self):
        if self.__class__ == TestGFFClass:
            kwargs = dict(block_featuretype='exon', thick_featuretype='CDS', name_field='ID')
        if self.__class__ == TestGTFClass:
            kwargs = dict(block_featuretype='exon', thick_featuretype='CDS', name_field='transcript_id')
        obs = self.db.bed12('FBtr0300689', **kwargs)
        exp = "chr2L	7528	9484	FBtr0300689	0	+	7679	8610	0,0,0	2	588,1292	0,664"
        assert obs == exp


        obs = self.db.bed12('FBtr0300690', **kwargs)
        exp = "chr2L	7528	9484	FBtr0300690	0	+	7679	9276	0,0,0	3	588,397,817	0,664,1139"
        assert obs == exp


class TestGFFClass(BaseDB):
    orig_fn = example_filename('FBgn0031208.gff')

class TestGTFClass(BaseDB):
    orig_fn = example_filename('FBgn0031208.gtf')


def test_random_chr():
    """
    Test on GFF files with random chromosome events.
    """
    gff_fname = gffutils.example_filename("random-chr.gff")
    db = helpers.get_gff_db(gff_fname)
    # Test that we can get children of only a selected type
    gene_id = \
        "chr1_random:165882:165969:-@chr1_random:137473:137600:-@chr1_random:97006:97527:-"
    mRNAs = db.children(gene_id, featuretype="mRNA")
    for mRNA_entry in mRNAs:
        assert (mRNA_entry.featuretype == "mRNA"), \
               "Not all entries are of type mRNA! %s" \
               %(",".join([entry.featuretype for entry in mRNAs]))
    print("Parsed random chromosome successfully.")


def test_gffwriter():
    """
    Test GFFWriter.
    """
    print("Testing GFF writer..")
    fn = gffutils.example_filename("unsanitized.gff")
    # Make a copy of it as temporary named file
    temp_f = tempfile.NamedTemporaryFile(delete=False)
    temp_fname_source = temp_f.name
    shutil.copy(fn, temp_fname_source)
    # Now write file in place
    source_first_line = open(temp_fname_source, "r").readline().strip()
    assert (not source_first_line.startswith("#GFF3")), \
           "unsanitized.gff should not have a gffutils-style header."
    db_in = gffutils.create_db(fn, ":memory:", keep_order=True)
    # Fetch first record
    rec = six.next(db_in.all_features())
    ##
    ## Write GFF file in-place test
    ##
    print("Testing in-place writing")
    gff_out = gffwriter.GFFWriter(temp_fname_source,
                                  in_place=True,
                                  with_header=True)
    gff_out.write_rec(rec)
    gff_out.close()
    # Ensure that the file was written with header
    rewritten = open(temp_fname_source, "r")
    new_header = rewritten.readline().strip()
    assert new_header.startswith("#GFF3"), \
           "GFFWriter serialized files should have a #GFF3 header."
    print("  - Wrote GFF file in-place successfully.")
    ##
    ## Write GFF file to new file test
    ##
    print("Testing writing to new file")
    new_file = tempfile.NamedTemporaryFile(delete=False)
    gff_out = gffwriter.GFFWriter(new_file.name)
    gff_out.write_rec(rec)
    gff_out.close()
    new_line = open(new_file.name, "r").readline().strip()
    assert new_line.startswith("#GFF3"), \
           "GFFWriter could not write to a new GFF file."
    print("  - Wrote to new file successfully.")



# def test_attributes_modify():
#     """
#     Test that attributes can be modified in a GFF record.

#     TODO: This test case fails?
#     """
#     # Test that attributes can be modified
#     db = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'), testdbfn_gff,
#                             verbose=False,
#                             force=True)
#     gene_id = "FBgn0031208"
#     gene_childs = list(db.children(gene_id))
#     print("First child is not an mRNA")
#     print(gene_childs[0].featuretype)
#     assert str(gene_childs[0].attributes) == 'ID=FBtr0300689;Name=CG11023-RB;Parent=FBgn0031208;Dbxref=FlyBase_Annotation_IDs:CG11023-RB;score_text=Strongly Supported;score=11'
#     gene_childs[0].attributes["ID"] = "Modified"
#     assert str(gene_childs[0].attributes) == 'ID=Modified;Name=CG11023-RB;Parent=FBgn0031208;Dbxref=FlyBase_Annotation_IDs:CG11023-RB;score_text=Strongly Supported;score=11;ID=Modified'
#     ###
#     ### NOTE: Would be ideal if database checked that this
#     ### change leaves "dangling" children; i.e. children
#     ### GFF nodes that point to Parent that does not exist.
#     ###


def test_create_db_from_iter():
    """
    Test creation of FeatureDB from iterator.
    """
    print("Testing creation of DB from iterator")
    db_fname = gffutils.example_filename("gff_example1.gff3")
    db = gffutils.create_db(db_fname, ":memory:", keep_order=True)
    def my_iterator():
        for rec in db.all_features():
            yield rec
    new_db = gffutils.create_db(my_iterator(), ":memory:", keep_order=True)
    print(list(new_db.all_features()))
    gene_feats = new_db.all_features(featuretype="gene")
    assert (len(list(gene_feats)) != 0), "Could not load genes from GFF."


def test_sanitize_gff():
    """
    Test sanitization of GFF. Should be merged with GFF cleaning
    I believe unless they are intended to have different functionalities.
    """
    # Get unsanitized GFF
    fn = gffutils.example_filename("unsanitized.gff")
    # Get its database
    db = helpers.get_gff_db(fn)
    # Sanitize the GFF
    sanitized_recs = helpers.sanitize_gff_db(db)
    # Ensure that sanitization work, meaning all
    # starts must be less than or equal to stops
    for rec in sanitized_recs.all_features():
        assert (rec.start <= rec.stop), "Sanitization failed."
    print("Sanitized GFF successfully.")


def test_region():

    db_fname = gffutils.example_filename("FBgn0031208.gff")
    db = gffutils.create_db(db_fname, ":memory:", keep_order=True)

    def _check(item):
        kwargs, expected = item
        try:
            obs = list(db.region(**kwargs))
            assert len(obs) == expected, \
                'expected %s got %s' % (expected, len(obs))
        except expected:
            pass

    regions = [

        # previously failed, see issue #45
        (dict(seqid='chr2L', start=1, end=2e9, completely_within=True), 27),

        (dict(region='chr2L', start=0), ValueError),
        (dict(region='chr2L', end=0), ValueError),
        (dict(region='chr2L', seqid=0), ValueError),


        # these coords should catch everything
        (dict(region="chr2L:7529-12500"), 27),

        # stranded versions:
        (dict(region="chr2L:7529-12500", strand='.'), 0),
        (dict(region="chr2L:7529-12500", strand='+'), 21),
        (dict(region="chr2L:7529-12500", strand='-'), 6),

        # different ways of selecting only that last exon in the last gene:
        (dict(seqid='chr2L', start=11500, featuretype='exon'), 1),
        (dict(seqid='chr2L', start=9500, featuretype='exon', strand='+'), 1),

        # alternative method
        (dict(seqid='chr2L', start=7529, end=12500), 27),

        # since default completely_within=False, this catches anything that
        # falls after 7680.  So it only excludes the 5'UTR, which ends at 7679.
        (dict(seqid='chr2L', start=7680), 26),

        # but completely_within=True will exclude the gene and mRNAs, first
        # exon and the 5'UTR
        (dict(seqid='chr2L', start=7680, completely_within=True), 22),

        # similarly, this will *exclude* anything before 7680
        (dict(seqid='chr2L', end=7680), 5),

        # and also similarly, this will only get us the 5'UTR which is the only
        # feature falling completely before 7680
        (dict(seqid='chr2L', end=7680, completely_within=True), 1),

        # and there's only features from chr2L in this file, so this catches
        # everything too
        (dict(region="chr2L"), 27),

        # using seqid should work similarly to `region` with only chromosome
        (dict(seqid='chr2L'), 27),

        # nonexistent
        (dict(region='nowhere'), 0),

    ]

    for item in regions:
        yield _check, item


def test_nonascii():
    # smoke test (prev. version returned Unicode)
    #
    db = gffutils.create_db(gffutils.example_filename('nonascii'), ":memory:",
                            keep_order=True)
    for i in db.all_features():
        # this works in IPython, or using nosetests --with-doctest...
        try:
            print(i)

        # ...but fails using plain nosetests or when using regular Python
        # interpreter
        except UnicodeEncodeError:
            print(six.text_type(i))


def test_feature_merge():
    # both "n" attribute and "source" field should be merged, since
    # force_merge_fields=['source'].
    gtfdata = dedent("""
    chr1	a	testing	1	10	.	+	.	gene_id "fake"; n "2";
    chr1	b	testing	1	10	.	+	.	gene_id "fake"; n "1";
    """)
    db = gffutils.create_db(gtfdata, ":memory:", from_string=True,
                            merge_strategy='merge', id_spec='gene_id',
                            force_merge_fields=['source'], keep_order=True,
                            sort_attribute_values=True)
    assert db.dialect['fmt'] == 'gtf'
    assert len(list(db.all_features())) == 1
    x = db['fake']
    x.keep_order = True
    x.attributes['n'].sort()
    assert str(x) == 'chr1	a,b	testing	1	10	.	+	.	gene_id "fake"; n "1,2";', str(x)

    gffdata = dedent("""
    chr1	a	testing	1	10	.	+	.	gene_id="fake"; n="2";
    chr1	b	testing	1	10	.	+	.	gene_id="fake"; n="1";
    """)
    db = gffutils.create_db(gffdata, ":memory:", from_string=True,
                            merge_strategy='merge', id_spec='gene_id',
                            force_merge_fields=['source'], keep_order=True)
    assert db.dialect['fmt'] == 'gff3'
    assert len(list(db.all_features())) == 1
    x = db['fake']
    x.attributes['n'].sort()
    x.keep_order = True
    assert str(x) == 'chr1	a,b	testing	1	10	.	+	.	gene_id="fake"; n="1,2";', str(x)


    # But when not using force_merge_fields, there should be separate entries;
    # accessing fake and fake_1 should not give FeatureNotFound errors.
    db = gffutils.create_db(gtfdata, ':memory:', from_string=True,
                            merge_strategy='merge', id_spec='gene_id',
                            keep_order=True)
    assert len(list(db.all_features())) == 2
    x = db['fake']
    y = db['fake_1']

    db = gffutils.create_db(gffdata, ':memory:', from_string=True,
                            merge_strategy='merge', id_spec='gene_id',
                            keep_order=True)
    assert len(list(db.all_features())) == 2
    x = db['fake']
    y = db['fake_1']


    assert_raises(ValueError, gffutils.create_db, gtfdata, ":memory:",
                  from_string=True, merge_strategy='merge', id_spec='gene_id',
                  force_merge_fields=['start'], keep_order=True)

    # test that warnings are raised because of strand and frame
    with warnings.catch_warnings(record=True) as w:
        gffdata = dedent("""
        chr1	a	testing	1	10	.	+	.	gene_id="fake"; n="2";
        chr1	a	testing	1	10	.	-	1	gene_id="fake"; n="1";
        """)
        db = gffutils.create_db(gffdata, ":memory:", from_string=True,
                                merge_strategy='merge', id_spec='gene_id',
                                force_merge_fields=['strand', 'frame'],
                                keep_order=True, sort_attribute_values=True)
        assert db.dialect['fmt'] == 'gff3'
        assert len(list(db.all_features())) == 1
        x = db['fake']
        x.keep_order = True
        x.attributes['n'].sort()
        assert str(x) == 'chr1	a	testing	1	10	.	+,-	.,1	gene_id="fake"; n="1,2";', str(x)
        assert len(w) == 2

def test_add_relation():
    db = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'), ':memory:', keep_order=True)
    L = len(list(db.children('FBgn0031208:3')))
    assert L == 0, L


    def func(parent, child):
        child['Parent'] = child['Parent'] + [parent.id]
        child['exon_parent'] = [parent.id]
        return child

    db.add_relation('FBgn0031208:3', 'CDS_FBgn0031208:1_737', 1, child_func=func)
    L = len(list(db.children('FBgn0031208:3')))
    assert L == 1, L

    L = list(db.children('FBgn0031208:3'))
    x = L[0]
    assert 'FBgn0031208:3' in x['Parent']
    assert x['exon_parent'] == ['FBgn0031208:3']


def test_create_db_from_url():
    """
    Test creation of FeatureDB from URL iterator.
    """
    print("Testing creation of DB from URL iterator")
    # initially run SimpleHTTPServer at port 0 and os will take first available
    Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
    httpd = SocketServer.TCPServer(("", 0), Handler)
    port = str(httpd.socket.getsockname()[1])
    print("serving at port", port)

    # Serving test/data folder
    served_folder = gffutils.example_filename('')
    savedir = os.getcwd()
    os.chdir(served_folder)

    print("Starting SimpleHTTPServer in thread")
    server_thread = threading.Thread(target=httpd.serve_forever)
    server_thread.deamon = True
    server_thread.start()
    try:
        url = ''.join(['http://localhost:', port, '/gff_example1.gff3'])
        db = gffutils.create_db(url, ":memory:", keep_order=True)
        def my_iterator():
            for rec in db.all_features():
                yield rec
        new_db = gffutils.create_db(my_iterator(), ":memory:", keep_order=True)
        print(list(new_db.all_features()))
        gene_feats = new_db.all_features(featuretype="gene")
        assert (len(list(gene_feats)) != 0), "Could not load genes from GFF."

        url = ''.join(['http://localhost:', port, '/gff_example1.gff3.gz'])
        db = gffutils.create_db(url, ":memory:", keep_order=True)
        def my_iterator():
            for rec in db.all_features():
                yield rec
        new_db = gffutils.create_db(my_iterator(), ":memory:", keep_order=True)
        print(list(new_db.all_features()))
        gene_feats = new_db.all_features(featuretype="gene")
        assert (len(list(gene_feats)) != 0), "Could not load genes from GFF."



    finally:
        print('Server shutdown.')
        httpd.shutdown()
        server_thread.join()
        os.chdir(savedir)


def test_empty_files():
    fn = tempfile.NamedTemporaryFile(delete=False).name
    a = open(fn, 'w')
    a.close()
    assert_raises(ValueError, gffutils.create_db, fn, fn + '.db')


def test_false_function():
    # smoke test: before commit ce4b7671f, this would raise "TypeError: object
    # of type 'function' has no len()"
    db = gffutils.create_db(
        gffutils.example_filename('FBgn0031208.gff'),
        ':memory:',
        keep_order=True,
        id_spec=lambda x: False,
        merge_strategy='create_unique'
    )


def test_inspect():
    file_results = inspect.inspect(gffutils.example_filename('FBgn0031208.gff'), verbose=False)
    db_results = inspect.inspect(
        gffutils.create_db(
            gffutils.example_filename('FBgn0031208.gff'),
            ':memory:'),
        verbose=False
    )
    expected =  {

        'featuretype': {
            'intron': 3,
            'five_prime_UTR': 1,
            'exon': 6,
            'mRNA': 4,
            'CDS': 5,
            'pcr_product': 1,
            'three_prime_UTR': 2,
            'protein': 2,
            'gene': 3,
        },

        'feature_count': 27,

        'chrom': {
            'chr2L': 27,
        },

        'attribute_keys': {
            u'': 3,
            'Dbxref': 6,
            'Name': 19,
            'Parent': 20,
            ' Parent': 1,
            'score_text': 2,
            'gbunit': 1,
            'derived_computed_cyto': 1,
            'Derives_from': 2,
            'derived_molecular_weight': 2,
            'score': 2,
            'ID': 25,
            'derived_isoelectric_point': 2,
            'Ontology_term': 1,
        }
    }
    assert file_results == db_results == expected


    # file and db work because db is created from

    kwargs = dict(
        look_for=['chrom', 'strand', 'attribute_keys', 'featuretype'],
        verbose=False,
        limit=10,
    )

    file_results = inspect.inspect(
        gffutils.example_filename('FBgn0031208.gff'),
        **kwargs
    )
    iter_results = inspect.inspect(
        iter(iterators._FileIterator(gffutils.example_filename('FBgn0031208.gff'))),
        **kwargs
    )
    db_results = inspect.inspect(
        gffutils.create_db(
            gffutils.example_filename('FBgn0031208.gff'),
            ':memory:'),
        **kwargs
    )

    expected = {
        'attribute_keys': {
            u'Name': 9,
            u'Parent': 9,
            u'score_text': 2,
            u'gbunit': 1,
            u'derived_computed_cyto': 1,
            u'score': 2,
            u'Dbxref': 3,
            u'ID': 8,
            u'Ontology_term': 1,
        },

        'feature_count': 10,

        'chrom': {u'chr2L': 10},

        'strand': {u'+': 10},

        'featuretype': {
            u'five_prime_UTR': 1,
            u'exon': 3,
            u'mRNA': 2,
            u'CDS': 1,
            'intron': 2,
            u'gene': 1}
    }
    assert file_results == db_results == iter_results == expected


def test_delete():
    db_fname = gffutils.example_filename("gff_example1.gff3")

    # incrementally delete all features
    db = gffutils.create_db(db_fname, ':memory:')
    ids = [i.id for i in db.all_features()]
    current = set(ids)
    for _id in ids:
        db.delete(_id)
        expected = current.difference([_id])
        current = set([i.id for i in db.all_features()])
        assert current == expected, (current, expected)
    assert len(current) == 0

    # same thing, but as a list of Feature objects rather than string IDs
    db = gffutils.create_db(db_fname, ':memory:')
    features = list(db.all_features())
    current = set(features)
    for feature in features:
        db.delete(feature)
        expected = current.difference([feature])
        current = set(list(db.all_features()))
        assert current == expected, (current, expected)
    assert len(current) == 0, current

    # same thing, but use a FeatureDB.
    db1 = gffutils.create_db(db_fname, ':memory:')
    db2 = gffutils.create_db(db_fname, ':memory:')
    db1.delete(db2)
    assert len(list(db1.all_features())) == 0


    db = gffutils.create_db(db_fname, ':memory:')
    db.delete('nonexistent')


def test_iterator_update():
    db_fname = gffutils.example_filename("gff_example1.gff3")
    db = gffutils.create_db(db_fname, ':memory:')
    assert len(list(db.all_features())) == 12
    orig_exon_coords = set([(i.start, i.stop) for i in db.features_of_type('exon')])


    # reset all features to have the same coords of start=1, stop=100
    def gen():
        for i in db.features_of_type('gene'):
            i.start = 1
            i.stop = 100
            yield i

    db.update(gen(), merge_strategy='replace')
    assert len(list(db.all_features())) == 12
    assert len(list(db.features_of_type('gene'))) == 1
    g = six.next(db.features_of_type('gene'))
    assert g.start == 1, g.start
    assert g.stop == 100, g.stop

    # exons should have remained unchanged.
    assert orig_exon_coords == set([(i.start, i.stop) for i in db.features_of_type('exon')])


    def _transform(f):
        f.start = 1
        f.stop = 100
        return f

    db_fname = gffutils.example_filename("gff_example1.gff3")
    db = gffutils.create_db(db_fname, ':memory:')
    db.update(db.features_of_type('gene'), merge_strategy='replace', transform=_transform)
    assert len(list(db.all_features())) == 12
    assert len(list(db.features_of_type('gene'))) == 1
    g = six.next(db.features_of_type('gene'))
    assert g.start == 1, g.start
    assert g.stop == 100, g.stop

    # exons should have remained unchanged.
    assert orig_exon_coords == set([(i.start, i.stop) for i in db.features_of_type('exon')])


def test_tempfiles():

    # specifiy a writeable temp dir for testing
    tempdir = '/tmp/gffutils-test'

    def clean_tempdir():
        tempfile.tempdir = tempdir
        if os.path.exists(tempdir):
            shutil.rmtree(tempdir)
        os.makedirs(tempdir)

    clean_tempdir()

    # default keep_tempfiles=False should give us nothing.
    db = gffutils.create_db(
        gffutils.example_filename('FBgn0031208.gtf'), ':memory:')
    assert len(os.listdir(tempdir)) == 0

    # adding keep_tempfiles=True should give us 1 tempfile for gtf...
    db = gffutils.create_db(
        gffutils.example_filename('FBgn0031208.gtf'), ':memory:', _keep_tempfiles=True)
    filelist = os.listdir(tempdir)
    assert len(filelist) == 1, filelist
    assert filelist[0].endswith('.gffutils')

    #...and another one for gff. This time, make sure the suffix
    db = gffutils.create_db(
        gffutils.example_filename('FBgn0031208.gff'), ':memory:', _keep_tempfiles=True)
    filelist = os.listdir(tempdir)
    assert len(filelist) == 2, filelist
    for i in filelist:
        assert i.endswith('.gffutils')

    # OK, now delete what we have so far...
    clean_tempdir()

    # Make sure that works for custom suffixes
    db = gffutils.create_db(
        gffutils.example_filename('FBgn0031208.gtf'), ':memory:', _keep_tempfiles='.GTFtmp')
    filelist = os.listdir(tempdir)
    assert len(filelist) == 1, filelist
    assert filelist[0].endswith('.GTFtmp')

    clean_tempdir()
    db = gffutils.create_db(
        gffutils.example_filename('FBgn0031208.gtf'), ':memory:', _keep_tempfiles='.GFFtmp')
    filelist = os.listdir(tempdir)
    assert len(filelist) == 1, filelist
    assert filelist[0].endswith('.GFFtmp')

    # Test n parallel instances of gffutils across PROCESSES processes.
    #
    # Note that travis-ci doesn't like it when you use multiple cores, so the
    # .travis.yml file sets this to 1.  This also means that
    #   1) `n` shouldn't be too large because travis-ci will run one at a time,
    #      but more importantly,
    #   2) this will only truly test parallel processes on a local machine with
    #      multiple cpus.
    clean_tempdir()


    # .travis.yml sets the PROCESSES env var; otherwise use all available.
    PROCESSES = int(os.environ.get("PROCESSES", multiprocessing.cpu_count()))
    pool = multiprocessing.Pool(PROCESSES)
    n = 100
    res = pool.map(make_db, range(n))
    assert sorted(list(res)) == list(range(n))
    filelist = os.listdir(tempdir)
    assert len(filelist) == n, len(filelist)

    expected = dedent("""\
        FBtr0300689	chr2L	7529	9484	+	transcript	4681	{"transcript_id":["FBtr0300689"],"gene_id":["FBgn0031208"]}
        FBgn0031208	chr2L	7529	9484	+	gene	4681	{"gene_id":["FBgn0031208"]}
        FBtr0300690	chr2L	7529	9484	+	transcript	4681	{"transcript_id":["FBtr0300690"],"gene_id":["FBgn0031208"]}
        transcript_Fk_gene_1	chr2L	10000	11000	-	transcript	4681	{"transcript_id":["transcript_Fk_gene_1"],"gene_id":["Fk_gene_1"]}
        Fk_gene_1	chr2L	10000	11000	-	gene	4681	{"gene_id":["Fk_gene_1"]}
        transcript_Fk_gene_2	chr2L	11500	12500	-	transcript	4681	{"transcript_id":["transcript_Fk_gene_2"],"gene_id":["Fk_gene_2"]}
        Fk_gene_2	chr2L	11500	12500	-	gene	4681	{"gene_id":["Fk_gene_2"]}
        """)


    def matches_expected(fn):
        """
        Python 3 has unpredictable dictionary ordering. This function checks
        the *semantic* similarity of lines by parsing the attributes into
        a dictonary.
        """
        exp_features = expected.splitlines(True)
        new_features = list(open(fn))
        assert len(exp_features) == len(new_features)
        for expline, newline in zip(exp_features, new_features):
            exp_toks = expline.split()
            new_toks = newline.split()
            assert exp_toks[:-1] == new_toks[:-1]
            assert json.loads(exp_toks[-1]) == json.loads(new_toks[-1])


    # make sure that each of the `n` files matches the expected output.
    for fn in filelist:
        fn = os.path.join(tempdir, fn)
        try:
            matches_expected(fn)
        except AssertionError:
            print(''.join(difflib.ndiff(expected.splitlines(True), this.splitlines(True))))
            raise

    clean_tempdir()

def test_disable_infer():
    """
    tests the new semantics for disabling gene/transcript inference
    """
    # To start, we construct a GTF db by inferring genes and transcripts
    db = gffutils.create_db(gffutils.example_filename('FBgn0031208.gtf'), ':memory:')

    # Then create a file missing transcripts, and another missing genes.
    import tempfile
    tempfile.tempdir = None
    no_transcripts = open(tempfile.NamedTemporaryFile(delete=False).name, 'w')
    no_genes = open(tempfile.NamedTemporaryFile(delete=False).name, 'w')
    for feature in db.all_features():
        if feature.featuretype != 'transcript':
            no_transcripts.write(str(feature) + '\n')
        if feature.featuretype != 'gene':
            no_genes.write(str(feature) + '\n')
    no_genes.close()
    no_transcripts.close()

    no_tx_db = gffutils.create_db(no_transcripts.name, ':memory:', disable_infer_transcripts=True)
    no_gn_db = gffutils.create_db(no_genes.name, ':memory:', disable_infer_genes=True)
    no_xx_db = gffutils.create_db(
        gffutils.example_filename('FBgn0031208.gtf'),
        ':memory:',
        disable_infer_genes=True,
        disable_infer_transcripts=True
    )

    # no transcripts but 3 genes
    assert len(list(no_tx_db.features_of_type('transcript'))) == 0
    assert len(list(no_tx_db.features_of_type('gene'))) == 3

    # no genes but 4 transcripts
    assert len(list(no_gn_db.features_of_type('gene'))) == 0
    assert len(list(no_gn_db.features_of_type('transcript'))) == 4

    # no genes or transcripts
    assert len(list(no_xx_db.features_of_type('gene'))) == 0
    assert len(list(no_xx_db.features_of_type('transcript'))) == 0


def test_deprecation_handler():
    return

    # TODO: when infer_gene_extent actually gets deprecated, test here.
    assert_raises(ValueError, gffutils.create_db,
            gffutils.example_filename('FBgn0031208.gtf'),
            ':memory:',
            infer_gene_extent=False)

def test_nonsense_kwarg():
    assert_raises(TypeError,
                  gffutils.create_db,
                  gffutils.example_filename('FBgn0031208.gtf'),
                  ":memory:",
                  asdf=True)

def test_infer_gene_extent():
    # Before we deprecate this, make sure it still works but emits a warning.
    with warnings.catch_warnings(record=True) as w:
        gffutils.create_db(
            gffutils.example_filename('FBgn0031208.gtf'),
            ':memory:',
            infer_gene_extent=False)
        assert len(w) == 1


# From #79
def test_issue_79():
    gtf = gffutils.example_filename('keep-order-test.gtf')
    db = gffutils.create_db(gtf, 'tmp.db',
                       disable_infer_genes=False,
                       disable_infer_transcripts=False,
                       id_spec={"gene": "gene_id", "transcript": "transcript_id"},
                       merge_strategy="create_unique",
                       keep_order=True,
                            force=True)

    exp = open(gtf).read()
    obs = '\n'.join([str(i) for i in db.all_features()])
    exp_1 = exp.splitlines(True)[0].strip()
    obs_1 = obs.splitlines(True)[0].strip()
    print('EXP')
    print(exp_1)
    print('OBS')
    print(obs_1)
    print('DIFF')
    print(''.join(difflib.ndiff([exp_1], [obs_1])))
    assert obs_1 == exp_1

def test_for_analyze():
    db = gffutils.create_db(
            gffutils.example_filename('FBgn0031208.gtf'),
            'deleteme',
            force=True
    )
    assert db._analyzed()
    db.execute('DROP TABLE sqlite_stat1')
    assert not db._analyzed()

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        db2 = gffutils.FeatureDB('deleteme')
        assert len(w) == 1
        assert "analyze" in str(w[-1].message)
    db.analyze()
    assert db._analyzed()
    os.unlink('deleteme')


def test_issue_82():
    # key-val separator is inside an unquoted attribute value
    x = (
        'Spenn-ch12\tsgn_markers\tmatch\t2621812\t2622049\t.\t+\t.\t'
        'Alias=SGN-M1347;ID=T0028;Note=marker name(s): T0028 SGN-M1347 |identity=99.58|escore=2e-126'
    )
    y = feature.feature_from_line(x)
    assert y.attributes['Note'] == ['marker name(s): T0028 SGN-M1347 |identity=99.58|escore=2e-126']

    gffutils.create_db(gffutils.example_filename('keyval_sep_in_attrs.gff'), ':memory:')

def test_sequence():
    fasta = gffutils.example_filename('dm6-chr2L.fa')
    f = feature.feature_from_line(
        'chr2L	FlyBase	gene	154	170	.	+	.	ID=one;')
    seq = f.sequence(fasta)
    assert seq == 'aCGAGATGATAATATAT'
    assert len(seq) == len(f)
    f.strand = '-'
    seq = f.sequence(fasta)
    assert seq == 'ATATATTATCATCTCGt'
    assert len(seq) == len(f)

def test_issue_85():
    # when start or stop was empty, #85 would fail Should now work with
    # blank fields
    f = feature.feature_from_line('\t'.join([''] * 9))

    # or with "." placeholders
    f = feature.feature_from_line('\t'.join(['.'] * 9))


def test_unquoting():
    # incoming is encoded
    s = (
        'chr1\tAUGUSTUS\tgene\t6950084\t6951407\t0.26\t-\t.\t'
        'ID=INIL01g00009;GeneSymbol=Ndufaf6;Note=NADH dehydrogenase '
        '(ubiquinone) complex I%2C assembly factor 6;GO_Terms=GO:0005743|'
        'GO:0016740|GO:0009058|GO:0032981;PFam=PF00494'
    )
    f = feature.feature_from_line(s, keep_order=True)

    # string representation should be identical
    assert str(f) == s

    # accessing attribute should be decoded
    n = f['Note']
    assert n == ['NADH dehydrogenase (ubiquinone) complex I, assembly factor 6']


def test_unreasonable_unquoting():
    s = (
        'chr1\t.\t.\t1\t2\t0.26\t-\t.\t'
        'newline=%0A;'
        'percent=%25;'
        'null=%00;'
        'comma=%2C;'

        # The first parent is "A," (A with a comma), the second is "B%"
        'Parent=A%2C,B%25,C;'
    )
    f = feature.feature_from_line(s, keep_order=True)
    assert f.attributes['newline'][0] == '\n'
    assert f.attributes['percent'][0] == '%'
    assert f.attributes['null'][0] == '\x00'
    assert f.attributes['comma'][0] == ','

    # Commas indicate
    assert f.attributes['Parent'] == ['A,', 'B%', 'C']
    assert str(f) == s


def test_unquoting_iter():
    s = 'chr1\t.\tgene\t1\t2\t.\t-\t.\tID=%2C;'
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, 'w') as fout:
        fout.write(s + '\n')
    assert list(gffutils.iterators.DataIterator(tmp))[0]['ID'][0] == ','

def test_db_unquoting():
    s = dedent(
        '''
        chr1\t.\tgene\t1\t2\t.\t-\t.\tID=a;Note=%2C;
        chr1\t.\tgene\t1\t2\t.\t-\t.\tID=b;Note=%2C;
        chr1\t.\tgene\t1\t2\t.\t-\t.\tID=c;Note=%2C;
        chr1\t.\tgene\t1\t2\t.\t-\t.\tID=d;Note=%2C;
        chr1\t.\tgene\t1\t2\t.\t-\t.\tID=e;Note=%2C;
        chr1\t.\tgene\t1\t2\t.\t-\t.\tID=f;Note=%2C;
        ''')
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, 'w') as fout:
        fout.write(s + '\n')
    db = gffutils.create_db(tmp, ':memory:', checklines=1)
    assert db['a']['Note'] == [',']
    assert db['b']['Note'] == [',']
    assert db['c']['Note'] == [',']
    assert db['d']['Note'] == [',']
    assert db['e']['Note'] == [',']
    assert db['f']['Note'] == [',']


def test_issue_105():
    fn = gffutils.example_filename('FBgn0031208.gtf')
    home = os.path.expanduser('~')
    newfn = os.path.join(home, '.gffutils.test')
    with open(newfn, 'w') as fout:
        fout.write(open(fn).read())
    f = gffutils.iterators.DataIterator(newfn)
    for i in f:
        pass
    os.unlink(newfn)


def test_issue_107():
    s = dedent(
        '''
        chr1\t.\tgene\t10\t15\t.\t+\t.\tID=b;
        chr1\t.\tgene\t1\t5\t.\t-\t.\tID=a;
        chr2\t.\tgene\t25\t50\t.\t-\t.\tID=c;
        chr2\t.\tgene\t55\t60\t.\t-\t.\tID=d;
        ''')
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, 'w') as fout:
        fout.write(s + '\n')
    db = gffutils.create_db(tmp, ':memory:')
    interfeatures = list(db.interfeatures(
        db.features_of_type('gene', order_by=('seqid', 'start'))))
    assert [str(i) for i in interfeatures] == [
            'chr1\tgffutils_derived\tinter_gene_gene\t6\t9\t.\t.\t.\tID=a,b;',
            'chr2\tgffutils_derived\tinter_gene_gene\t16\t54\t.\t-\t.\tID=c,d;',
    ]


def test_issue_119():
    # First file has these two exons with no ID:
    #
    #   chr2L FlyBase exon  8193  8589  .  +  .  Parent=FBtr0300690
    #   chr2L FlyBase exon  7529  8116  .  +  .  Name=CG11023:1;Parent=FBtr0300689,FBtr0300690
    #
    db0  = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'),':memory:')

    # And this one, a bunch of reads with no IDs anywhere
    db1 = gffutils.create_db(gffutils.example_filename('F3-unique-3.v2.gff'),':memory:')

    # When db1 is updated by db0
    db2 = db1.update(db0)
    assert (
        db2._autoincrements
        == db1._autoincrements
        == {'exon': 2, 'read': 112}
    ), db2._autoincrements


    assert len(list(db0.features_of_type('exon'))) == 6

    # Now we update that with db0 again
    db3 = db2.update(db0, merge_strategy='replace')

    # Using the "replace" strategy, we should have only gotten another 2 exons
    assert len(list(db3.features_of_type('exon'))) == 8

    # Make sure that the autoincrements for exons jumped by 2
    assert (
        db2._autoincrements
        == db3._autoincrements
        == {'exon': 4, 'read': 112}
    ), db2._autoincrements

    # More isolated test, merging two databases each created from the same file
    # which itself contains only a single feature with no ID.
    tmp = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp, 'w') as fout:
        fout.write('chr1\t.\tgene\t10\t15\t.\t+\t.\t\n')

    db4 = gffutils.create_db(tmp, tmp + '.db')
    db5 = gffutils.create_db(tmp, ':memory:')

    assert db4._autoincrements == {'gene': 1}
    assert db5._autoincrements == {'gene': 1}

    db6 = db4.update(db5)

    db7 = gffutils.FeatureDB(db4.dbfn)

    # both db4 and db6 should now have the same, updated autoincrements because
    # they both point to the same db.
    assert db6._autoincrements == db4._autoincrements == {'gene': 2}

    # But db5 was created independently and should have unchanged autoincrements
    assert db5._autoincrements == {'gene': 1}

    # db7 was created from the database pointed to by both db4 and db6. This
    # tests that when a FeatureDB is created it should have the
    # correctly-updated autoincrements read from the db
    assert db7._autoincrements == {'gene': 2}

def test_pr_131():
    db  = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'),':memory:')

    # previously would raise ValueError("No lines parsed -- was an empty
    # file provided?"); now just does nothing
    db2 = db.update([])


def test_pr_133():
    # Previously, merge_attributes would not deep-copy the values from the
    # second dict, and when the values are then modified, the second dict is
    # unintentionally modified.
    d1 = {'a': [1]}
    d2 = {'a': [2]}
    d1a = {'a': [1]}
    d2a = {'a': [2]}
    d3 = gffutils.helpers.merge_attributes(d1, d2)
    assert d1 == d1a, d1
    assert d2 == d2a, d2


def test_pr_139():
    db  = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'),':memory:')
    exons = list(db.features_of_type('exon'))
    inter = list(db.interfeatures(exons))

    # previously, the first exon's attributes would show up in subsequent merged features
    assert exons[0].attributes['Name'][0] not in inter[1].attributes['Name']
    assert exons[0].attributes['Name'][0] not in inter[2].attributes['Name']
    assert exons[0].attributes['Name'][0] not in inter[3].attributes['Name']


def test_pr_144():
    # previously this would fail with:
    #   UnboundLocalError: local variable 'part' referenced before assignment
    f = gffutils.Feature(attributes={'a': ['']})

    # Make sure everything got converted correctly
    assert f.attributes['a'] == ['']
    assert str(f) == ".	.	.	.	.	.	.	.	a"
    g = gffutils.feature.feature_from_line(str(f))
    assert g == f

if __name__ == "__main__":
    # this test case fails
    #test_attributes_modify()
    #test_sanitize_gff()
    #test_random_chr()
    #test_nonascii()
    test_iterator_update()



