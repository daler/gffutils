import warnings
import expected
from gffutils import example_filename, create, parser, feature
import gffutils
import gffutils.helpers as helpers
import gffutils.gffwriter as gffwriter
import sys
import os
import shutil
import sqlite3
import nose.tools as nt
import difflib
import pprint
import copy
import tempfile
from textwrap import dedent
from nose.tools import assert_raises


testdbfn_gtf = ':memory:'
testdbfn_gff = ':memory:'

def test_update():
    # check both in-memory and file-based dbs
    db = create.create_db(
        example_filename('FBgn0031208.gff'), ':memory:', verbose=False,
        force=True)

    orig_num_features = len(list(db.all_features()))

    f = feature.feature_from_line(
        'chr2L . testing 1 10 . + . ID=testing_feature;n=1',
        dialect=db.dialect, strict=False)

    # no merge strategy required because we're adding a new feature
    db.update([f])
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    assert str(x[0]) == "chr2L	.	testing	1	10	.	+	.	ID=testing_feature;n=1", str(x[0])

    # ought to be one more now . . .
    num_features = len(list(db.all_features()))
    assert num_features == orig_num_features + 1, num_features

    # Now try updating with the same feature, but using merge_strategy="merge",
    # which appends items to attributes ( n=1 --> n=1,2 )
    f = feature.feature_from_line(
        'chr2L . testing 1 10 . + . ID=testing_feature;n=1',
        dialect=db.dialect, strict=False)
    f.attributes['n'] = ['2']
    db.update([f], merge_strategy='merge')
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    assert str(x[0]) == "chr2L	.	testing	1	10	.	+	.	ID=testing_feature;n=1,2"

    # still should have the same number of features as before (still 2)
    num_features = len(list(db.all_features()))
    assert num_features == orig_num_features + 1, num_features


    # Merging while iterating.  e.g., if you're updating children with gene
    # IDs.
    db = create.create_db(example_filename('FBgn0031208.gff'), ':memory:',
                          verbose=False, force=True)
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
            db.update([child], merge_strategy='replace')

    print "\n\nafter\n\n"
    for child in db.children(gene):
        print child.id
        assert child.attributes['gene_id'] == ['FBgn0031208'], (child, child.attributes)

    num_entries = 0
    for gene_recs in list(db.iter_by_parent_childs()):
        # Add attribute to each gene record
        rec = gene_recs[0]
        rec.attributes["new"] = ["new_value"]
        db.update([rec])
        num_entries += 1
    print list(db.all_features())


    assert (num_entries > 1), "Only %d left after update" % (num_entries)


    # Replace
    f = feature.feature_from_line(
        'chr2L . testing 1 10 . + . ID=testing_feature;n=1',
        dialect=db.dialect, strict=False)
    f.attributes['n'] = ['3']
    db.update([f], merge_strategy='replace')
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    assert str(x[0]) == "chr2L	.	testing	1	10	.	+	.	ID=testing_feature;n=3"
    # still should have the same number of features as before (still 2)
    num_features = len(list(db.all_features()))
    assert num_features == orig_num_features + 1, num_features


    # Same thing, but GTF instead of GFF.
    db = create.create_db(
        example_filename('FBgn0031208.gtf'), ':memory:', verbose=False,
        force=True)
    f = feature.feature_from_line('chr2L . testing 1 10 . + . gene_id "fake"; n "1"', strict=False)
    db.update([f], merge_strategy='merge')
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    x = x[0]

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
            verbose=False
        )
        self.c = self.db.conn.cursor()
        self.dialect = self.db.dialect

    def table_test(self):
        expected_tables = ['features', 'relations', 'meta', 'directives', 'autoincrements', 'duplicates']
        self.c.execute('select name from sqlite_master where type="table"')
        observed_tables = [i[0] for i in self.c.execute('select name from sqlite_master where type="table"')]
        assert set(expected_tables) == set(observed_tables), observed_tables

    def _count1(self,featuretype):
        """Count using SQL"""
        self.c.execute('select count() from features where featuretype = ?',(featuretype,))
        results = self.c.fetchone()[0]
        print 'count1("%s") says: %s' % (featuretype,results)
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
        print 'count2("%s") says: %s' % (featuretype, cnt)
        return cnt

    def _count3(self,featuretype):
        """Count with the count_features_of_type method"""
        results = self.db.count_features_of_type(featuretype)
        print 'count3("%s") says: %s' % (featuretype, results)
        return results

    def _count4(self,featuretype):
        """Count by iterating over all features of this type"""
        cnt = 0
        for i in self.db.features_of_type(featuretype):
            cnt += 1
        print 'count4("%s") says: %s' % (featuretype,cnt)
        return cnt

    def featurecount_test(self):
        #  Right number of each featuretype, using multiple different ways of
        #  counting?
        print 'format:', self.dialect['fmt']
        expected_feature_counts = expected.expected_feature_counts[self.dialect['fmt']]
        for featuretype, expected_count in expected_feature_counts.items():
            rawsql_cnt = self._count1(featuretype)
            fileparsed_cnt = self._count2(featuretype)
            count_feature_of_type_cnt = self._count3(featuretype)
            iterator_cnt = self._count4(featuretype)
            print "expected count:", expected_count
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
            print 'observed parents for %s:' % child, set(observed_parents)
            print 'expected parents for %s:' % child, set(expected_parents)
            assert set(observed_parents) == set(expected_parents)


    def test_parents_level_2(self):
        parents1, parents2 = self._expected_parents()
        for child, expected_parents in parents2.items():
            observed_parents = [i.id for i in self.db.parents(child, level=2)]
            print self.db[child]
            print 'observed parents for %s:' % child, set(observed_parents)
            print 'expected parents for %s:' % child, set(expected_parents)
            assert set(observed_parents) == set(expected_parents)


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
    print "Parsed random chromosome successfully."


def test_gffwriter():
    """
    Test GFFWriter.
    """
    print "Testing GFF writer.."
    fn = gffutils.example_filename("unsanitized.gff")
    # Make a copy of it as temporary named file
    temp_f = tempfile.NamedTemporaryFile(delete=False)
    temp_fname_source = temp_f.name
    shutil.copy(fn, temp_fname_source)
    # Now write file in place
    source_first_line = open(temp_fname_source, "r").readline().strip()
    assert (not source_first_line.startswith("#GFF3")), \
           "unsanitized.gff should not have a gffutils-style header."
    db_in = gffutils.create_db(fn, ":memory:")
    # Fetch first record
    rec = db_in.all_features().next()
    ##
    ## Write GFF file in-place test
    ##
    print "Testing in-place writing"
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
    print "  - Wrote GFF file in-place successfully."
    ##
    ## Write GFF file to new file test
    ##
    print "Testing writing to new file"
    new_file = tempfile.NamedTemporaryFile(delete=False)
    gff_out = gffwriter.GFFWriter(new_file.name)
    gff_out.write_rec(rec)
    gff_out.close()
    new_line = open(new_file.name, "r").readline().strip()
    assert new_line.startswith("#GFF3"), \
           "GFFWriter could not write to a new GFF file."
    print "  - Wrote to new file successfully."



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
#     print "First child is not an mRNA"
#     print gene_childs[0].featuretype
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
    print "Testing creation of DB from iterator"
    db_fname = gffutils.example_filename("gff_example1.gff3")
    db = gffutils.create_db(db_fname, ":memory:")
    def my_iterator():
        for rec in db.all_features():
            yield rec
    new_db = gffutils.create_db(my_iterator(), ":memory:")
    print list(new_db.all_features())
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
    print "Sanitized GFF successfully."


def test_region():
    db_fname = gffutils.example_filename("gff_example1.gff3")
    db = gffutils.create_db(db_fname, ":memory:")
    all_in_region = list(db.region("chr1:4000000-5000000"))
    all_minus = list(db.region("chr1:4000000-5000000:-"))
    all_plus = list(db.region("chr1:4000000-5000000:+"))
    all_unstranded = list(db.region("chr1:4000000-5000000:."))

    out_of_range = list(db.region("nowhere:1-100"))

    assert len(all_in_region) == 12
    assert len(all_minus) == 12
    assert len(all_plus) == 0
    assert len(all_unstranded) == 0
    assert len(out_of_range) == 0


def test_nonascii():
    # smoke test (prev. version returned Unicode)
    #
    db = gffutils.create_db(gffutils.example_filename('nonascii'), ":memory:")
    for i in db.all_features():
        # this works in IPython, or using nosetests --with-doctest...
        try:
            print i

        # ...but fails using plain nosetests or when using regular Python
        # interpreter
        except UnicodeEncodeError:
            print unicode(i)


def test_feature_merge():
    # both "n" attribute and "source" field should be merged, since
    # force_merge_fields=['source'].
    gtfdata = dedent("""
    chr1	a	testing	1	10	.	+	.	gene_id "fake"; n "2";
    chr1	b	testing	1	10	.	+	.	gene_id "fake"; n "1";
    """)
    db = gffutils.create_db(gtfdata, ":memory:", from_string=True, merge_strategy='merge', id_spec='gene_id', force_merge_fields=['source'])
    assert db.dialect['fmt'] == 'gtf'
    assert len(list(db.all_features())) == 1
    x = db['fake']
    assert str(x) == 'chr1	a,b	testing	1	10	.	+	.	gene_id "fake"; n "1,2";', str(x)

    gffdata = dedent("""
    chr1	a	testing	1	10	.	+	.	gene_id="fake"; n="2";
    chr1	b	testing	1	10	.	+	.	gene_id="fake"; n="1";
    """)
    db = gffutils.create_db(gffdata, ":memory:", from_string=True, merge_strategy='merge', id_spec='gene_id', force_merge_fields=['source'])
    assert db.dialect['fmt'] == 'gff3'
    assert len(list(db.all_features())) == 1
    x = db['fake']
    assert str(x) == 'chr1	a,b	testing	1	10	.	+	.	gene_id="fake"; n="1,2";', str(x)


    # But when not using force_merge_fields, there should be separate entries;
    # accessing fake and fake_1 should not give FeatureNotFound errors.
    db = gffutils.create_db(gtfdata, ':memory:', from_string=True, merge_strategy='merge', id_spec='gene_id')
    assert len(list(db.all_features())) == 2
    x = db['fake']
    y = db['fake_1']

    db = gffutils.create_db(gffdata, ':memory:', from_string=True, merge_strategy='merge', id_spec='gene_id')
    assert len(list(db.all_features())) == 2
    x = db['fake']
    y = db['fake_1']


    assert_raises(ValueError, gffutils.create_db, gtfdata, ":memory:", from_string=True, merge_strategy='merge', id_spec='gene_id', force_merge_fields=['start'])

    # test that warnings are raised because of strand and frame
    with warnings.catch_warnings(True) as w:
        gffdata = dedent("""
        chr1	a	testing	1	10	.	+	.	gene_id="fake"; n="2";
        chr1	a	testing	1	10	.	-	1	gene_id="fake"; n="1";
        """)
        db = gffutils.create_db(gffdata, ":memory:", from_string=True, merge_strategy='merge', id_spec='gene_id', force_merge_fields=['strand', 'frame'])
        assert db.dialect['fmt'] == 'gff3'
        assert len(list(db.all_features())) == 1
        x = db['fake']
        assert str(x) == 'chr1	a	testing	1	10	.	+,-	1,.	gene_id="fake"; n="1,2";', str(x)
        assert len(w) == 2

def test_add_relation():
    db = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'), ':memory:')
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



if __name__ == "__main__":
    # this test case fails
    #test_attributes_modify()
    #test_sanitize_gff()
    #test_random_chr()
    test_nonascii()
