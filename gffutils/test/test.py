import expected
from .. import create
from .. import interface
from .. import parser
from .. import feature
from ..__init__ import example_filename
import gffutils
import gffutils.helpers as helpers
import sys
import os
import sqlite3
import nose.tools as nt
import difflib
import pprint
import copy
import tempfile


testdbfn_gtf = ':memory:'
testdbfn_gff = ':memory:'

def test_update():
    # check both in-memory and file-based dbs
    db = create.create_db(
        example_filename('FBgn0031208.gff'), ':memory:', verbose=False,
        force=True)
    f = feature.feature_from_line('chr2L . testing 1 10 . + . ID="testing_feature"', dialect=db.dialect)
    db.update([f])
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    assert str(x[0]) == "chr2L	.	testing	1	10	.	+	.	ID=testing_feature"

    db = create.create_db(
        example_filename('FBgn0031208.gtf'), ':memory:', verbose=False,
        force=True)
    f = feature.feature_from_line('chr2L . testing 1 10 . + . gene_id "fake"', dialect=db.dialect)
    db.update([f], merge_strategy='merge')
    x = list(db.features_of_type('testing'))
    assert len(x) == 1
    assert str(x[0]) == 'chr2L	.	testing	1	10	.	+	.	gene_id "fake";', str(x[0])

class BaseDB(object):
    """
    Generic test class.  Run different versions by subclassing and overriding orig_fn.
    """
    orig_fn = None
    def setup(self):

        def gff_id_func(f):
            if 'ID' in f['attributes']:
                return f['attributes']['ID'][0]
            elif 'Name' in f['attributes']:
                return f['attributes']['Name'][0]
            else:
                return '{0.featuretype}:{0.seqid}:{0.start}-{0.end}:{0.strand}'.format(f)

        def gtf_id_func(f):
            if f['featuretype'] == 'gene':
                if 'gene_id' in f['attributes']:
                    return f['attributes']['gene_id'][0]
            elif f['featuretype'] == 'transcript':
                if 'transcript_id' in f['attributes']:
                    return f['attributes']['transcript_id'][0]
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
        expected_tables = ['features', 'relations', 'meta', 'directives', 'autoincrements']
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
    db_fname = helpers.get_db_fname(gff_fname)
    db = gffutils.FeatureDB(db_fname)
    # Test that we can get children of only a selected type
    gene_id = \
        "chr1_random:165882:165969:-@chr1_random:137473:137600:-@chr1_random:97006:97527:-"
    mRNAs = db.children(gene_id, featuretype="mRNA")
    for mRNA_entry in mRNAs:
        assert (mRNA_entry.featuretype == "mRNA"), \
               "Not all entries are of type mRNA! %s" \
               %(",".join([entry.featuretype for entry in mRNAs]))
    print "Parsed random chromosome successfully."
            

def test_attributes_modify():
    """
    Test that attributes can be modified in a GFF record.

    TODO: This test case fails?
    """
    # Test that attributes can be modified
    db = gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'), testdbfn_gff,
                            verbose=False,
                            force=True)
    gene_id = "FBgn0031208"
    gene_childs = list(db.children(gene_id))
    print "First child is not an mRNA"
    print gene_childs[0].featuretype
    assert str(gene_childs[0].attributes) == 'ID=FBtr0300689;Name=CG11023-RB;Parent=FBgn0031208;Dbxref=FlyBase_Annotation_IDs:CG11023-RB;score_text=Strongly Supported;score=11'
    gene_childs[0].attributes["ID"] = "Modified"
    assert str(gene_childs[0].attributes) == 'ID=Modified;Name=CG11023-RB;Parent=FBgn0031208;Dbxref=FlyBase_Annotation_IDs:CG11023-RB;score_text=Strongly Supported;score=11;ID=Modified'
    ###
    ### NOTE: Would be ideal if database checked that this
    ### change leaves "dangling" children; i.e. children
    ### GFF nodes that point to Parent that does not exist.
    ###
    

    
def test_sanitize_gff():
    """
    Test sanitization of GFF. Should be merged with GFF cleaning
    I believe unless they are intended to have different functionalities.
    """
    # Get unsanitized GFF
    fn = gffutils.example_filename("unsanitized.gff")
    # Get its database
    db_fname = helpers.get_db_fname(fn)
    # Sanitize the GFF
    sanitized_recs = helpers.sanitize_gff(db_fname)
    # Ensure that sanitization work, meaning all
    # starts must be less than or equal to stops
    for rec in sanitized_recs:
        assert (rec.start <= rec.stop), "Sanitization failed."
    # Remove temporary db file
    os.unlink(db_fname)
    print "Sanitized GFF successfully."


if __name__ == "__main__":
    # this test case fails
    #test_attributes_modify()
    test_sanitize_gff()
    test_random_chr()
