import expected
from .. import create
from .. import interface
from .. import parser
from .. import feature
from ..__init__ import example_filename
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
                return '{featuretype}:{seqid}:{start}-{end}:{strand}'.format(**f)

        def gtf_id_func(f):
            if f['featuretype'] == 'gene':
                if 'gene_id' in f['attributes']:
                    return f['attributes']['gene_id'][0]
            elif f['featuretype'] == 'transcript':
                if 'transcript_id' in f['attributes']:
                    return f['attributes']['transcript_id'][0]
            else:
                return '{featuretype}:{seqid}:{start}-{end}:{strand}'.format(**f)

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


