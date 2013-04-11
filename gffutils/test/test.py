import expected
from .. import create
from .. import interface
from .. import parser
from ..__init__ import example_filename
import sys
import os
import sqlite3
import nose.tools as nt
import difflib
import pprint
import copy


testdbfn_gtf = ':memory:'
testdbfn_gff = ':memory:'

def setup():
    create.create_db(example_filename('FBgn0031208.gff'), testdbfn_gff, verbose=False, force=True)
    create.create_db(example_filename('FBgn0031208.gtf'), testdbfn_gtf, verbose=False, force=True)


class BaseDB(object):
    """
    Generic test class.  Run different versions by subclassing and overriding orig_fn.
    """
    orig_fn = None
    def setup(self):
        self.db = create.create_db(self.orig_fn, ':memory:', verbose=False)
        self.c = self.db.conn.cursor()
        self.dialect = self.db.dialect

    def table_test(self):
        """
        Do the right tables exist?
        """
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






class TestGFFClass(BaseDB):
    orig_fn = example_filename('FBgn0031208.gff')

class TestGTFClass(BaseDB):
    orig_fn = example_filename('FBgn0031208.gtf')


