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

def EXPECTED_DATA():
    # list the children and their expected first-order parents for the GFF test file.
    GFF_parent_check_level_1 = {'FBtr0300690':['FBgn0031208'],
                                'FBtr0300689':['FBgn0031208'],
                                'CG11023:1':['FBtr0300689','FBtr0300690'],
                                'five_prime_UTR_FBgn0031208:1_737':['FBtr0300689','FBtr0300690'],
                                'CDS_FBgn0031208:1_737':['FBtr0300689','FBtr0300690'],
                                'intron_FBgn0031208:1_FBgn0031208:2':['FBtr0300690'],
                                'intron_FBgn0031208:1_FBgn0031208:3':['FBtr0300689'],
                                'FBgn0031208:3':['FBtr0300689'],
                                'CDS_FBgn0031208:3_737':['FBtr0300689'],
                                'CDS_FBgn0031208:2_737':['FBtr0300690'],
                                'exon:chr2L:8193-8589:+':['FBtr0300690'],
                                'intron_FBgn0031208:2_FBgn0031208:4':['FBtr0300690'],
                                'three_prime_UTR_FBgn0031208:3_737':['FBtr0300689'],
                                'FBgn0031208:4':['FBtr0300690'],
                                'CDS_FBgn0031208:4_737':['FBtr0300690'],
                                'three_prime_UTR_FBgn0031208:4_737':['FBtr0300690'],
                               }

    # and second-level . . . they should all be grandparents of the same gene.
    GFF_parent_check_level_2 = {
                                'CG11023:1':['FBgn0031208'],
                                'five_prime_UTR_FBgn0031208:1_737':['FBgn0031208'],
                                'CDS_FBgn0031208:1_737':['FBgn0031208'],
                                'intron_FBgn0031208:1_FBgn0031208:2':['FBgn0031208'],
                                'intron_FBgn0031208:1_FBgn0031208:3':['FBgn0031208'],
                                'FBgn0031208:3':['FBgn0031208'],
                                'CDS_FBgn0031208:3_737':['FBgn0031208'],
                                'CDS_FBgn0031208:2_737':['FBgn0031208'],
                                'FBgn0031208:2':['FBgn0031208'],
                                'intron_FBgn0031208:2_FBgn0031208:4':['FBgn0031208'],
                                'three_prime_UTR_FBgn0031208:3_737':['FBgn0031208'],
                                'FBgn0031208:4':['FBgn0031208'],
                                'CDS_FBgn0031208:4_737':['FBgn0031208'],
                                'three_prime_UTR_FBgn0031208:4_737':['FBgn0031208'],
                               }

    # Same thing for GTF test file . . .
    GTF_parent_check_level_1 = {
                                'exon:chr2L:7529-8116:+':['FBtr0300689','FBtr0300690'],
                                'exon:chr2L:8193-9484:+':['FBtr0300689'],
                                'exon:chr2L:8193-8589:+':['FBtr0300690'],
                                'exon:chr2L:8668-9484:+':['FBtr0300690'],
                                'exon:chr2L:10000-11000:-':['transcript_Fk_gene_1'],
                                'exon:chr2L:11500-12500:-':['transcript_Fk_gene_2'],
                                'CDS:chr2L:7680-8116:+':['FBtr0300689','FBtr0300690'],
                                'CDS:chr2L:8193-8610:+':['FBtr0300689'],
                                'CDS:chr2L:8193-8589:+':['FBtr0300690'],
                                'CDS:chr2L:8668-9276:+':['FBtr0300690'],
                                'CDS:chr2L:10000-11000:-':['transcript_Fk_gene_1'],
                                'FBtr0300689':['FBgn0031208'],
                                'FBtr0300690':['FBgn0031208'],
                                'transcript_Fk_gene_1':['Fk_gene_1'],
                                'transcript_Fk_gene_2':['Fk_gene_2'],
                                'start_codon:chr2L:7680-7682:+':['FBtr0300689','FBtr0300690'],
                                'start_codon:chr2L:10000-11002:-':['transcript_Fk_gene_1'],
                                'stop_codon:chr2L:8611-8613:+':['FBtr0300689'],
                                'stop_codon:chr2L:9277-9279:+':['FBtr0300690'],
                                'stop_codon:chr2L:11001-11003:-':['transcript_Fk_gene_1'],
                                }


    GTF_parent_check_level_2 = {
            'exon:chr2L:7529-8116:+':['FBgn0031208'],
                                'exon:chr2L:8193-9484:+':['FBgn0031208'],
                                'exon:chr2L:8193-8589:+':['FBgn0031208'],
                                'exon:chr2L:8668-9484:+':['FBgn0031208'],
                                'exon:chr2L:10000-11000:-':['Fk_gene_1'],
                                'exon:chr2L:11500-12500:-':['Fk_gene_2'],
                                'CDS:chr2L:7680-8116:+':['FBgn0031208'],
                                'CDS:chr2L:8193-8610:+':['FBgn0031208'],
                                'CDS:chr2L:8193-8589:+':['FBgn0031208'],
                                'CDS:chr2L:8668-9276:+':['FBgn0031208'],
                                'CDS:chr2L:10000-11000:-':['Fk_gene_1'],
                                'FBtr0300689':['FBgn0031208'],
                                'FBtr0300690':['FBgn0031208'],
                                'transcript_Fk_gene_1':['Fk_gene_1'],
                                'transcript_Fk_gene_2':['Fk_gene_2'],
                                'start_codon:chr2L:7680-7682:+':['FBgn0031208'],
                                'start_codon:chr2L:10000-11002:-':['Fk_gene_1'],
                                'stop_codon:chr2L:8611-8613:+':['FBgn0031208'],
                                'stop_codon:chr2L:9277-9279:+':['FBgn0031208'],
                                'stop_codon:chr2L:11001-11003:-':['Fk_gene_1'],
                               }

    expected_feature_counts = {
                'GFF':{'gene':3,
                       'mRNA':4,
                       'exon':6,
                       'CDS':5,
                       'five_prime_UTR':1,
                       'intron':3,
                       'pcr_product':1,
                       'protein':2,
                       'three_prime_UTR':2},
                'GTF':{'gene':3,
                       'mRNA':4,
                       'CDS':5,
                       'exon':6,
                       'start_codon':2,
                       'stop_codon':3}
                }

    expected_features = {'GFF':['gene',
                                'mRNA',
                                'protein',
                                'five_prime_UTR',
                                'three_prime_UTR',
                                'pcr_product',
                                'CDS',
                                'exon',
                                'intron'],
                        'GTF':['gene',
                               'mRNA',
                               'CDS',
                               'exon',
                               'start_codon',
                               'stop_codon']}

    return GFF_parent_check_level_1,GFF_parent_check_level_2,GTF_parent_check_level_1,GTF_parent_check_level_2,expected_feature_counts,expected_features

(
GFF_parent_check_level_1,
GFF_parent_check_level_2,
GTF_parent_check_level_1,
GTF_parent_check_level_2,
expected_feature_counts,
expected_features,
) = EXPECTED_DATA()

def test_clean_gff():
    # test the "full" cleaning -- remove some featuretypes, do sanity-checking,
    # add chr
    fn = gffutils.example_filename('dirty.gff')
    gffutils.clean_gff(fn, newfn='cleaned.tmp',featuretypes_to_remove=['pcr_product','protein'],addchr=True)
    observed = open('cleaned.tmp').readlines()
    expected = open(gffutils.example_filename('fully-cleaned.gff')).readlines()
    assert observed==expected
    os.unlink('cleaned.tmp')
    gffutils.clean_gff(fn, featuretypes_to_remove=None, sanity_check=False)
    observed = open(gffutils.example_filename('dirty.gff.cleaned')).read()
    expected = open(gffutils.example_filename('basic-cleaned.gff')).read()
    assert observed == expected
    os.unlink(gffutils.example_filename('dirty.gff.cleaned'))


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


def test_inspect_featuretypes():
    observed = gffutils.inspect_featuretypes(gffutils.example_filename('FBgn0031208.gff'))
    observed.sort()
    expected = ['CDS', 'exon', 'five_prime_UTR', 'gene', 'intron', 'mRNA', 'pcr_product', 'protein', 'three_prime_UTR']
    print observed
    print expected
    assert observed == expected

class GenericDBClass(object):
    featureclass = None

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
    gene_id = "chr1_random:165882:165969:-@chr1_random:137473:137600:-@chr1_random:97006:97527:-"
    mRNAs = db.children(gene_id, featuretype="mRNA")
    for mRNA_entry in mRNAs:
        assert (mRNA_entry.featuretype == "mRNA"), \
               "Not all entries are of type mRNA! %s" \
               %(",".join([entry.featuretype for entry in mRNAs]))
    print "Parsed random chromosome successfully."
            

def __test_attributes_modify():
    """
    Test that attributes can be modified in a GFF record.

    TODO: This test case fails?
    """
    # Test that attributes can be modified
    gffutils.create_db(gffutils.example_filename('FBgn0031208.gff'), testdbfn_gff,
                       verbose=False,
                       force=True)
    db = gffutils.FeatureDB(testdbfn_gff)
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
    

if __name__ == "__main__":
    # this test case fails
    #test_attributes_modify()
    test_sanitize_gff()
    test_random_chr()
