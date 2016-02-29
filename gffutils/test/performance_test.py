import random
from IPython import embed
from nose.plugins.attrib import attr
import unittest
import os

import gffutils

class PerformanceTestFeatureDB(object):
    @classmethod
    def setUpClass(cls):
        # Chromosome sizes for testing fetching features from regions
        cls.chromosome_sizes = {}
        with open(cls.chromsizes_file) as chromosome_sizes:
            for chromosome_size in chromosome_sizes:
                chromosome, size = chromosome_size.split()
                size = int(size)
                cls.chromosome_sizes[chromosome] = size
        random.seed(1842346386)

    def test_make_db(self):
        gffutils.create_db(self.gff_file, ':memory:', merge_strategy="merge")

    def test_iterate_features(self):
        '''
            Walk through all records of a particular featuretype.
        '''
        for featuretype in self.db.featuretypes():
            for feature in self.db.features_of_type(featuretype):
                pass

    def test_find_genes(self):
        '''
            Given a gene id find it's features in db.
        '''
        with open(self.gene_list) as gene_list:
            for gene_id in gene_list:
                gene_id = gene_id.strip()
                gene_from_db = self.db[gene_id]

    def test_find_trainscripts(self):
        '''
            Given a transcript find a gene it belongs to.
        '''
        found_parent = False
        with open(self.transcript_list) as transript_list:
            for transcript in transript_list:
                transcript = transcript.strip()
                found_parent = False
                for gene in self.db.parents(transcript, featuretype=('gene', 'tRNA_gene', 'snoRNA_gene',
                                                                     'snRNA_gene', 'ncRNA_gene', 'rRNA_gene',
                                                                     'pseudogene')):
                    found_parent = True
                assert found_parent, "parent gene not found for transcript %s " % transcript
        assert found_parent, "No single parent gene found"

    def region_fetch_helper(self, number_of_attempts, strand = None, featuretype=None):
        fetched_at_least_one = False
        for fetch in range(number_of_attempts):
            chromosome, size = random.choice(self.chromosome_sizes.items())
            start = random.randint(1, size)
            end = random.randint(start, size)
            for feature in self.db.region(seqid=chromosome, start=start, end=end, strand=strand, featuretype=featuretype):
                fetched_at_least_one = True
        assert fetched_at_least_one, "Didn't fetch any features. Either unlucky or a bug"

    def test_fetch_from_regoins_simple(self):
        '''
            Testing fetching features from random part of random chromosome, simple case.
        '''
        self.region_fetch_helper(number_of_attempts=1000)

    def test_fetch_from_regions_strand(self):
        '''
            Testing fetching features from random part of random chromosome, case with specific strand.
        '''
        self.region_fetch_helper(number_of_attempts=500, strand='+')
        self.region_fetch_helper(number_of_attempts=500, strand='-')

    def test_fetch_from_regions_genes_only(self):
        '''
            Testing fetching features from random part of random chromosome, case with 'gene' as featuretype.
        '''
        self.region_fetch_helper(number_of_attempts=1000, featuretype='gene')

    def test_fetch_from_regions_genes_and_transcripts(self):
        '''
            Testing fetching features from random part of random chromosome, case with ['gene', 'transcript'] as featuretype.
        '''
        self.region_fetch_helper(number_of_attempts=1000, featuretype=['gene', 'transcript'])


@attr('slow')
class TestPerformanceOnSacCer(PerformanceTestFeatureDB, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gff_file = gffutils.example_filename('Saccharomyces_cerevisiae.R64-1-1.83.gff3')
        cls.db = gffutils.create_db(cls.gff_file, 'test.db', merge_strategy="merge")
        cls.chromsizes_file = gffutils.example_filename('Saccharomyces_cerevisiae.R64-1-1.83.chromsizes.txt')
        cls.gene_list = gffutils.example_filename('Saccharomyces_cerevisiae.R64-1-1.83.5000_gene_ids.txt')
        cls.transcript_list = gffutils.example_filename('Saccharomyces_cerevisiae.R64-1-1.83.5000_transcript_ids.txt')
        super(TestPerformanceOnSacCer, cls).setUpClass()
    @classmethod
    def tearDownClass(cls):
        os.remove('test.db')
        super(TestPerformanceOnSacCer, cls).tearDownClass()
