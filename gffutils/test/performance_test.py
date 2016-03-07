'''
Performance testing. Run them with https://github.com/mahmoudimus/nose-timer:

```
nosetests --nocapture -a slow --with-timer
```

WARNING: These tests can take about 1.5 hours to run!

Download required annotation files by running
gffutils/test/data/download-large-annotation-files.sh
'''
import tempfile
import random
import unittest
import os

from nose.plugins import attrib

import gffutils


class PerformanceTestFeatureDB(object):
    '''
    Common scenarios on which to test performance.  Subclass both this class
    and unittest.TestCase to provide it with neeeded files.
    '''
    # All these variables should be given file names in subclasses
    #
    # gff or gtf file to build database from
    gff_file = None

    # each line of chromsizes_file is
    # <chromosome_name><space><chromosome_length>
    chromsizes_file = None

    # each line is gene_id
    gene_list = None

    # each line is transcript_id
    transcript_list = None

    # these args are passed to gffutils.create_db
    # can be overrided in subclasses
    create_db_kwargs = {'disable_infer_transcripts': True,
                        'disable_infer_genes': True,
                        'merge_strategy': 'merge'}

    @classmethod
    def setUpClass(cls):
        '''
        Prepearing for all tests.
        '''
        dbfile_handle, cls.dbfilename = tempfile.mkstemp(
            suffix='.db', prefix='test')
        os.close(dbfile_handle)
        # The annotation files are large so they are not inclided in repo,
        # download it manualy with
        # gffutils/test/data/download-large-annotation-files.sh
        try:
            gffutils.create_db(
                cls.gff_file, cls.dbfilename, **cls.create_db_kwargs)
        except ValueError:
            raise EnvironmentError(
                "Annotation files not found. Download them manualy by "
                "running "
                "gffutils/test/data/download-large-annotation-files.sh")
        cls.db = gffutils.FeatureDB(cls.dbfilename)

        # Chromosome sizes for testing fetching features from regions
        cls.chromosome_sizes = {}
        with open(cls.chromsizes_file) as chromosome_sizes:
            for chromosome_size in chromosome_sizes:
                chromosome, size = chromosome_size.split()
                size = int(size)
                cls.chromosome_sizes[chromosome] = size
        random.seed(1842346386)
        print('Preparation finished')

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.dbfilename)

    def test_make_db(self):
        '''
        Measure time of creating new FeatureDB.
        '''
        gffutils.create_db(self.gff_file, ':memory:', **self.create_db_kwargs)

    def test_iterate_features(self):
        '''
        Walk through all records of a particular featuretype.
        '''
        for featuretype in self.db.featuretypes():
            for dummy_feature in self.db.features_of_type(featuretype):
                pass

    def test_find_genes(self):
        '''
        Given a gene id find its features in db.
        '''
        with open(self.gene_list) as gene_list:
            for gene_id in gene_list:
                gene_id = gene_id.strip()
                dummy_gene_from_db = self.db[gene_id]

    def test_find_trainscripts(self):
        '''
        Given a transcript find the gene it belongs to.
        '''
        found_parent = False
        with open(self.transcript_list) as transript_list:
            for transcript in transript_list:
                transcript = transcript.strip()
                found_parent = False
                for dummy_gene in self.db.parents(
                        transcript, featuretype=('gene', 'tRNA_gene',
                                                 'snoRNA_gene', 'snRNA_gene',
                                                 'ncRNA_gene', 'rRNA_gene',
                                                 'pseudogene')):
                    found_parent = True
                assert found_parent, (
                    "parent gene not found for transcript %s " % transcript)
        assert found_parent, "No single parent gene found"

    def region_fetch_helper(self, number_of_repeats, strand=None,
                            featuretype=None):
        '''
        Helper for fetching features from random region number_of_repeats
        times.
        '''
        fetched_at_least_one = False
        for dummy_repeat in range(number_of_repeats):
            chromosome, size = random.choice(self.chromosome_sizes.items())
            start = random.randint(1, size)
            end = random.randint(start, size)
            for dummy_feature in self.db.region(seqid=chromosome, start=start,
                                                end=end, strand=strand,
                                                featuretype=featuretype):
                fetched_at_least_one = True
        assert fetched_at_least_one, \
            "Didn't fetch any features. Either unlucky or a bug"

    def test_fetch_from_regoins_simple(self):
        '''
        Testing fetching features from random part of random chromosome, simple
        case.
        '''
        self.region_fetch_helper(number_of_repeats=500)

    def test_fetch_from_regions_strand(self):
        '''
        Testing fetching features from random part of random chromosome,
        specific strand.
        '''
        self.region_fetch_helper(number_of_repeats=200, strand='+')
        self.region_fetch_helper(number_of_repeats=200, strand='-')

    def test_fetch_from_regions_genes_only(self):
        '''
        Testing fetching features from random part of random chromosome, 'gene'
        as featuretype.
        '''
        self.region_fetch_helper(number_of_repeats=400, featuretype='gene')

    def test_fetch_from_regions_genes_and_transcripts(self):
        '''
            Testing fetching features from random part of random chromosome,
            ['gene', 'transcript'] as featuretype.
        '''
        self.region_fetch_helper(
            number_of_repeats=400, featuretype=['gene', 'transcript'])


@attrib.attr('slow')
class TestPerformanceOnSacCer(PerformanceTestFeatureDB, unittest.TestCase):
    '''
        Test frequent scenarios on medium size genome of yeast.
    '''
    gff_file = gffutils.example_filename(
        'Saccharomyces_cerevisiae.R64-1-1.83.gff3')
    chromsizes_file = gffutils.example_filename(
        'Saccharomyces_cerevisiae.R64-1-1.83.chromsizes.txt')
    gene_list = gffutils.example_filename(
        'Saccharomyces_cerevisiae.R64-1-1.83.5000_gene_ids.txt')
    transcript_list = gffutils.example_filename(
        'Saccharomyces_cerevisiae.R64-1-1.83.5000_transcript_ids.txt')


@attrib.attr('slow')
class TestPerformanceOnMouse(PerformanceTestFeatureDB, unittest.TestCase):
    '''
        Test frequent scenarios on large genome of mouse.
    '''
    gff_file = gffutils.example_filename('gencode.vM8.annotation.gff3')
    chromsizes_file = gffutils.example_filename('gencode.vM8.chromsizes.txt')
    gene_list = gffutils.example_filename('gencode.vM8.5000_gene_ids.txt')
    transcript_list = gffutils.example_filename(
        'gencode.vM8.5000_transcript_ids.txt')
