from unittest import TestCase
import gffutils
synthetic_path = gffutils.example_filename('synthetic.gff3')
num_synthetic_features = 18
num_synthetic_overlap = 13


class TestWithSynthDB(TestCase):
    def setUp(self):
        self.db = gffutils.create_db(synthetic_path, ":memory:", merge_strategy='create_unique')  # type: gffutils.FeatureDB
        self.assertEqual(num_synthetic_features, self.db.count_features_of_type())

    def _dump_db(self):
        return '\n'.join(map(str, self.db.all_features()))
