from copy import deepcopy

from gffutils import merge_criteria as mc
from gffutils.test.synth_test_base import TestWithSynthDB, num_synthetic_features, num_synthetic_overlap


class TestMerge(TestWithSynthDB):
    def _merge_and_compare(self, ids, expected_merge_children, **kwargs):
        features = [self.db[id] for id in ids]
        features_backup = deepcopy(features)
        merged = list(self.db.merge(features, **kwargs))
        dump = self._dump_db()
        m_dump = '\n'.join(map(str, merged))

        # Ensure db wasn't modified
        self.assertEqual(num_synthetic_features, self.db.count_features_of_type(), dump)
        for expected, actual in zip(features_backup, ids):
            self.assertEqual(expected, self.db[actual], '\n'.join(map(str, [expected, self.db[actual]])) + '\n')

        # Ensure input Features not modified
        for expected, actual in zip(features_backup, features):
            self.assertEqual(expected, actual, '\n'.join(map(str, [expected, actual])) + '\n')

        # Ensure output matches expected
        self.assertEqual(len(expected_merge_children), len(merged), m_dump)
        for expected, actual in zip(expected_merge_children, merged):
            c_dump = '\n'.join(map(str, actual.children)) + '\n' + str(actual) + '\n'
            self.assertEqual(len(expected), len(actual.children), '\n'.join(map(str, [expected, actual])) + '\n')
            self.assertTrue(all(child.id in expected for child in actual.children), c_dump)

            if actual.strand != '.':
                self.assertTrue(all(f.strand == actual.strand for f in actual.children), c_dump)

            if actual.frame != '.':
                self.assertTrue(all(f.frame == actual.frame for f in actual.children), c_dump)

            if expected:
                self.assertEqual(actual.start, min(f.start for f in actual.children), c_dump)
                self.assertEqual(actual.end, max(f.end for f in actual.children), c_dump)
                if all(f.seqid == actual.children[0].seqid for f in actual.children):
                    self.assertEqual(actual.seqid, actual.children[0].seqid, c_dump)
                else:
                    self.assertEqual(set(f.seqid for f in actual.children), set(actual.seqid.split(',')), c_dump)

                if all(f.featuretype == actual.children[0].featuretype for f in actual.children):
                    self.assertTrue(all(f.featuretype == actual.featuretype for f in actual.children), c_dump)
                else:
                    self.assertEqual('sequence_feature', actual.featuretype, c_dump)

        return merged

    def test_defaults(self):
        merged = list(self.db.merge(self.db.all_features(order_by=('seqid', 'featuretype', 'strand', 'start'))))
        dump = '\n'.join(map(str, merged))
        self.assertEqual(num_synthetic_features, self.db.count_features_of_type(), dump)
        self.assertEqual(6, len(merged), dump)
        merged_count = 0
        for f in merged:
            if f.children:
                self.assertEqual(num_synthetic_overlap, len(f.children), dump)
                merged_count += 1

        self.assertEqual(1, merged_count)

    def test_none(self):
        merged = list(self.db.merge([]))
        dump = self._dump_db()
        self.assertEqual(num_synthetic_features, self.db.count_features_of_type(), dump)
        self.assertEqual(0, len(merged))

    def test_one(self):
        self._merge_and_compare(['basic1'], [[]])

    def test_no_overlap(self):
        self._merge_and_compare(['basic1', 'no_overlap1'], [[],[]])

    def test_perfect_overlap(self):
        self._merge_and_compare(['basic1', 'perfect_overlap1'], [['basic1', 'perfect_overlap1']])

    def test_overlap(self):
        self._merge_and_compare(['basic1', 'adjacent1'], [['basic1', 'adjacent1']])

    def test_any_overlap(self):
        self._merge_and_compare(['basic1', 'overlap_start1'], [['basic1', 'overlap_start1']], merge_criteria=(mc.seqid, mc.overlap_any_inclusive, mc.strand, mc.feature_type))
        self._merge_and_compare(['overlap_start1', 'basic1'], [['basic1', 'overlap_start1']])
        self._merge_and_compare(['adjacent1', 'basic1'], [['basic1', 'adjacent1']], merge_criteria=(mc.seqid, mc.overlap_any_inclusive, mc.strand, mc.feature_type))
        self._merge_and_compare(['adjacent4', 'basic1'], [['basic1', 'adjacent4']], merge_criteria=(mc.seqid, mc.overlap_any_inclusive, mc.strand, mc.feature_type))

    def test_adjacent_len1(self):
        self._merge_and_compare(['basic1', 'adjacent2'], [['basic1', 'adjacent2']])
        self._merge_and_compare(['adjacent2', 'basic1'], [['basic1', 'adjacent2']], merge_criteria=(mc.seqid, mc.overlap_any_inclusive, mc.strand, mc.feature_type))
        self._merge_and_compare(['basic1', 'adjacent3'], [['basic1', 'adjacent3']], merge_criteria=(mc.seqid, mc.overlap_any_inclusive, mc.strand, mc.feature_type))
        self._merge_and_compare(['adjacent3', 'basic1'], [['basic1', 'adjacent3']])

    def test_end_overlap(self):
        self._merge_and_compare(['basic1', 'overlap_end1'], [['basic1', 'overlap_end1']])

    # TODO test various feature orders

    # TODO test various merge criteria

    # TODO test multiline, this doesn't really work until gffutils better supports multiline features

    # TODO add test for PR #152
