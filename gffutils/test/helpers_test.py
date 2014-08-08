import unittest
from gffutils import helpers
from gffutils import feature


class Test(unittest.TestCase):
    
    def test_merge_attributes(self):
        """
        Tests all possible cases of merging two dictionaries together
        """
        x = {'foo': [1], "baz": 1, "buz": [1], "biz": 1, "boo": [1]}
        y = {'bar': [2], "baz": 2, "buz": [2], "biz": 1, "boo": [1]}
        test = helpers.merge_attributes(x, y)

        for k, v in list(test.items()):
            test[k] = sorted(v)

        true = {'foo': [1],
                'bar': [2],
                "baz": [1, 2],
                "boo": [1],
                "buz": [1, 2],
                "biz": [1]}
        self.assertDictEqual(test, true) 

    def test_merge_Attributes(self):
        f1 = feature.feature_from_line('chr2L . testing 1 10 . + . foo=1; baz=1; buz=1; biz=1; boo=1;', strict=False)
        f2 = feature.feature_from_line('chr2L . testing 1 10 . + . bar=2; baz=2; buz=2; biz=1; boo=1;', strict=False)
        test = helpers.merge_attributes(f1.attributes, f2.attributes)

        for k, v in list(test.items()):
            test[k] = sorted(v)

        true = {'foo': ['1'],
                'bar': ['2'],
                "baz": ['1', '2'],
                "boo": ['1'],
                "buz": ['1', '2'],
                "biz": ['1']}
        self.assertDictEqual(test, true)

