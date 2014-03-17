import unittest
from gffutils import helpers


class Test(unittest.TestCase):
    
    def test_merge_attributes(self):
        """
        Tests all possible cases of merging two dictionaries together
        """
        
        x = {'foo': [1], "baz": 1, "buz": [1], "biz": 1, "boo": [1]}
        y = {'bar': [2], "baz": 2, "buz": [2], "biz": 1, "boo": [1]}
        test = helpers.merge_attributes(x, y)
        true = {'foo': [1],
                'bar': [2],
                "baz": [1, 2],
                "boo": [1],
                "buz": [1, 2],
                "biz": [1]}
        self.assertDictEqual(test, true) 
