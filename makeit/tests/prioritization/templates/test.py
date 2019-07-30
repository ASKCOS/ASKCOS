import makeit.prioritization.templates.relevance as rel
import unittest
import os
import sys
is_py2 = sys.version[0] == '2'
if is_py2:
    import cPickle as pickle
else:
    import pickle as pickle

class TestTemplateRelevance(unittest.TestCase):
    def setUp(self):
        self.model = rel.RelevanceTemplatePrioritizer(use_tf=True)
        self.model.load_model()

    def test_01_get_topk_from_smi(self):
        result = self.model.get_topk_from_smi('CCCOCCC')
        with open(os.path.join(os.path.dirname(__file__), 'expected/relevance_01.pkl'), 'rb') as t:
            expected = pickle.load(t)
        self.assertEqual(expected, result)

    def test_02_get_topk_from_smi(self):
        result = self.model.get_topk_from_smi('CCCNc1ccccc1')
        with open(os.path.join(os.path.dirname(__file__), 'expected/relevance_02.pkl'), 'rb') as t:
            expected = pickle.load(t)
        self.assertEqual(expected, result)


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
