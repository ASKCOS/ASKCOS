import os
import pickle
import unittest

import numpy as np

import makeit.prioritization.templates.relevance as rel


class TestTemplateRelevance(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """This method is run once before all tests in this class."""
        cls.model = rel.RelevanceTemplatePrioritizer()
        cls.model.load_model()

    def test_01_get_topk_from_smi(self):
        """Test that the template relevance model returns the expected result for CCCOCCC"""
        scores, indices = self.model.predict('CCCOCCC', 100, 0.995)

        with open(os.path.join(os.path.dirname(__file__), 'test_data/relevance_01.pkl'), 'rb') as t:
            expected = pickle.load(t)

        self.assertEqual(len(expected[0]), len(scores))
        self.assertEqual(len(expected[1]), len(indices))
        self.assertTrue(np.allclose(expected[0], scores))
        self.assertTrue(np.array_equal(expected[1], indices))

    def test_02_get_topk_from_smi(self):
        """Test that the template relevance model returns the expected result for CCCNc1ccccc1"""
        scores, indices = self.model.predict('CCCNc1ccccc1', 100, 0.995)

        with open(os.path.join(os.path.dirname(__file__), 'test_data/relevance_02.pkl'), 'rb') as t:
            expected = pickle.load(t)

        self.assertEqual(len(expected[0]), len(scores))
        self.assertEqual(len(expected[1]), len(indices))
        self.assertTrue(np.allclose(expected[0], scores))
        self.assertTrue(np.array_equal(expected[1], indices))


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
