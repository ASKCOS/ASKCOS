import os
import pickle
import sys
import unittest

import numpy as np

import makeit.synthetic.evaluation.rexgen_release.CoreFinder.mol_graph as mg


class TestCFMolGraph(unittest.TestCase):

    def test_01_smiles2graph_batch(self):
        """Test makeit.synthetic.evaluation.rexgen_release.CoreFinder.mol_graph.smiles2graph_batch"""
        np.set_printoptions(threshold=sys.maxsize)
        result = mg.smiles2graph_batch(["c1cccnc1", 'c1nccc2n1ccc2'])
        with open(os.path.join(os.path.dirname(__file__), 'test_data/CF_smiles2graph_batch.pkl'), 'rb') as t:
            expected = pickle.load(t, encoding='iso-8859-1')
        for e, r in zip(expected, result):
            self.assertTrue((e == r).all())


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
