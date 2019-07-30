import makeit.synthetic.evaluation.rexgen_release.CoreFinder.mol_graph as mg
import unittest
import numpy as np
import os
import sys
is_py2 = sys.version[0] == '2'
if is_py2:
    import cPickle as pickle
else:
    import pickle as pickle

class TestCFMolGraph(unittest.TestCase):
    def test_01_smiles2graph_batch(self):
        np.set_printoptions(threshold=sys.maxsize)
        result = mg.smiles2graph_batch(["c1cccnc1",'c1nccc2n1ccc2'])
        with open(os.path.join(os.path.dirname(__file__), 'expected/CF_smiles2graph_batch.pkl'), 'rb') as t:
            expected = pickle.load(t)
        for e, r in zip(expected, result):
            self.assertTrue((e == r).all())


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
