import makeit.synthetic.evaluation.rexgen_release.CandRanker.mol_graph as mg
import unittest
import os
import sys
is_py2 = sys.version[0] == '2'
if is_py2:
    import cPickle as pickle
else:
    import pickle as pickle

class TestCRMolGraph(unittest.TestCase):
    def test_01_smiles2graph(self):
        result = mg.smiles2graph("[OH:1][CH3:2]", "[O:1]=[CH2:2]", [(0,1)])
        with open(os.path.join(os.path.dirname(__file__), 'expected/CR_smiles2graph.pkl'), 'rb') as t:
            expected = pickle.load(t)
        for e, r in zip(expected[0], result[0]):
            self.assertTrue((e == r).all())
        self.assertEqual(sorted(expected[1]), sorted(result[1]))


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
