import os
import pickle
import unittest

import makeit.synthetic.evaluation.rexgen_release.CandRanker.mol_graph as mg


class TestCRMolGraph(unittest.TestCase):

    @unittest.skip('Non-deterministic')
    def test_01_smiles2graph(self):
        """Test makeit.synthetic.evaluation.rexgen_release.CandRanker.mol_graph.smiles2graph"""
        result = mg.smiles2graph("[OH:1][CH3:2]", "[O:1]=[CH2:2]", [(0, 1)])
        with open(os.path.join(os.path.dirname(__file__), 'expected/CR_smiles2graph.pkl'), 'rb') as t:
            expected = pickle.load(t, encoding='iso-8859-1')
        for e, r in zip(expected[0], result[0]):
            self.assertTrue((e == r).all())
        self.assertEqual(sorted(expected[1]), sorted(result[1]))


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
