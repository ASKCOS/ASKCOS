import makeit.prioritization.precursors.scscore as sc
import makeit.prioritization.precursors.mincost as mc
import unittest

class TestSCScore(unittest.TestCase):
    def setUp(self):
        self.model = sc.SCScorePrecursorPrioritizer()

    def test_01_get_score_from_smiles_1024bool(self):
        self.model.load_model(model_tag='1024bool')
        result = self.model.get_score_from_smiles('CCCOCCC', noprice=True)
        expected = 1.43226093752
        self.assertLess(abs(expected - result), 1e-4)

    def test_02_get_score_from_smiles_1024bool(self):
        self.model.load_model(model_tag='1024bool')
        result = self.model.get_score_from_smiles('CCCC', noprice=True)
        expected = 1.32272217424
        self.assertLess(abs(expected - result), 1e-4)

    def test_03_get_priority_2048bool(self):
        self.model.load_model(model_tag='2048bool', FP_len=2048)
        result = self.model.get_priority('CCCOCCC')
        expected = -0.02
        self.assertLess(abs(expected - result), 1e-4)

    def test_04_get_priority_2048bool(self):
        self.model.load_model(model_tag='2048bool', FP_len=2048)
        result = self.model.get_priority('CCCNc1ccccc1')
        expected = -0.45
        self.assertLess(abs(expected - result), 1e-4)

    def test_05_get_priority_1024uint8(self):
        self.model.load_model(model_tag='1024uint8')
        result = self.model.get_priority('CCCOCCC')
        expected = -0.02
        self.assertLess(abs(expected - result), 1e-4)

    def test_06_get_priority_1024uint8(self):
        self.model.load_model(model_tag='1024uint8')
        result = self.model.get_priority('CCCNc1ccccc1')
        expected = -0.45
        self.assertLess(abs(expected - result), 1e-4)

class TestMincost(unittest.TestCase):
    def setUp(self):
        self.model = mc.MinCostPrecursorPrioritizer()
        self.model.load_model()

    def test_01_get_priority(self):
        result = self.model.get_priority('CC(=O)N1C=C(C=C2N=C(N(N=CC3C=CC=CC=3)C(C)=O)N(C(C)=O)C2=O)C2=CC=CC=C21')
        expected = 4.01056028509
        self.assertLess(abs(expected - result), 1e-4)

    def test_02_get_priority(self):
        result = self.model.get_priority('CCCNc1ccccc1')
        expected = 0.0
        self.assertLess(abs(expected - result), 1e-4)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
