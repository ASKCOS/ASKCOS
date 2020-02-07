import unittest
import makeit.utilities.buyable.pricer as p

class TestPricer(unittest.TestCase):
    def setUp(self):
        self.pricer = p.Pricer()
        self.pricer.load()

    def test_01_lookup_smiles(self):
        result = self.pricer.lookup_smiles('CCCCCO')
        expected = 1.0
        self.assertAlmostEqual(expected, result)

    def test_02_lookup_smiles(self):
        result = self.pricer.lookup_smiles('CCCCXCCO')
        expected = 0.0
        self.assertAlmostEqual(expected, result)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
