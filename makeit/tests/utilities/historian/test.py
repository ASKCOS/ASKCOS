import unittest
import makeit.utilities.historian.chemicals as c

class TestChemHistorian(unittest.TestCase):
    def setUp(self):
        self.chemhistorian = c.ChemHistorian()
        self.chemhistorian.load_from_file(refs=False, compressed=True)

    def test_01_lookup_smiles(self):
        result = self.chemhistorian.lookup_smiles('CCCCO')
        expected = {'as_reactant_refs': [], 'as_product_refs': [], 'as_reactant': 17450, 'as_product': 2726}
        self.assertEqual(expected, result)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
