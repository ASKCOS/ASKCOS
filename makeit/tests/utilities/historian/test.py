import unittest

from makeit.utilities.historian.chemicals import ChemHistorian


class TestChemHistorian(unittest.TestCase):

    def setUp(self):
        """This method is run once before each test in this class."""
        self.chemhistorian = ChemHistorian(hashed=True)
        self.chemhistorian.load()

    def test_01_lookup_smiles(self):
        """Test that we can look up a SMILES string in chemhistorian."""
        result = self.chemhistorian.lookup_smiles('CCCCO')
        expected = {'as_product': 2726, 'as_reactant': 17450, 'template_set': 'reaxys'}
        self.assertEqual(expected, result)


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
