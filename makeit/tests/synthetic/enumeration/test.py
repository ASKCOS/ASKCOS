import makeit.synthetic.enumeration.transformer as transformer
from makeit.synthetic.enumeration.results import ForwardResult
import makeit.global_config as gc
import unittest
import os
import sys
is_py2 = sys.version[0] == '2'
if is_py2:
    import cPickle as pickle
else:
    import pickle as pickle


class TestTransformer(unittest.TestCase):
    def setUp(self):
        self.ft = transformer.ForwardTransformer()
        self.ft.load()

    def test_01_get_outcomes(self):
        smiles = 'NC(=O)[C@H](CCC=O)N1C(=O)c2ccccc2C1=O'
        res = self.ft.get_outcomes(smiles)
        self.assertEqual(type(res[0]), str)
        self.assertEqual(len(res[1]), 181)


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
