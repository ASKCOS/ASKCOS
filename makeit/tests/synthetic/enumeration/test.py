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

# BUG: This fails if the body of the outer loop is placed into a separate function
#      it seems like some state is being carried over between batch sizes...

def test_batch(ft, smiles, template_count, size):
    outcomes = []
    for start_at in range(0, template_count, size):
        outcomes.append(ft.get_outcomes(smiles, 100, start_at=start_at,
                                        end_at=start_at+size, template_prioritization=gc.popularity))
    unique_res = ForwardResult(smiles)
    for smiles, result in outcomes:
        unique_res.add_products(result.products)
    with open(os.path.join(os.path.dirname(__file__), 'expected/'+str(size)+'.pkl'), 'rb') as f:
        expected = pickle.load(f).get_products()
    result = unique_res.get_products()
    if len(result) != len(expected):
        return False
    for i in range(len(expected)):
        if expected[i].as_dict() != result[i].as_dict():
            return False
    return True

class TestTransformer(unittest.TestCase):
    def setUp(self):
        self.ft = transformer.ForwardTransformer(mincount=10)
        self.ft.load()

    def test_01_get_outcomes(self):
        template_count = self.ft.template_count()
        smiles = 'NC(=O)[C@H](CCC=O)N1C(=O)c2ccccc2C1=O'

        for batch_size in range(100, 1000, 100):
            test_batch(self.ft, smiles, template_count, batch_size)


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
