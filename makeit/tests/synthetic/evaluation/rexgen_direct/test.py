import makeit.synthetic.evaluation.rexgen_direct.predict as p
import unittest
import os
import sys
is_py2 = sys.version[0] == '2'
if is_py2:
    import cPickle as pickle
else:
    import pickle as pickle

class TestPredict(unittest.TestCase):
    def test_01_predict(self):
        tffp = p.TFFP()
        result = tffp.predict('CCCO.CCCBr')
        with open(os.path.join(os.path.dirname(__file__), 'expected/predict.pkl'), 'rb') as t:
            expected = pickle.load(t)
        self.assertEqual(expected, result)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
