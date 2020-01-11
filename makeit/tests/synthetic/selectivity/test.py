from makeit.synthetic.selectivity.site_selectivity import Site_Predictor

import unittest

class SiteSelectivity(unittest.TestCase):
    def test_01_predict(self):
        react = 'Cc1ccccc1'
        predictor = Site_Predictor()
        res = predictor.predict(react)
        self.assertEqual(len(res), 123)
        self.assertEqual(type(res[0]), dict)
        self.assertEqual(len(res[0].get('atom_scores')), 7)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
