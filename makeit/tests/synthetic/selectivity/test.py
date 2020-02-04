import unittest

import tensorflow.compat.v1 as tf

from makeit.synthetic.selectivity.site_selectivity import Site_Predictor


class SiteSelectivity(unittest.TestCase):

    def setUp(self):
        """This method is run once before each test in this class."""
        # Clear tensorflow sessions - leftover sessions from other models cause unexpected errors
        tf.keras.backend.clear_session()

    def test_01_predict(self):
        """Test that the Site_Predictor works as expected."""
        react = 'Cc1ccccc1'
        predictor = Site_Predictor()
        res = predictor.predict(react)
        self.assertEqual(len(res), 123)
        self.assertEqual(type(res[0]), dict)
        self.assertEqual(len(res[0].get('atom_scores')), 7)


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
