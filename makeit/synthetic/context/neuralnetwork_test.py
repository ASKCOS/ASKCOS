import os
import pickle
import unittest

import numpy as np

import makeit.global_config as gc
import makeit.synthetic.context.neuralnetwork as nn


class TestNeuralNetwork(unittest.TestCase):

    def test_01_get_n_conditions(self):
        """Test that NeuralNetContextRecommender.get_n_conditions gives the expected result."""
        cont = nn.NeuralNetContextRecommender()
        cont.load_nn_model(
            model_path=gc.NEURALNET_CONTEXT_REC['model_path'],
            info_path=gc.NEURALNET_CONTEXT_REC['info_path'],
            weights_path=gc.NEURALNET_CONTEXT_REC['weights_path']
        )
        result = cont.get_n_conditions('CC1(C)OBOC1(C)C.Cc1ccc(Br)cc1>>Cc1cccc(B2OC(C)(C)C(C)(C)O2)c1', 10,
                                       with_smiles=False, return_scores=True)

        with open(os.path.join(os.path.dirname(__file__), 'test_data/get_n_conditions.pkl'), 'rb') as t:
            expected = pickle.load(t, encoding='iso-8859-1')

        for e, r in zip(expected[0], result[0]):
            self.assertAlmostEqual(e[0], r[0], places=4)
            self.assertEqual(e[1:-2], r[1:-2])

        self.assertTrue(np.allclose(expected[1:], result[1:]))


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
