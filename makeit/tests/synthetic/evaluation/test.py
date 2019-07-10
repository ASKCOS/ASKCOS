import makeit.synthetic.evaluation.template_free as tf
import makeit.synthetic.evaluation.fast_filter as ff
import makeit.synthetic.evaluation.evaluator as ev
import makeit.global_config as gc
import unittest
import os
import sys
is_py2 = sys.version[0] == '2'
if is_py2:
    import cPickle as pickle
else:
    import pickle as pickle

class TestTemplateFree(unittest.TestCase):
    def test_01_evaluate(self):
        react = 'CCCCO.CCCCBr'
        scorer = tf.TemplateFreeNeuralNetScorer()
        result = scorer.evaluate(react)
        with open(os.path.join(os.path.dirname(__file__), 'expected/template_free.pkl'), 'rb') as t:
            expected = pickle.load(t)
        self.assertEqual(expected, result)

class TestFastFilter(unittest.TestCase):
    def setUp(self):
        self.model = ff.FastFilterScorer()
        self.model.load(model_path=gc.FAST_FILTER_MODEL['trained_model_path'])

    def test_01_evaluate(self):
        result = self.model.evaluate('CCO.CC(=O)O', 'CCOC(=O)C')
        expected = [[{'outcome': {'smiles': 'CCOC(=O)C', 'template_ids': [], 'num_examples': 0}, 'score': 0.9789425730705261, 'rank': 1.0, 'prob': 0.9789425730705261}]]
        self.assertEqual(expected, result)

    def test_02_evaluate(self):
        result = self.model.evaluate('[CH3:1][C:2](=[O:3])[O:4][CH:5]1[CH:6]([O:7][C:8]([CH3:9])=[O:10])[CH:11]([CH2:12][O:13][C:14]([CH3:15])=[O:16])[O:17][CH:18]([O:19][CH2:20][CH2:21][CH2:22][CH2:23][CH2:24][CH2:25][CH2:26][CH2:27][CH2:28][CH3:29])[CH:30]1[O:31][C:32]([CH3:33])=[O:34].[CH3:35][O-:36].[CH3:38][OH:39].[Na+:37]', 'CCCCCCCCCCOC1OC(CO)C(O)C(O)C1O')
        expected = [[{'outcome': {'smiles': 'CCCCCCCCCCOC1OC(CO)C(O)C(O)C1O', 'template_ids': [], 'num_examples': 0}, 'score': 0.9983256459236145, 'rank': 1.0, 'prob': 0.9983256459236145}]]
        self.assertEqual(expected, result)

    def test_03_evaluate(self):
        result = self.model.evaluate('CNC.Cc1ccc(S(=O)(=O)OCCOC(c2ccccc2)c2ccccc2)cc1', 'CN(C)CCOC(c1ccccc1)c2ccccc2')
        expected = [[{'outcome': {'smiles': 'CN(C)CCOC(c1ccccc1)c2ccccc2', 'template_ids': [], 'num_examples': 0}, 'score': 0.9968607425689697, 'rank': 1.0, 'prob': 0.9968607425689697}]]
        self.assertEqual(expected, result)

    def test_04_filter_with_threshold(self):
        flag_result, score_result = self.model.filter_with_threshold('CCO.CC(=O)O', 'CCOC(=O)C', 0.75)
        expected_flag = [[ True]]
        expected_score = 0.978942573071
        self.assertEqual(expected_flag, flag_result)
        self.assertAlmostEqual(expected_score, score_result)

class TestEvaluator(unittest.TestCase):
    def test_01_evaluate(self):
        evaluator = ev.Evaluator(celery=False)
        result = evaluator.evaluate('CCCCO.CCCCBr', 'O=C1CCCCCCCO1', [(20,'','','','','')], forward_scorer=gc.templatefree,
            return_all_outcomes=True)
        with open(os.path.join(os.path.dirname(__file__), 'expected/evaluator.pkl'), 'rb') as t:
            expected = pickle.load(t)
        self.assertEqual(expected, result)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
