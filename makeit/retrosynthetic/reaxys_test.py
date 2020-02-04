import time
import unittest
from makeit import global_config as gc
from makeit.retrosynthetic.transformer import RetroTransformer
from makeit.prioritization.templates.relevance import RelevanceTemplatePrioritizer

class TestTransformer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        template_prioritizer = RelevanceTemplatePrioritizer()
        template_prioritizer.load_model(
            gc.RELEVANCE_TEMPLATE_PRIORITIZATION['reaxys']['model_path']
        )
        cls.transformer = RetroTransformer(
            load_all=False, use_db=False, template_set='reaxys',
            template_prioritizer=template_prioritizer
        )

    def setUp(self):
        self._started_at = time.time()

    def tearDown(self):
        elapsed = time.time() - self._started_at
        print('{} ({}s)'.format(self.id(), round(elapsed, 2)))

    def test_01_load(self):
        self.transformer.load()
        self.assertEqual(len(self.transformer.templates), gc.RELEVANCE_TEMPLATE_PRIORITIZATION['reaxys']['output_size'])
        self.assertIsInstance(self.transformer.templates[0], dict)
        self.assertNotEqual(self.transformer.templates[0].get('index'), None)
        self.assertNotEqual(self.transformer.templates[0].get('reaction_smarts'), None)
        self.assertNotEqual(self.transformer.templates[0].get('count'), None)

    def test_02_get_outcomes(self):
        outcomes = self.transformer.get_outcomes('CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1')
        self.assertNotEqual(outcomes[0].get('smiles_split'), None)

    def test_03_apply_one_template_by_idx(self):
        worker_id = 1
        template_idx = 109659
        smiles = 'CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1'
        max_cum_prob=0.995
        max_num_templates=100
        outcomes_result = self.transformer.apply_one_template_by_idx(
            worker_id, smiles, template_idx, template_set='reaxys',
            max_cum_prob=max_cum_prob, max_num_templates=max_num_templates
        )
        self.assertEqual(worker_id, outcomes_result[0][0])
        self.assertEqual(smiles, outcomes_result[0][1])
        self.assertEqual(template_idx, outcomes_result[0][2])
        result = outcomes_result[0][3][0]
        self.assertIsInstance(result[0], str)
        self.assertIsInstance(result[1][0], float)
        self.assertIsInstance(result[2][0], int)
        self.assertLessEqual(len(result[1]), max_num_templates)
        self.assertLessEqual(sum(result[1]), max_cum_prob)
        self.assertIsInstance(outcomes_result[0][4], float)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
