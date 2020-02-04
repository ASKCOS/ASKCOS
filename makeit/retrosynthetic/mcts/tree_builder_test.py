import unittest
import makeit.retrosynthetic.mcts.tree_builder as mcts_tree
import makeit.retrosynthetic.mcts.rl_model
import makeit.global_config as gc
import rdkit.Chem as Chem
import os
import sys
is_py2 = sys.version[0] == '2'
if is_py2:
    import cPickle as pickle
else:
    import pickle as pickle

@unittest.skip("Non-deterministic")
class TestMCTSTreeBuilder(unittest.TestCase):
    def setUp(self):
        self.NCPUS = 4
        self.Tree = mcts_tree.MCTS(nproc=self.NCPUS, mincount=gc.RETRO_TRANSFORMS_CHIRAL['mincount'],
            mincount_chiral=gc.RETRO_TRANSFORMS_CHIRAL['mincount_chiral'],
            celery=False)

    def test_01_scopolamine_test(self):
        smiles = 'Cc1ncc([N+](=O)[O-])n1CC(C)O'
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), True)
        status_result, paths_result = self.Tree.get_buyable_paths(smiles,
                                            nproc=self.NCPUS,
                                            expansion_time=30,
                                            max_cum_template_prob=0.995,
                                            template_count=100,
                                            # min_chemical_history_dict={'as_reactant':5, 'as_product':5,'logic':'none'},
                                            soft_reset=False,
                                            soft_stop=True)
        with open(os.path.join(os.path.dirname(__file__), 'test_data/test_01_scopolamine_test.pkl'), 'rb') as t:
            expected_status, expected_paths = pickle.load(t)
        self.assertEqual(expected_status, status_result)
        self.assertEqual(expected_paths, paths_result)

    def test_02_get_buyable_paths(self):
        smiles = 'CCCCCN(CCCCC)CCCC(=O)OCCC'
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), True)
        status_result, paths_result = self.Tree.get_buyable_paths(smiles,
                                            nproc=self.NCPUS,
                                            expansion_time=30,
                                            max_cum_template_prob=0.995,
                                            template_count=100,
                                            soft_reset=False,
                                            soft_stop=True)
        with open(os.path.join(os.path.dirname(__file__), 'test_data/test_02_get_buyable_paths.pkl'), 'rb') as t:
            expected_status, expected_paths = pickle.load(t)
        self.assertEqual(expected_status, status_result)
        self.assertEqual(expected_paths, paths_result)

# class TestRlModel(unittest.TestCase):
#     def test_01_(self):
#         model = RelevanceTemplatePrioritizer(use_tf=True)
#         model.load_model()
#         smis = ['CCCOCCC', 'CCCNc1ccccc1']
#         result = []
#         for smi in smis:
#             lst = model.get_topk_from_smi(smi)
#             result.append((smi, lst))
#
#         expected = []
#         self.assertEqual(expected, result)

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
