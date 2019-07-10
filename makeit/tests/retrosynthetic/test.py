import unittest
import makeit.retrosynthetic.tree_builder as retro_tree
import makeit.retrosynthetic.transformer as retro_trans
import makeit.global_config as gc

class TestTreeBuilder(unittest.TestCase):
    def test_01_get_buyable_paths(self):
        treeBuilder = retro_tree.TreeBuilder(celery=False, mincount=25, mincount_chiral=10)
        status_result, paths_result = treeBuilder.get_buyable_paths('CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2', max_depth=4, template_prioritization=gc.relevance,
                                            precursor_prioritization=gc.relevanceheuristic, nproc=2, expansion_time=60, max_trees=500, max_ppg=10,
                                            max_branching=25, apply_fast_filter=True, filter_threshold=0.75,
                                            min_chemical_history_dict={'as_reactant':5, 'as_product':1, 'logic':'none'})
        expected_status = (28, 25, {0: 1, 0.5: 25, 1: 27})
        expected_paths = [{'smiles': 'CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2', 'ppg': 6.0, 'as_product': 17, 'as_reactant': 36, 'children': [], 'is_chemical': True, 'id': 1}, {'smiles': 'CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2', 'ppg': 6.0, 'as_product': 17, 'as_reactant': 36, 'children': [{'smiles': 'CN1C2CCC1CC(O)C2.O=C(O)C(CO)c1ccccc1>>CN1C2CCC1CC(OC(=O)C(CO)c1ccccc1)C2', 'is_reaction': True, 'necessary_reagent': u'', 'children': [{'smiles': 'CN1C2CCC1CC(O)C2', 'ppg': 1.0, 'as_product': 61, 'as_reactant': 112, 'children': [], 'is_chemical': True, 'id': 3}, {'smiles': 'O=C(O)C(CO)c1ccccc1', 'ppg': 3.0, 'as_product': 75, 'as_reactant': 265, 'children': [], 'is_chemical': True, 'id': 4}], 'plausibility': 0.9999294281005859, 'template_score': 0.01923198811709881, 'score': -0.2079868173610025, 'num_examples': 19578, 'tforms': ['59c5118c05581eb9f5753c8c', '59c5118c05581eb9f5753c9b', '59c511b905581eb9f5756604'], 'id': 2, 'template': '59c5118c05581eb9f5753c8c'}], 'is_chemical': True, 'id': 1}]
        self.assertEqual(expected_status, status_result)
        self.assertEqual(expected_paths, paths_result)

class TestTransformer(unittest.TestCase):
    def test_01_get_outcomes(self):
        t = retro_trans.RetroTransformer()
        t.load(chiral=True, refs=False, rxns=True)
        outcomes = t.get_outcomes('CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1', 100, (gc.relevanceheuristic, gc.relevance))
        precursors = outcomes.precursors
        precursors_result = [precursor.smiles_list for precursor in precursors]
        expected_precursors = [['C=CC(=O)N1[C@@H](c2ccccc2)CC[C@@H]1c1ccccc1', 'CCOC(=O)C/N=C/c1ccccc1']]
        self.assertEqual(expected_precursors, precursors_result)

    def test_02_apply_one_template_by_idx(self):
        t = retro_trans.RetroTransformer()
        t.load(chiral=True, refs=False, rxns=True)
        outcomes_result = t.apply_one_template_by_idx(1, 'CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1', 109659)
        expected_outcomes = [(1, 'CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1', 109659, [('C=CC(=O)N1[C@@H](c2ccccc2)CC[C@@H]1c1ccccc1', [0.2523617446422577, 0.1969166249036789, 0.15234193205833435, 0.13455800712108612, 0.030202925205230713, 0.024427464231848717, 0.023061120882630348, 0.018903391435742378, 0.017541231587529182, 0.01372396107763052, 0.013607553206384182, 0.012942434288561344, 0.01270198542624712, 0.012236274778842926, 0.00939953699707985, 0.008472910150885582, 0.0059701986610889435, 0.0048681264743208885, 0.004287674557417631, 0.0042765019461512566, 0.004114005248993635, 0.003224626649171114, 0.003134896280243993, 0.00310404016636312, 0.002972274785861373, 0.002696158131584525, 0.0026850011199712753, 0.0021682260558009148, 0.002018335508182645, 0.001683939597569406, 0.0014018624788150191, 0.0012026543263345957, 0.001027433667331934, 0.0010042842477560043, 0.0009036235860548913, 0.0008185449405573308, 0.0006995245930738747, 0.0005884007550776005, 0.0005666478537023067, 0.000491276616230607, 0.0004684898303821683, 0.00045587049680761993, 0.00043441925663501024, 0.0003655038308352232, 0.0003575180599000305, 0.0003510901879053563, 0.00034910428803414106, 0.00030314744799397886, 0.00029094854835420847, 0.000254554208368063, 0.00024915405083447695, 0.00022937390895094723, 0.0001904226082842797, 0.0001790589012671262, 0.00017418322386220098, 0.00016860727919265628, 0.00015959763550199568, 0.00015663796511944383, 0.00015368156891781837, 0.0001446626556571573, 0.00014138488040771335, 0.00013653609494213015], [922, 33433, 11500, 7263, 4005, 3581, 5541, 150347, 62514, 54396, 10503, 115523, 26842, 67342, 39406, 7421, 121814, 44003, 104591, 54496, 12802, 17425, 29358, 40183, 89253, 17054, 17142, 14952, 4909, 155144, 73901, 95665, 3419, 17528, 155365, 47886, 31303, 974, 70550, 50432, 97960, 36674, 114062, 108, 6776, 21165, 130595, 81282, 7635, 52, 82107, 95674, 1065, 74925, 88017, 92158, 37880, 9606, 6613, 49535, 16476, 94670], 1), ('CCOC(=O)C/N=C/c1ccccc1', [0.9659425020217896, 0.01252610795199871, 0.005059113260358572, 0.0038746828213334084, 0.0022534450981765985, 0.001429172931239009, 0.0008012531325221062, 0.0007077240734361112, 0.000661780999507755, 0.0006321603432297707, 0.0006210291176103055, 0.0005571309011429548], [1287, 186, 128290, 2120, 129604, 1910, 17747, 70579, 303, 713, 113217, 8087], 1)], 0.9996535778045654)]
        self.assertEqual(expected_outcomes, outcomes_result)


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
