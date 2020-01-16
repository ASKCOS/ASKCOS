# import os, sys
# project_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from itertools import combinations
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, rdmolops
import time
from tqdm import tqdm
import socket, pickle


class ImpurityPredictor:

    def __init__(self, predictor, inspector, mapper,
                 topn_outcome=3, insp_threshold=0.75, celery_task=None, check_mapping=True):
        """

        :param predictor:
        :param inspector:
        :param mapper:
        :param topn_outcome: check top-n results
        :param insp_threshold:
        """
        self.host = 'localhost'
        self.predictor = predictor
        self.inspector = inspector
        self.mapper = mapper
        self.insp_threshold = insp_threshold  # probability threshold
        self.topn_outcome = topn_outcome
        self.n_top = 30  # this parameter is used for forward predictor, different from topn_outcome
        self.check_mapping = check_mapping
        self.task = celery_task

    def predict(self, reactants, reagents='', products='', solvents=''):
        """
        predict all possible products
        :param reactants: required inputs
        :param reagents: reagents will always be added in the prediction,
        if reagent is not specified, single reactant reaction prediction will not be considered
        reagents should not contribute atoms
        :param products:
        :param solvents:
        :return:
        """
        self.update_task(0, 'Impurity prediction started.')
        rct = reactants.split('.')
        #
        if not products:
            prd = []
        else:
            prd = products.split('.')

        # reagents are always considered a whole smiles
        if not reagents:
            rea = []
        else:
            rea = [reagents]

        if not solvents:
            sol = []
        else:
            sol = solvents.split('.')

        outcome_all = []  # predicted possible products
        # normal prediction, used as a baseline
        normal_outcomes = self.predictor(self.merge_smiles(rct + rea))
        # Mode 1: normal prediction
        outcome_all.extend(self.predict_M1(rct, rea, prd, sol))
        # Mode 2: over-reaction prediction
        outcome_all.extend(self.predict_M2(rct, rea, sol, outcome_all))
        # Mode 3: dimer prediction
        outcome_all.extend(self.predict_M3(rct, rea, sol, outcome_all))
        # Mode 4: solvent adduct prediction
        outcome_all.extend(self.predict_M4(rct, rea, sol, outcome_all))
        # Mode 5: subset of reactants available
        outcome_all.extend(self.predict_M5(rct, rea, sol, outcome_all))

        return {'predict_expand': self.merge_duplicates(outcome_all),
                'predict_normal': normal_outcomes}

    def predict_M1(self, rct, rea, prd, sol):
        """
        mode 1: normal prediction
        solvent is not included in the prediction since solvent will be considered invidually
        if prd is supplied, give prd as output
        """
        t0 = time.time()
        print('-----------Mode 1 normal prediction started-----------')
        output = []

        if not prd:
            self.update_task(0, 'Mode 1: normal prediction')
            outcomes = self.pred_insp(self.merge_smiles(rct + rea))
            if outcomes:
                for outcome in outcomes:
                    output.append({'rct_smiles': self.merge_smiles(rct),
                                   'prd_smiles': outcome['smiles'],
                                   'rea_smiles': rea,
                                   'sol_smiles': sol,
                                   'mode': 1,
                                   'insp_score': outcome['insp_score']})
        else:
            for item in prd:
                output.append({'rct_smiles': self.merge_smiles(rct),
                               'prd_smiles': Chem.MolToSmiles(Chem.MolFromSmiles(item)),  # canonize the smiles
                               'rea_smiles': rea,
                               'sol_smiles': sol,
                               'mode': 1,
                               'insp_score': 1})
        print('Mode 1 prediction finished: {} s.'.format(time.time() - t0))

        # use the rank 1 prediction from Mode 1 or first provided major product for similarity calculations
        self.major_product = output[0]['prd_smiles']
        self.update_task(1, 'Mode 1: normal prediction')
        return output

    def predict_M2(self, rct, rea, sol, outcome_all):
        """
        mode 2: over-reaction of products with all or partial reactants
        this step will make the predicted products are over-reaction products by using atommapping
        some dimer formation can also be included in this category if the dimer is formed by in this mode: 2A + B --> A-B-A
        """
        t0 = time.time()
        print('-----------Mode 2 over-reaction started-----------')
        output = []

        prd_for_overreaction = []
        # get a list of products for examinating overreaction
        if outcome_all:
            for outcome in outcome_all:
                if outcome['mode'] == 1:
                    prd_for_overreaction.append(outcome['prd_smiles'])

        # if reagent is empty, consider at least one component from reactants
        # if reagent isn't empty, can consider product and reagent with anything from reactants
        if not rea:
            rct_comb = self.smiles_comb(rct, 1, len(rct))
        else:
            rct_comb = self.smiles_comb(rct, 0, len(rct))

        total = len(rct_comb) * len(prd_for_overreaction)
        count = 0
        for item in tqdm(rct_comb):
            for prd in prd_for_overreaction:
                rsmi = self.merge_smiles(item + rea + [prd])
                self.update_task(count/total, 'Mode 2: over-reaction')
                for outcome in self.pred_insp(rsmi):
                    pred_psmi = outcome['smiles']
                    # check whether the prd is present in the overreaction prediction
                    # only consider the prediction that is consistent with prediction mode
                    if self.check_mode_outcome(rsmi, pred_psmi, [prd]):
                        output.append({'rct_smiles': self.merge_smiles(item + [prd]),
                                       'prd_smiles': pred_psmi,
                                       'rea_smiles': rea,
                                       'sol_smiles': sol,
                                       'mode': 2,
                                       'insp_score': outcome['insp_score']})
                count += 1
        print('Mode 2 prediction finished: {} s.'.format(time.time() - t0))
        self.update_task(1, 'Mode 2: over-reaction')
        return output

    def predict_M3(self, rct, rea, sol, outcome_all):
        """
        mode 3: dimer formations
        1. solvents should be excluded from dimer formaitons
        2. products dimer can also be considered
        here, only products from mode 1 are considered for dimer formation

        inspector is problematic when checking dimer formation, since inputing dimer and mononer give
        yield the same fingerprints
        """
        t0 = time.time()
        print('-----------Mode 3 dimerization started-----------')
        output = []
        monomer_list = []
        # create the monomer list for checking dimer formation
        for item in rct:
            monomer_list.append(item)
        # only consider the dimer formation from normal predicitions
        for outcome in outcome_all:
            if outcome['mode'] == 1:
                monomer_list.append(outcome['prd_smiles'])

        total = len(monomer_list)
        count = 0
        for monomer in tqdm(monomer_list):
            rsmi = self.merge_smiles([monomer, monomer] + rea)
            self.update_task(count/total, 'Mode 3: dimerization')
            for outcome in self.pred_insp(rsmi):
                pred_psmi = outcome['smiles']
                # check whether two monomers are present in the dimer prediction
                # only consider the prediction that is consistent with prediction mode
                if self.check_mode_outcome(rsmi, pred_psmi, [monomer, monomer]):
                    output.append({'rct_smiles': self.merge_smiles([monomer, monomer]),
                                   'prd_smiles': pred_psmi,
                                   'rea_smiles': rea,
                                   'sol_smiles': sol,
                                   'mode': 3,
                                   'insp_score': outcome['insp_score']})
            count += 1
        print('Mode 3 prediction finished: {} s.'.format(time.time() - t0))
        self.update_task(1, 'Mode 3: dimerization')
        return output

    def predict_M4(self, rct, rea, sol, outcome_all):
        """
        mode 4: solvent adduct
        only consider solvent adduct for individual
        product solvent adduct are also considered
        """
        t0 = time.time()
        print('-----------Mode 4 solvent adduct started-----------')
        output = []

        if sol:
            smiles_to_check = []
            for item in rct:
                smiles_to_check.append(item)

            for outcome in outcome_all:
                if outcome['mode'] == 1:
                    smiles_to_check.append(outcome['prd_smiles'])

            total = len(smiles_to_check) * len(sol)
            count = 0
            for smiles in tqdm(smiles_to_check):
                for sol_smiles in sol:
                    self.update_task(count / total, 'Mode 4: solvent adduct')
                    rsmi = self.merge_smiles([smiles, sol_smiles] + rea)
                    for outcome in self.pred_insp(rsmi):
                        pred_psmi = outcome['smiles']
                        # check whether the solvent is present in the solvent adduct prediction
                        # only consider the prediction that is consistent with prediction mode
                        if self.check_mode_outcome(rsmi, pred_psmi, [sol_smiles]):
                            output.append({'rct_smiles': self.merge_smiles([smiles, sol_smiles]),
                                           'prd_smiles': outcome['smiles'],
                                           'rea_smiles': rea,
                                           'sol_smiles': sol,
                                           'mode': 4,
                                           'insp_score': outcome['insp_score']})
                    count += 1
        print('Mode 4 prediction finished: {} s.'.format(time.time() - t0))
        self.update_task(1, 'Mode 4: solvent adduct')
        return output

    def predict_M5(self, rct, rea, sol, outcome_all):
        """
        mode 5: subset of reactants
        this is a little tricky since it might generate many duplicated outcomes as M1
        solvent adducts are not considered here
        only give non-duplicated results as the output
        """
        t0 = time.time()
        print('-----------Mode 5 subset prediction started-----------')
        output = []
        # if reagent is empty, consider at least two component from reactants
        # if reagent isn't empty, can consider one component from reactants
        # this is actually a harsh requirement for this, since some pistachio data is recorded as two reactants, but
        # without reagents
        if not rea:
            rct_comb = self.smiles_comb(rct, 2, len(rct) - 1)  # don't consider all reactants
            # rct_comb = self.smiles_comb(rct, 2, len(rct) - 1)  # temporary changed to this to include
        else:
            rct_comb = self.smiles_comb(rct, 1, len(rct) - 1)  # don't consider all reactants

        previous_outcome_list = [previous_outcome['prd_smiles'] for previous_outcome in outcome_all]

        total = len(rct_comb)
        count = 0
        for item in tqdm(rct_comb):
            self.update_task(count / total, 'Mode 5: subset prediction')
            for outcome in self.pred_insp(self.merge_smiles(item + rea)):
                # remove duplicated outcomes
                if outcome['smiles'] not in previous_outcome_list:
                    output.append({'rct_smiles': self.merge_smiles(item),
                                   'prd_smiles': outcome['smiles'],
                                   'rea_smiles': rea,
                                   'sol_smiles': sol,
                                   'mode': 5,
                                   'insp_score': outcome['insp_score']})
            count += 1

        print('Mode 5 prediction finished: {} s.'.format(time.time() - t0))
        self.update_task(1, 'Mode 5: subset prediction')
        return output

    def pred_insp(self, rsmi):
        # predict the outcomes with forward predictor
        outcomes = self.predictor(rsmi)
        # check the outcomes with reaction inspector
        inspected_outcomes = []
        if outcomes:
            for i in range(min(self.topn_outcome, len(outcomes[0]))):
                psmi = outcomes[0][i]['outcome']['smiles']

                insp_score = self.inspector(rsmi + '>>' + psmi)
                if insp_score > self.insp_threshold:
                    inspected_outcomes.append({'smiles': psmi,
                                               'pred_rank': outcomes[0][i]['rank'],
                                               'insp_score': insp_score})
        return inspected_outcomes

    def check_mode_outcome(self, rsmi, psmi, sub_rct_to_check):
        """
        This is used to make sure that the predicted outcome is whthin the specified mode of the reaction
        for example, the dimer formation reaction should have two monomers present in the molecule.
        :param rsmi:
        :param psmi: only one product allowed
        :param sub_rct_to_check: it is a list of smiles that we want to check
        :return: whether the desired compounds are present in the products
        """

        def extract_am_num(smi_am):
            """return a set of atom mapping numbers, smi must have atom mapping numbers"""
            mol = Chem.MolFromSmiles(smi_am)
            amnum_list = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
            return amnum_list

        def remove_am(smi_am):
            mol = Chem.MolFromSmiles(smi_am)
            for atom in mol.GetAtoms():
                # get rid of clearing atom mapping function, so the output will have atom maps
                atom.ClearProp('molAtomMapNumber')
            return Chem.MolToSmiles(mol)

        if self.check_mapping:
            # find the atom mapping
            # print(rsmi + '>>' + psmi)
            [rsmi_am, psmi_am] = self.mapper(rsmi + '>>' + psmi).split('>>')
            if psmi_am:
                sub_rct_am = []
                rsmi_am_list = rsmi_am.split('.')

                # find the atom mapping for the interested smiles
                for sub_rct in sub_rct_to_check:
                    for i, rct_am in enumerate(rsmi_am_list):
                        if sub_rct == remove_am(rct_am):
                            sub_rct_am.append(rct_am)
                            rsmi_am_list.pop(i)
                            break

                # get all the atom numbers in the products as a set
                prd_amnum_list = extract_am_num(psmi_am)
                # check whether the atoms are present in the product
                check_present = True
                for smi_am in sub_rct_am:
                    count_present = 0
                    amnum_list = extract_am_num(smi_am)
                    for amnum in amnum_list:
                        if amnum in prd_amnum_list:
                            count_present += 1
                    # if 30% atoms are present in the products, then we consider this substructure is present in the product
                    if count_present / len(amnum_list) > 0.3:
                        check_present = check_present and True
                    else:
                        check_present = check_present and False
                        break
                return check_present
            else:
                return False
        else:
            return True

    def smiles_comb(self, smiles_list, min_num_comb, max_num_comb):
        """
        provide a complete combination of smiles list combos
        """
        # max_len = len(smiles_list)
        smiles_comb_list = []
        for i in range(min_num_comb, max_num_comb + 1):
            smiles_comb_list.extend(combinations(smiles_list, i))
        output = []
        for comb in smiles_comb_list:
            output.append(list(comb))
        return output

    def merge_smiles(self, smiles_list):
        merged_smiles = ''
        for smiles in smiles_list:
            merged_smiles = merged_smiles + '.' + smiles
        return merged_smiles[1::]

    def merge_duplicates(self, outcome_all):
        reaction_mode_name = {1: 'Normal reaction',
                              2: 'Over-reaction',
                              3: 'Dimerization',
                              4: 'Solvent adducts',
                              5: 'Subset of reactants'}
        record = {}
        for i, outcome in enumerate(outcome_all):
            if outcome['prd_smiles'] in record:
                record[outcome['prd_smiles']].append(i)
            else:
                record[outcome['prd_smiles']] = [i]

        new_outcome_all = []
        count = 1

        # calculate similarity and sort
        all_prd_smiles = [prd_smiles for prd_smiles in record]
        all_prd_similarity = [self.calculate_similarity(smi) for smi in all_prd_smiles]
        sorted_inedx = sorted(range(len(all_prd_similarity)),
                              key=lambda k: all_prd_similarity[k],
                              reverse=True)
        for i in sorted_inedx:
            prd_smiles = all_prd_smiles[i]
            rct_rea_sol = []
            modes = []
            insp_scores = []
            for index in record[prd_smiles]:
                rct_rea_sol.append({'rct_smiles': outcome_all[index]['rct_smiles'],
                                    'rea_smiles': outcome_all[index]['rea_smiles'],
                                    'sol_smiles': outcome_all[index]['sol_smiles']})
                modes.append(outcome_all[index]['mode'])
                insp_scores.append(outcome_all[index]['insp_score'])
            modes = list(set(modes))
            modes.sort()
            new_outcome_all.append({'no': count,
                                    'prd_smiles': prd_smiles,
                                    'rct_rea_sol': rct_rea_sol,
                                    'insp_scores': insp_scores,
                                    'avg_insp_score': sum(insp_scores)/len(insp_scores),
                                    'similarity_to_major': all_prd_similarity[i],
                                    'modes': modes,
                                    'modes_name': ', '.join([reaction_mode_name[i] for i in modes])})
            count += 1


        return new_outcome_all

    def calculate_similarity(self, smiles):
        fp_major_prod = Chem.RDKFingerprint(Chem.MolFromSmiles(self.major_product))
        fp = Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))
        return DataStructs.FingerprintSimilarity(fp_major_prod, fp)

    def update_task(self, percent, message):
        if self.task:
            self.task.update_state(
                state='running',
                meta={
                    'percent': percent,
                    'message': message
                })


if __name__ == '__main__':
    from makeit.synthetic.evaluation.template_free import TemplateFreeNeuralNetScorer
    from makeit.synthetic.evaluation.fast_filter import FastFilerImpurityInspector
    from makeit.synthetic.atom_mapper.wln_mapper import WLN_AtomMapper

    #
    predictor = TemplateFreeNeuralNetScorer().evaluate
    inspector = FastFilerImpurityInspector().evaluate
    mapper = WLN_AtomMapper().evaluate
    # %%
    pe = ImpurityPredictor(predictor, inspector, mapper, 1, 0.75, check_mapping=True)

    rct_smi = 'CC(C)(C)OC(=O)N1CCC(N)CC1.O=C(O)c1cccc2oc(-c3ccccc3)nc12'

    prd_smi = ''
    sol_smi = ''
    rea_smi = ''
    outcomes = pe.predict(rct_smi,
                          products=prd_smi,
                          reagents=rea_smi,
                          solvents=sol_smi)
