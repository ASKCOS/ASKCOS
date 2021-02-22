from rdkit.Chem import rdmolops, Descriptors
import makeit.global_config as gc
import rdkit.Chem as Chem
import numpy as np
from makeit.utilities.io.logger import MyLogger
template_free_scorer_loc = 'template_free_scorer'


class WLN_AtomMapper:
    """Template-free neural net evaluator.

        Attributes:
            model (makeit.synthetic.evaluation.rexgen_direct.predict.TFFP):
                Template-free forward predictor.
        """

    def __init__(self, **kwargs):
        """Initializes TemplateFreeNeuralNetScorer.

        Args:
            **kwargs: Unused.
        """
        from makeit.synthetic.evaluation.rexgen_direct.predict import TFFP
        self.model = TFFP()

    def get_wln_resutls(self, reactants_smiles, contexts=[(20,'','','','','')], **kwargs):
        """Evaluates possible reaction outcomes for given reactants.

        Evaluation does not use context, but left as dummy pos var.

        Args:
            reactants_smiles (list of str??): SMILES string of reactants.
            contexts (list, optional): Unused.
                (default: {[(20,'','','','','')]})
            **kwargs: Additional optional parameters.

        Returns:
            list: Predicted reaction outcomes.
        """

        all_outcomes = []
        all_mapped_rct_smiles = []
        for (T1, slvt1, rgt1, cat1, t1, y1) in contexts:
            this_reactants_smiles = reactants_smiles
            if slvt1:
                this_reactants_smiles += '.' + slvt1
            if rgt1:
                this_reactants_smiles += '.' + rgt1
            if cat1:
                this_reactants_smiles += '.' + cat1

            mapped_rct_smiles, outcomes = self.model.predict(this_reactants_smiles, top_n=1e10, atommap=True)
            if not outcomes:
                all_outcomes.append([{
                    'rank': 1,
                    'outcome': {'smiles': '', 'template_ids':[], 'num_examples':0},
                    'score': 0,
                    'prob': 0,
                    }])
                continue

            outcomes_to_ret = {}
            reactants_smiles_split = this_reactants_smiles.split('.')
            for outcome in outcomes:
                smiles_list = set(self.remove_am(outcome['smiles']).split('.'))

                # Canonicalize
                smiles_canonical = set()
                for smi in smiles_list:
                    mol = Chem.MolFromSmiles(smi)
                    if not mol:
                        MyLogger.print_and_log('Template free evaluator could not reparse product {}'.format('.'.join(smiles_list)),
                        template_free_scorer_loc, 1)
                        continue
                    smiles_canonical.add(Chem.MolToSmiles(mol))


                # Remove unreacted frags
                smiles_canonical = smiles_canonical - set(reactants_smiles_split)
                if not smiles_canonical:
                    continue # no reaction?

                smiles = max(smiles_canonical, key=len) # NOTE: this is not great...byproducts may be longer

                if not smiles:
                    continue

                smiles_w_am = ''
                for item in outcome['smiles'].split('.'):
                    # print(item)
                    try:
                        mol = Chem.MolFromSmiles(self.remove_am(item))
                        # print(Chem.MolToSmiles(mol))
                        if Chem.MolToSmiles(mol) == smiles:
                            smiles_w_am = item
                            break
                    except Exception as e:
                        print(e)
                        pass

                if smiles in outcomes_to_ret:
                    outcomes_to_ret[smiles]['rank'] = min(outcomes_to_ret[smiles]['rank'], outcome['rank'])
                    outcomes_to_ret[smiles]['score'] = np.log(np.exp(outcomes_to_ret[smiles]['score']) + np.exp(outcome['score']))
                    outcomes_to_ret[smiles]['prob'] += outcome['prob']
                else:
                    # Append outcome information
                    outcomes_to_ret[smiles] = {
                        'rank': outcome['rank'],
                        'outcome': {
                            'smiles': smiles,
                            'smiles_w_am': smiles_w_am,
                            'template_ids': [],
                            'num_examples': 0,
                        },
                        'score': float(outcome['score']),
                        'prob': float(outcome['prob']),
                        'mol_wt': float(Descriptors.MolWt(Chem.MolFromSmiles(smiles)))
                    }

            # Renormalize and re-rank
            outcomes = sorted(outcomes_to_ret.values(), key=lambda x: x['prob'], reverse=True)
            total_prob = sum([outcome['prob'] for outcome in outcomes])
            for i, outcome in enumerate(outcomes):
                outcomes[i]['rank'] = i + 1
                outcomes[i]['prob'] = outcome['prob'] / total_prob

            all_outcomes.append(outcomes)
            all_mapped_rct_smiles.append(mapped_rct_smiles)
        # Return in lists as if we received a list of contexts
        return all_mapped_rct_smiles, all_outcomes

    def evaluate(self, rxnsmi):
        """only one product allowed at this point"""
        # [rsmi, psmi] = rxnsmi.split('>>')
        rsmi, resmi, psmi = rxnsmi.split('>')
        rsmi_canonical = self.remove_am(rsmi)
        psmi_canonical = self.remove_am(psmi)
        rsmi_am, outcomes = self.get_wln_resutls(rsmi_canonical)

        psmi_am = ''
        for outcome in outcomes[0]:
            # print(outcome['outcome']['smiles'])
            # print(self.remove_am(outcome['outcome']['smiles']))
            if outcome['outcome']['smiles'] == psmi_canonical:
                psmi_am = outcome['outcome']['smiles_w_am']

        if not psmi_am:
            print('Failed to find the atom mapping.')
            return '>>'
        else:
            return rsmi_am[0]+'>' + resmi + '>'+psmi_am

    def remove_am(self, smi_am):
        mol = Chem.MolFromSmiles(smi_am)
        for atom in mol.GetAtoms():
            # get rid of clearing atom mapping function, so the output will have atom maps
            atom.ClearProp('molAtomMapNumber')
        # since wln doesn't have stereochemistry, we remove stereochemistry here
        rdmolops.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol)

if __name__ == '__main__':
    mapper = WLN_AtomMapper()

    rsmi = '[C:8]([O:9][C:10](=[O:11])[N:15]1[CH2:16][CH2:17][CH:18]([CH2:21][O:22][C:23](=[O:24])[CH:25]2[N:26]3[C:27](=[O:38])[N:28]([O:33][S:34](=[O:35])(=[O:36])[OH:37])[CH:29]([CH2:30][CH2:31]2)[CH2:32]3)[CH2:19][CH2:20]1)([CH3:12])([CH3:13])[CH3:14].[F:1][C:2]([F:3])([F:4])[C:5]([OH:6])=[O:7]'
    psmi = '[NH:15]1[CH2:16][CH2:17][CH:18]([CH2:21][O:22][C:23](=[O:24])[CH:25]2[N:26]3[C:27](=[O:38])[N:28]([O:33][S:34](=[O:35])(=[O:36])[OH:37])[CH:29]([CH2:30][CH2:31]2)[CH2:32]3)[CH2:19][CH2:20]1 10-15'

    rsmi = mapper.remove_am(rsmi)
    psmi = mapper.remove_am(psmi)
    rxnsmi = rsmi+'>>'+psmi
    # rxnsmi = 'CC(C)(C)OC(=O)N1CCC(N)CC1.CC(C)(C)OC(=O)N1CCC(NC(=O)c2cccc3oc(-c4ccccc4)nc23)CC1>>CC(C)(C)OC(=O)NC1CCN(C(=O)OC(C)(C)C)CC1'
    rxnsmi_am = mapper.evaluate(rxnsmi)
    # all_mapped_rct_smiles, all_outcomes = mapper.get_wln_resutls(rsmi)