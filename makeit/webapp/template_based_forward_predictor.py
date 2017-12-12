from __future__ import absolute_import, unicode_literals, print_function

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

import rdkit.Chem as Chem 
import rdkit.Chem.AllChem as AllChem
from makeit.embedding.descriptors import edits_to_vectors
from makeit.predict.summarize_reaction_outcome import summarize_reaction_outcome

from collections import defaultdict
import numpy as np 
import os 

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()

class TemplateBasedForwardPredictor():
    '''A class to organize everything needed for template-based forward 
    prediction (slow). This method of evaluation should be phased out in 
    favor of more scalable and accurate methods.'''

    def __init__(self):
        self.templates = []
        self.solvent_name_to_smiles = {}
        self.solvent_smiles_to_params = {}
        self.F_atom = None 
        self.F_bond = None

    def load_templates(self, TEMPLATE_DB, mincount=25, countsonly=False):
        '''
        Load synthetic templates for forward prediction. Does not use a separate
        forward enumeration class.
        '''
        print('Loading synthetic templates for forward predictor')
        self.templates = []
        for doc in TEMPLATE_DB.find({'count': {'$gte': mincount}}, ['_id', 'count', 'reaction_smarts']):
            if 'reaction_smarts' not in doc: continue
            reaction_smarts = doc['reaction_smarts']
            if not reaction_smarts: continue
            template = {
                'reaction_smarts':      reaction_smarts,
                'count':                doc['count'] if 'count' in doc else 0,
                '_id':                  doc['_id'] if '_id' in doc else -1,
            }
            try:
                rxn_f = AllChem.ReactionFromSmarts(str('(' + reaction_smarts.replace('>>', ')>>(') + ')'))
                if rxn_f.Validate()[1] != 0:
                    print('Could not validate {}'.format(reaction_smarts))
                    continue
                template['rxn_f'] = rxn_f
            except Exception as e:
                print('Couldnt load forward: {}: {}'.format(reaction_smarts, e))
                continue
            self.templates.append(template)

        self.num_templates = len(self.templates)
        self.templates = sorted(self.templates, key=lambda z: z['count'], reverse=True)
        print('Loaded {} templates'.format(self.num_templates))
        self.template_counts = [x['count'] for x in self.templates]
        
        # Only need template counts? Free up memory by deleting templates
        if countsonly:
            del self.templates

    def load_solvents(self, SOLVENT_DB):
        '''Get Abraham parameters from solvent DB'''
        for doc in SOLVENT_DB.find():
            try:
                self.solvent_name_to_smiles[doc['name']] = doc['_id']
                self.solvent_smiles_to_params[doc['_id']] = doc 
            except KeyError:
                print('solvent doc {} missing a name'.format(doc))

    def load_model(self, folder):
        '''Load the trained Keras model'''
        # Feature sizes
        print('Getting feature sizes')
        mol = Chem.MolFromSmiles('[C:1][C:2]')
        (a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)
        self.F_atom = len(a[0])
        self.F_bond = len(b[0])

            # Get model args
        ARGS_FPATH = os.path.join(folder, 'args.json')
        with open(ARGS_FPATH, 'r') as fid:
            import json
            args = json.load(fid)

        N_h2 = int(args['Nh2'])
        N_h1 = int(args['Nh1'])
        N_h3 = int(args['Nh3'])
        N_hf = int(args['Nhf'])
        l2v = float(args['l2'])
        lr = float(args['lr'])
        context_weight = float(args['context_weight'])
        enhancement_weight = float(args['enhancement_weight'])
        optimizer          = args['optimizer']
        inner_act          = args['inner_act']
        TARGET_YIELD       = False

        print('building model')
        from makeit.predict.reaxys_score_candidates_from_edits_compact_v2 import build
        self.model = build(F_atom=self.F_atom, F_bond=self.F_bond, N_h1=N_h1, 
                N_h2=N_h2, N_h3=N_h3, N_hf=N_hf, l2v=l2v, inner_act=inner_act,
                context_weight=context_weight, enhancement_weight=enhancement_weight, TARGET_YIELD=TARGET_YIELD,
                absolute_score=True)

        print('Loading weights')
        WEIGHTS_FPATH = os.path.join(folder, 'weights.h5')
        self.model.load_weights(WEIGHTS_FPATH, by_name=True)

    def context_to_mats(self, reagents='', solvent='toluene', T='20'):
        '''Convert a context to a matrix as input to the trained model'''
        if not self.solvent_smiles_to_params:
            raise ValueError('Need to load solvents first!')

        # Temperature is easy
        try:
            T = float(T)
        except TypeError:
            print('Cannot convert temperature {} to float'.format(T))
            return None

        # Look up solvent from saved dicts
        try:
            solvent_mol = Chem.MolFromSmiles(solvent)
            doc = self.solvent_smiles_to_params[Chem.MolToSmiles(solvent_mol)]
        except Exception as e: # KeyError or Boost.Python.ArgumentError
            try:
                doc = self.solvent_smiles_to_params[solvent_name_to_smiles[solvent]]
            except KeyError:
                print('Could not parse solvent {}'.format(solvent))
                return None
        solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
        solvent = doc['name']

        # Unreacting reagents
        reagents = [Chem.MolFromSmiles(reagent) for reagent in reagents.split('.')]
        if None in reagents:
            print('Could not parse all reagents!')
            return None
        reagent_fp = np.zeros((1, 256))
        for reagent in reagents:
            reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
        
        # Return list
        return [reagent_fp, np.reshape(np.array(solvent_vec), (1, 6)), np.reshape(np.array(T), (1, 1))]

    def prepare_reactants_smiles(self, reactant_smiles):
        '''Parse, set mapping, and return SMILES'''

        reactants = Chem.MolFromSmiles(reactant_smiles)
        if not reactants: 
            print('Could not parse reactants {}'.format(reactant_smiles))
            return [[] for i in range(len(contexts))]
        print('Number of reactant atoms: {}'.format(len(reactants.GetAtoms())))
        # Report current reactant SMILES string
        [a.ClearProp(str('molAtomMapNumber')) for a in reactants.GetAtoms() if a.HasProp(str('molAtomMapNumber'))]
        reactants_smiles_no_map = Chem.MolToSmiles(reactants, isomericSmiles=True)
        print('Reactants w/o map: {}'.format(reactants_smiles_no_map))
        # Add new atom map numbers
        [a.SetProp(str('molAtomMapNumber'), str(i+1)) for (i, a) in enumerate(reactants.GetAtoms())]
        # Report new reactant SMILES string
        print('Reactants w/ map: {}'.format(Chem.MolToSmiles(reactants, isomericSmiles=True)))
        return Chem.MolToSmiles(reactants)

    def get_candidate_edits_range(self, reactant_smiles, start_at=0, end_at=1e9):
        '''Get candidate edits for templates with indices between start_at
        and end_at

        reactants_smiles should be atom-mapped'''

        reactants = Chem.MolFromSmiles(reactant_smiles)
        candidate_list = []
        for i in range(start_at, end_at):
            try:

                outcomes = self.templates[i]['rxn_f'].RunReactants([reactants])
                if not outcomes: continue # no match

                for j, outcome in enumerate(outcomes):
                    outcome = outcome[0] # all products represented as single mol by transforms

                    try:
                        outcome.UpdatePropertyCache()
                        Chem.SanitizeMol(outcome)
                        [a.SetProp(str('molAtomMapNumber'), a.GetProp(str('old_molAtomMapNumber'))) \
                            for a in outcome.GetAtoms() \
                            if str('old_molAtomMapNumber') in a.GetPropsAsDict()]
                    
                        # Reduce to largest (longest) product only
                        candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles=True)
                        candidate_smiles = max(candidate_smiles.split('.'), key=len)
                        outcome = Chem.MolFromSmiles(candidate_smiles)
                            
                        # Find what edits were made
                        edits = summarize_reaction_outcome(reactants, outcome)

                        # Remove mapping before matching
                        [x.ClearProp(str('molAtomMapNumber')) for x in outcome.GetAtoms() \
                            if x.HasProp(str('molAtomMapNumber'))] # remove atom mapping from outcome

                        # Overwrite candidate_smiles without atom mapping numbers
                        candidate_smiles = Chem.MolToSmiles(outcome, isomericSmiles=True)

                        # Add to ongoing list
                        candidate_list.append((candidate_smiles, edits))
                    except Exception as e: # other RDKit error?
                        #print(e) # fail quietly
                        continue
            except IndexError: # out of range w/ templates
                print('INDEX ERROR!')
                break 
            except Exception as e: # other RDKit error?
                print(e)
                continue    
        return candidate_list

    def get_template_range(self, reactants_smiles, mincount=25):
        '''Get candidate edits for templates with a certain mincount'''
        end_at = len(self.template_counts)
        for (i, count) in enumerate(self.template_counts):
            if count < mincount:
                end_at = i 
                break
        print('Range to apply is template index {} through {}'.format(0, end_at))
        return (0, end_at)

    def score_candidate_edits(self, reactant_smiles, contexts, candidate_edits, top_n=50):
        '''Get and score outcomes'''
        
        # Pre-calc descriptors for this set of reactants
        reactants = Chem.MolFromSmiles(reactant_smiles)
        atom_desc_dict = edits_to_vectors([], reactants, return_atom_desc_dict=True)

        # Convert to tensors
        Nc = len(candidate_edits)
        if Nc == 0:
            print('No candidate products found at all...?')
            return [[] for i in range(len(contexts))]
        Ne1 = 1
        Ne2 = 1
        Ne3 = 1
        Ne4 = 1
        for (candidate_smiles, edits) in candidate_edits:
            Ne1 = max(Ne1, len(edits[0]))
            Ne2 = max(Ne2, len(edits[1]))
            Ne3 = max(Ne3, len(edits[2]))
            Ne4 = max(Ne4, len(edits[3]))
        x_h_lost = np.zeros((1, Nc, Ne1, self.F_atom))
        x_h_gain = np.zeros((1, Nc, Ne2, self.F_atom))
        x_bond_lost = np.zeros((1, Nc, Ne3, self.F_bond))
        x_bond_gain = np.zeros((1, Nc, Ne4, self.F_bond))

        smiles = []
        for c, (candidate_smiles, edits) in enumerate(candidate_edits):
            smiles.append(candidate_smiles)
            edit_h_lost_vec, edit_h_gain_vec, \
                edit_bond_lost_vec, edit_bond_gain_vec = edits_to_vectors(edits, reactants, 
                    atom_desc_dict=atom_desc_dict)

            for (e, edit_h_lost) in enumerate(edit_h_lost_vec):
                x_h_lost[0, c, e, :] = edit_h_lost
            for (e, edit_h_gain) in enumerate(edit_h_gain_vec):
                x_h_gain[0, c, e, :] = edit_h_gain
            for (e, edit_bond_lost) in enumerate(edit_bond_lost_vec):
                x_bond_lost[0, c, e, :] = edit_bond_lost
            for (e, edit_bond_gain) in enumerate(edit_bond_gain_vec):
                x_bond_gain[0, c, e, :] = edit_bond_gain

        # Get rid of NaNs/Infs
        x_h_lost[np.isnan(x_h_lost)] = 0.0
        x_h_gain[np.isnan(x_h_gain)] = 0.0
        x_bond_lost[np.isnan(x_bond_lost)] = 0.0
        x_bond_gain[np.isnan(x_bond_gain)] = 0.0
        x_h_lost[np.isinf(x_h_lost)] = 0.0
        x_h_gain[np.isinf(x_h_gain)] = 0.0
        x_bond_lost[np.isinf(x_bond_lost)] = 0.0
        x_bond_gain[np.isinf(x_bond_gain)] = 0.0
        x = [x_h_lost, x_h_gain, x_bond_lost, x_bond_gain]

        # Make prediction(s)
        all_outcomes = []
        for (T, reagents, solvent) in contexts:
            xc = self.context_to_mats(reagents=reagents, solvent=solvent, T=T)
            if xc is None: # unparseable
                all_outcomes.append([])
                continue
            scores = self.model.predict(x + xc)[0]
            probs = softmax(scores)

            this_outcome = defaultdict(float)
            this_scores  = defaultdict(lambda: -999)
            for (smi, sco, pro) in zip(smiles, scores, probs):
                this_outcome[smi] += pro 
                this_scores[smi] = max(this_scores[smi], sco)

            # Sort
            this_outcome = sorted(this_outcome.iteritems(), key=lambda x: x[1], reverse=True)
            outcomes = []
            for i, (smi, pro) in enumerate(this_outcome):
                if i == top_n: break
                outcomes.append({
                    'rank': i + 1,
                    'smiles': smi,
                    'score': '{:.2f}'.format(this_scores[smi]),
                    'prob': '{:.4f}'.format(pro*100.),
                })
            all_outcomes.append(outcomes)

        return all_outcomes

if __name__ == '__main__':
    from pymongo import MongoClient
    db_client = MongoClient('mongodb://guest:guest@askcos2.mit.edu/admin', 27017)
    TEMPLATE_DB = db_client['reaxys_v2']['transforms_forward_v2']
    SOLVENT_DB = db_client['reaxys']['solvents']
    MODEL_FOLDER = '/home/ccoley/Make-It/makeit/predict/output/01_23_2017'

    predictor = TemplateBasedForwardPredictor()
    predictor.load_templates(TEMPLATE_DB)
    predictor.load_solvents(SOLVENT_DB)
    predictor.load_model(MODEL_FOLDER)

    reactant_smiles = 'CCCCCC(=O)O.NCCCCCC'
    contexts = [('20', '[Na+].[OH-]', 'O'), ('100', '', 'c1ccccc1')]

    reactant_smiles = predictor.prepare_reactants_smiles(reactant_smiles)
    (start_at, end_at) = predictor.get_template_range(reactant_smiles)
    candidate_edits = predictor.get_candidate_edits_range(reactant_smiles, start_at, end_at)
    outcomes = predictor.score_candidate_edits(reactant_smiles, contexts, candidate_edits)
    print(outcomes)