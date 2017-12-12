import os
from pymongo import MongoClient
import rdkit.Chem as Chem
import cPickle as pickle 
from rdchiral.main import rdchiralRun, rdchiralReaction, rdchiralReactants
from time import time
from makeit.webapp.transformer_v3 import Transformer, apply_one_retrotemplate

class FocusedTransformer(Transformer):
    '''
    Uses a neural network model to find the template numbers that are predicted 
    to be the most relevant for a retrosynthetic expansion
    '''


    def load_prioritizer(self, path):
        from makeit.webapp.retrotemp.retrotemp.standalone_retrotemp_numpy import RetroTempPrioritizer
        self.prioritizer = RetroTempPrioritizer()    
        self.prioritizer.restore(weight_path=path)
        smis = ['CCCOCCC', 'CCCNc1ccccc1']

        # Try a simple one
        for smi in smis:
            lst = self.prioritizer.get_topk_from_smi(smi)
            # print('{} -> {}'.format(smi, lst))

    def perform_retro(self, smiles, mincount=0, k=100, breakif=50):
        '''
        Performs a one-step retrosynthesis given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.

        mincount is unused
        '''

        # Define mol to operate on
        rct = rdchiralReactants(smiles)
        
        # Also use normal RDKit to canonicalize
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), True)

        # Initialize results object
        # result = RetroResult(smiles)

        print('Performing retro on %s' % smiles)

        # Get list of template priorities
        probs, template_nums = self.prioritizer.get_topk_from_smi(smiles, k=k)
        print('Got list of template priorities')

        # Try each in turn - keep sorted by est. applicability
        prec_list = []
        # add_precursor = result.add_precursor
        for i, template_num in enumerate(template_nums):
            for precursor in apply_one_retrotemplate(rct, smiles, self.templates[template_num]):
                # add_precursor(precursor, Pricer=self.Pricer)
                prec_list.append({
                    'template_rank': i+1,
                    'template_prob': probs[i],
                    'smiles_list': precursor.smiles_list,
                    'template_id': self.templates[template_num]['_id'],
                    'num_examples': self.templates[template_num]['count']
                })
            if len(prec_list) > breakif:
                break

        return prec_list


if __name__ == '__main__':

    db_client = MongoClient('mongodb://guest:guest@askcos2.mit.edu/admin', 27017)
    reaction_db = db_client['reaxys_v2']['reactions']

    RETRO_TRANSFORMS_CHIRAL = {
        'database': 'reaxys_v2',
        'collection': 'transforms_retro_v9',
        'mincount': 25,
        'mincount_chiral': 10
    }

    mincount_retro = RETRO_TRANSFORMS_CHIRAL['mincount']
    mincount_retro_chiral = RETRO_TRANSFORMS_CHIRAL['mincount_chiral']
    RetroFocusedTransformerChiral = FocusedTransformer()


    RetroFocusedTransformerChiral.load_prioritizer(os.path.join(
        os.path.dirname(__file__), 'retrotemp', 'models', '5d4M_Reaxys', 'model.ckpt-105660.as_numpy.pickle'))
    print('Loaded prioritizer')

    print('Loading templates now...')
    database = db_client[RETRO_TRANSFORMS_CHIRAL['database']]
    RETRO_DB = database[RETRO_TRANSFORMS_CHIRAL['collection']]
    RetroFocusedTransformerChiral.load(RETRO_DB, mincount=mincount_retro, get_retro=True, 
        get_synth=False, refs=False, mincount_chiral=mincount_retro_chiral)
    print('Loaded {} template'.format(len(RetroFocusedTransformerChiral.templates)))

    prompt = ''
    while prompt != 'done':
        prompt = raw_input('Enter SMILES: ')
        if prompt == 'done': break
        try:
            start = time()
            prec_list = RetroFocusedTransformerChiral.perform_retro(prompt.strip(), k=100, breakif=100)
            print('RESULTS CALCULATED IN {:.3f} SECONDS'.format(time()-start))
            for prec in prec_list[:25]:
                print(prec)
        except Exception as e:
            print(e)