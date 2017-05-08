'''
The role of a forward predictor coordinator is to generate a plausibility
score for an intended reaction. It offloads the template-expansion
work onto a separate worker pool. Returned candidate edits are used to
build up the tensors used with the trained Keras model for forward
prediction.
'''

from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from celery.result import allow_join_result 
# NOTE: allow_join_result is only because the worker is separate
import numpy as np 
import time

CORRESPONDING_QUEUE = 'fp_coordinator'
SOLVENT_DB = None
model = None
template_counts = None
F_atom = None
F_bond = None

def load_model(folder):
    '''Load a trained Keras model'''

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

    from makeit.predict.reaxys_score_candidates_from_edits_compact_v2 import build
    model = build(F_atom = self.F_atom, F_bond = self.F_bond, N_h1 = N_h1, 
            N_h2 = N_h2, N_h3 = N_h3, N_hf = N_hf, l2v = l2v, inner_act = inner_act,
            context_weight = context_weight, enhancement_weight = enhancement_weight, TARGET_YIELD = TARGET_YIELD,
            absolute_score = True)

    WEIGHTS_FPATH = os.path.join(folder, 'weights.h5')
    model.load_weights(WEIGHTS_FPATH, by_name = True)
    return model

def context_to_mats(SOLVENT_DB, reagents='', solvent='toluene', T='20'):
    # Temperature is easy
    try:
        T = float(T)
    except TypeError:
        return 'Cannot convert temperature {} to float'.format(T)

    # Solvent needs a lookup
    solvent_mol = Chem.MolFromSmiles(solvent)
    if solvent_mol: 
        doc = SOLVENT_DB.find_one({'_id': Chem.MolToSmiles(solvent_mol)})
    else:
        doc = SOLVENT_DB.find_one({'name': solvent})
    if not doc:
        return 'Could not parse solvent {}'.format(solvent)
    solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
    solvent = doc['name']

    # Unreacting reagents
    reagents = [Chem.MolFromSmiles(reagent) for reagent in reagents.split('.')]
    if None in reagents:
        return 'Could not parse all reagents!'
    reagent_fp = np.zeros((1, 256))
    for reagent in reagents:
        reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
    
    # Return list
    return [reagent_fp, np.reshape(np.array(solvent_vec), (1, 6)), np.reshape(np.array(T), (1, 1))]

@celeryd_init.connect
def configure_coordinator(options={},**kwargs):
    if 'queues' not in options: 
        return 
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A FORWARD PREDICTOR COORDINATOR ###')

    global SOLVENT_DB
    global template_counts
    global model 
    global F_atom 
    global F_bond
    
    # Get Django settings
    from django.conf import settings

    # Database
    from database import db_client
    db = db_client[settings.SYNTH_TRANSFORMS['database']]
    SYNTH_DB = database[settings.SYNTH_TRANSFORMS['collection']]
    db = db_client[settings.SOLVENTS['database']]
    SOLVENT_DB = db[settings.SOLVENTS['collection']]

    # RDKit for layer sizes
    import rdkit.Chem as Chem 
    import rdkit.Chem.AllChem as AllChem 

    # Misc.
    from makeit.embedding.descriptors import edits_to_vectors
    from .common import load_templates
    from .worker import get_candidate_edits

    # Feature sizes
    mol = Chem.MolFromSmiles('[C:1][C:2]')
    (a, _, b, _) = edits_to_vectors((['1'],[],[('1','2',1.0)],[]), mol)
    F_atom = len(a[0])
    F_bond = len(b[0])

    # Intelligent predictor
    mincount_synth = settings.SYNTH_TRANSFORMS['mincount']
    template_counts = load_templates(SYNTH_DB=SYNTH_DB, mincount=mincount_synth, 
        countsonly=True)
    model = load_model(settings.PREDICTOR['trained_model_path'])
    print('Loaded trained Keras model')

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()

@shared_task
def get_outcomes(reactants, contexts, mincount=0, top_n=10, chunksize=10):
    '''Evaluate the plausibility of a proposed forward reaction

    reactants = SMILES of reactants
    reagents = smiles of reagents
    solvent = smiles or name of solvent
    T = temperature
    mincount = minimum count for forward synthetic templates'''

    print('Treebuilder coordinator was asked to expand {}'.format(smiles))

    global SOLVENT_DB
    global template_counts
    global model 
    global F_atom 
    global F_bond

    # Get atom descriptors
    reactants = Chem.MolFromSmiles(reactants)
    if not reactants: 
        print('Could not parse reactants {}'.format(reactants))
        raise ValueError('Could not parse reactants')
    print('Number of reactant atoms: {}'.format(len(reactants.GetAtoms())))
    # Report current reactant SMILES string
    [a.ClearProp('molAtomMapNumber') for a in reactants.GetAtoms() if a.HasProp('molAtomMapNumber')]
    print('Reactants w/o map: {}'.format(Chem.MolToSmiles(reactants)))
    # Add new atom map numbers
    [a.SetProp('molAtomMapNumber', str(i+1)) for (i, a) in enumerate(reactants.GetAtoms())]
    # Report new reactant SMILES string
    print('Reactants w/ map: {}'.format(Chem.MolToSmiles(reactants)))
    reactants_smiles = Chem.MolToSmiles(reactants)

    if intended_product:
        intended_product = Chem.MolFromSmiles(intended_product)
    if intended_product:
        intended_product = Chem.MolToSmiles(intended_product, isomericSmiles=True)

    # Pre-calc descriptors for this set of reactants
    atom_desc_dict = edits_to_vectors([], reactants, return_atom_desc_dict=True)

    # Figure out what templates we are going to use
    end_at = 0
    for (i, count) in enumerate(template_counts):
        if count < mincount:
            end_at = i
            break

    # Chunk and add to queue
    pending_results = []
    for start_at in range(0, end_at, chunksize):
        pending_results.append(
            get_candidate_edits.delay(reactants=reactants, start_at=start_at, 
                end_at=start_at + chunksize)
        )

    # Wait to collect all candidates
    candidate_edits = set()
    while len(pending_results) > 0:

        # Look for done results
        is_ready = [i for (i, res) in enumerate(pending_results) if res.ready()]
        with allow_join_result(): # required to use .get()
            for i in is_ready:
                for candidate_edit in pending_results[i].get(timeout=1):
                    candidate_edits.add(candidate_edit)
                pending_results[i].forget()
                
        # Update list
        pending_results = [res for (i, res) in enumerate(pending_results) if i not in is_ready]


    # Convert to tensors
    Nc = len(candidate_edits)
    if Nc == 0:
        raise ValueError('No candidate products found at all...?')
    Ne1 = 1
    Ne2 = 1
    Ne3 = 1
    Ne4 = 1
    for (candidate_smiles, edits) in candidate_edits:
        Ne1 = max(Ne1, len(edits[0]))
        Ne2 = max(Ne2, len(edits[1]))
        Ne3 = max(Ne3, len(edits[2]))
        Ne4 = max(Ne4, len(edits[3]))
    x_h_lost = np.zeros((1, Nc, Ne1, F_atom))
    x_h_gain = np.zeros((1, Nc, Ne2, F_atom))
    x_bond_lost = np.zeros((1, Nc, Ne3, F_bond))
    x_bond_gain = np.zeros((1, Nc, Ne4, F_bond))

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
    all_plausibilities = []
    for (T, reagents, solvent) in contexts:
        xc = context_to_mats(SOLVENT_DB, reagents=reagents, solvent=solvent, T=T)
        scores = np.flatten(model.predict(x + context)[0])
        probs = softmax(scores)

        # Sort
        this_outcome = sorted(zip(smiles, scores, probs), key=lambda x: x[1], reverse=True)
        outcomes = []
        for i, outcome in enumerate(this_outcome):
            if i == top_n: break
            outcomes.append({
                'rank': i + 1,
                'smiles': outcome[0],
                'score': '{:.2f}'.format(outcome[1]),
                'prob': '{:.2e}'.format(outcome[2]),
            })
        all_outcomes.append(outcomes)

    return all_outcomes
