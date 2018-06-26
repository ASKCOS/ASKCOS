import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
import makeit.global_config as gc
from pymongo import MongoClient
from makeit.utilities.io.logging import MyLogger

contexts_loc = 'contexts'

def clean_contexts(contexts):
    contexts_predictor = []
    for context in contexts:
        contexts_predictor.append(clean_context(context))
        
    return contexts_predictor

def clean_context(context):
    (T1, slvt1, rgt1, cat1, t1, y1) = context
    ##remove chemicals without parsible smiles
    slvs = slvt1.split('.')
    rgts = rgt1.split('.')
    cats = cat1.split('.')
    rgt_name = [rgt for rgt in rgts if 'Reaxys' not in rgt]
    slv_name = [slv for slv in slvs if 'Reaxys' not in slv]
    cat_name = [cat for cat in cats if 'Reaxys' not in cat]
    slvt1 = '.'.join(slvs)
    rgt1 = '.'.join(rgts)
    cat1 = '.'.join(cats)

    slvt1 = trim_trailing_period(slvt1)
    rgt1 = trim_trailing_period(rgt1)
    cat1 = trim_trailing_period(cat1)
    (rgt1, cat1, slvt1) = fix_rgt_cat_slvt(rgt1, cat1, slvt1)
    context_predictor = (T1, slvt1, rgt1, cat1, t1, y1)   
    return context_predictor

def fix_rgt_cat_slvt(rgt1, cat1, slvt1):
    # Merge cat and reagent for forward predictor
    if rgt1 and cat1:
        rgt1 = rgt1 + '.' + cat1
    elif cat1:
        rgt1 = cat1
    # Reduce solvent to single one
    if '.' in slvt1:
        slvt1 = slvt1.split('.')[0]
    return (rgt1, cat1, slvt1)

def trim_trailing_period(txt):
    if txt:
        if txt[-1] == '.':
            return txt[:-1]
    return txt

def context_to_edit(context, solvent_name_to_smiles, solvent_smiles_to_params):
    
    (T, solvent, reagents_str, cat1, t1, y1) = context
    # Temperature is easy
    try:
        T = float(T)
    except TypeError:
        MyLogger.print_and_log('Cannot convert temperature {} to float'.format(T), contexts_loc)
        return None

    # Look up solvent from saved dicts
    try:
        if solvent == '':
            solvent = 'default'
            
        solvent_mol = Chem.MolFromSmiles(solvent)
        doc = solvent_smiles_to_params[Chem.MolToSmiles(solvent_mol)]
    except Exception as e: # KeyError or Boost.Python.ArgumentError
        try:
            doc = solvent_smiles_to_params[solvent_name_to_smiles[solvent]]
        except KeyError:
            MyLogger.print_and_log('Could not parse solvent {}'.format(solvent), contexts_loc)
            doc = solvent_smiles_to_params[solvent_name_to_smiles['default']]
    solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]

    # Unreacting reagents
    reagents = [Chem.MolFromSmiles(reagent) for reagent in reagents_str.split('.')]
    if None in reagents:
        MyLogger.print_and_log('Could not parse all reagents!{}'.format(reagents_str), contexts_loc)
        return None
    reagent_fp = np.zeros((1, 256))
    for reagent in reagents:
        reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
    
    # Return list
    return [reagent_fp, np.reshape(np.array(solvent_vec), (1, 6)), np.reshape(np.array(T), (1, 1))]


def contexts_to_edits(contexts):
    '''Set multiple contexts to try at once'''
    db_client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    db = db_client[gc.SOLVENTS['database']]
    SOLVENT_DB = db[gc.SOLVENTS['collection']]
    contexts_edits = []
    context_labels = []
    errors = []
    for (T, reagents, solvent) in contexts:
        # Temperature is easy
        try:
            T = float(T)
        except TypeError:
            #return 'Cannot convert temperature {} to float'.format(T)
            MyLogger.print_and_log('Cannot convert temp {} to float'.format(T), contexts_loc, level= 2)
            continue
        # Solvent needs a lookup
        if solvent:
            try:
                solvent_mol = Chem.MolFromSmiles(solvent)
                try:
                    doc = SOLVENT_DB.find_one({'_id': Chem.MolToSmiles(solvent_mol)})
                    solvent = doc['name']
                except Exception as e:
                    MyLogger.print_and_log('Could not parse solvent {}'.format(solvent), contexts_loc, level = 2)
            except Exception as e:
                try:
                    doc = SOLVENT_DB.find_one({'name': solvent})
                    solvent = doc['name']
                except Exception as e:
                    MyLogger.print_and_log('Could not parse solvent {}'.format(solvent), contexts_loc, level = 2)
        else:
            try:
                doc = SOLVENT_DB.find_one({'_id': 'default'})
                solvent = doc['_id']
            except Exception as e:
                MyLogger.print_and_log('Could not parse solvent {}'.format(solvent), contexts_loc, level = 2)
        
        solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
        
        # Unreacting reagents
       
        reagents_mols = [Chem.MolFromSmiles(reagent) for reagent in reagents.split('.')]
        if None in reagents_mols:
            #return 'Could not parse all reagents!'
            MyLogger.print_and_log('Could not parse all reagents', contexts_loc)
            continue
        reagent_fp = np.zeros((1, 256))
        
        for reagent in reagents_mols:
            reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
    
        # Save list
        contexts_edits.append([reagent_fp, np.reshape(np.array(solvent_vec), (1, 6)), np.reshape(np.array(T), (1, 1))])
        context_labels.append('T:{}, rgt:{}, solv:{}'.format(T, reagents, solvent))
        errors.append(False)
        
    return contexts_edits
