import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
import global_config as gc
from pymongo import MongoClient
from i_o.logging import MyLogger

contexts_loc = 'contexts'

def clean_contexts(contexts):
    contexts_predictor = []
    for context in contexts:
        contexts_predictor.append(clean_context(context))
        
    return contexts_predictor

def clean_context(context):
    (T1, slvt1, rgt1, cat1, t1, y1) = context
    slvt1 = trim_trailing_period(slvt1)
    rgt1 = trim_trailing_period(rgt1)
    cat1 = trim_trailing_period(cat1)
    (rgt1, cat1, slvt1) = fix_rgt_cat_slvt(rgt1, cat1, slvt1)
    context_predictor = (T1, rgt1, slvt1)   
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

def context_to_edit(context):
    db_client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    db = db_client[gc.SOLVENTS['database']]
    SOLVENT_DB = db[gc.SOLVENTS['collection']]
    (T, reagents, solvent) = context
    # Temperature is easy
    try:
        T = float(T)
    except TypeError:
        MyLogger.print_and_log('Cannot convert temp {} to float'.format(T), contexts_loc, level = 2)
    
    # Solvent needs a lookup
    if solvent:
        try:
            solvent_mol = Chem.MolFromSmiles(solvent)
            try:
                doc = SOLVENT_DB.find_one({'_id': Chem.MolToSmiles(solvent_mol)})
                solvent = doc['name']
            except Exception as e:
                MyLogger.print_and_log('Could not parse solvent {}: {}'.format(solvent, e), contexts_loc, level = 2)
        except Exception as e:
            try:
                doc = SOLVENT_DB.find_one({'name': solvent})
                solvent = doc['name']
            except Exception as e:
                MyLogger.print_and_log('Could not parse solvent {}: {}'.format(solvent, e), contexts_loc, level = 2)
    else:
        try:
            doc = SOLVENT_DB.find_one({'_id': 'default'})
            solvent = doc['_id']
        except Exception as e:
            MyLogger.print_and_log('Could not parse solvent {}: {}'.format(solvent, e), contexts_loc, level = 2)
    
    solvent_vec = [doc['c'], doc['e'], doc['s'], doc['a'], doc['b'], doc['v']]
    
    # Unreacting reagents
    reagents_mols = [Chem.MolFromSmiles(reagent) for reagent in reagents.split('.')]
    if None in reagents_mols:
        #return 'Could not parse all reagents!'
        MyLogger.print_and_log('Could not parse all reagents', contexts_loc, level = 2)

    reagent_fp = np.zeros((1, 256))
    
    for reagent in reagents_mols:
        reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
    # Save list
    return [reagent_fp, np.reshape(np.array(solvent_vec), (1, 6)), np.reshape(np.array(T), (1, 1))]
    '''
    reagent_fp = [0 for i in range (256)]
    for reagent in reagents_mols:
        #reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
        reagent_fp = [sum(x) for x in zip(reagent_fp, AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))]
    return [np.array(reagent_fp), np.array(solvent_vec), np.array([[T]])]
    '''

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
            errors.append('cannot convert temp {} to float'.format(T))
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
            errors.append('could not parse all reagents')
            continue
        reagent_fp = np.zeros((1, 256))
        
        for reagent in reagents_mols:
            reagent_fp += np.array(AllChem.GetMorganFingerprintAsBitVect(reagent, 2, nBits = 256))
    
        # Save list
        contexts_edits.append([reagent_fp, np.reshape(np.array(solvent_vec), (1, 6)), np.reshape(np.array(T), (1, 1))])
        context_labels.append('T:{}, rgt:{}, solv:{}'.format(T, reagents, solvent))
        errors.append(False)
        
    return contexts_edits