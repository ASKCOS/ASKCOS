import global_config as gc
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from i_o.logging import MyLogger
from parsing import check_smiles
fingerprinting_loc = 'fingerprinting'

def create_rxn_Morgan2FP(rxn_smiles, fpsize=gc.fingerprint_bits, useFeatures=True):
    """Create a rxn Morgan (r=2) fingerprint as bit vector from a reaction SMILES string

        Modified from Schneider's code (2014)"""

    rsmi = rxn_smiles.split('>')[0].split('.')
    psmi = rxn_smiles.split('>')[2].split('.')

    rfp = None
    pfp = None
    for react in rsmi:
        mol = Chem.MolFromSmiles(react)
        fp = np.array(
            AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=fpsize, useFeatures=useFeatures))
        if rfp is None:
            rfp = fp
        else:
            rfp += fp
    for product in psmi:
        mol = Chem.MolFromSmiles(product)
        fp = np.array(
            AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=fpsize, useFeatures=useFeatures))
        if pfp is None:
            pfp = fp
        else:
            pfp += fp

    if pfp is not None and rfp is not None:
        pfp -= rfp
    return pfp

def get_condition_input_from_smiles(conditions_smiles , split = False):
    '''
    If split is used: first molecule in the conditions_smiles should be the solvent!
    '''
    conditions_mol = []
    inputlength = 0
    #If input is a single string: immediately extract molecule
    if (type(conditions_smiles) is str) or (type(conditions_smiles) is unicode):
        try:
            conditions_mol.append(Chem.MolFromSmiles(conditions_smiles))
            inputlength += 1
        #If any unparsable: skip this iteration
        except Exception as e:
            MyLogger.print_and_log('Unparsable conditions. Leaving out this reaction: {}.'.format(instance['_id']), fingerprinting_loc, level = 1)
    #otherwise: assume list of strings
    else:
        for smiles in conditions_smiles: 
            try:
                conditions_mol.append(Chem.MolFromSmiles(smiles))
                inputlength += 1
            #If any unparsable: skip this iteration
            except Exception as e:
                MyLogger.print_and_log('Unparsable conditions. Leaving out this reaction: {}.'.format(instance['_id']), fingerprinting_loc, level = 1)
                continue
        
    fps = []
    for mol in conditions_mol:
        try:
            p = AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=gc.fingerprint_bits, useFeatures=True)
        except Exception as e:
            p = [0]*gc.fingerprint_bits
        fps.append(p)
    if split:
        for fp in fps:
            fp = np.array(fp)
        return fps
    else:
        fps = np.array(fp).reshape(1,gc.fingerprint_bits*inputlength)
        input = fps
        return input
def get_condition_input_from_instance(instance, chemicals, asone = False, astwo = False, use_new = False, split = False):
    conditions_smiles = get_input_condition_as_smiles(instance, chemicals, asone = asone, astwo = astwo, use_new = use_new)
    return get_condition_input_from_smiles(conditions_smiles, split = split)

def get_reaction_input_from_instance(instance, reactions, chemicals):
    
    reaction_smiles = get_reaction_as_smiles(instance, reactions, chemicals)
        
    return get_reaction_input_from_smiles(reaction_smiles)

def get_reaction_input_from_smiles(reaction_smiles):

    react = reaction_smiles.split('>')
    reag_s = react[0]
    prod_s = react[len(react)-1]
    reag = Chem.MolFromSmiles(reag_s)
    prod = Chem.MolFromSmiles(prod_s)
    
    #factor for reactionfp length
    reactionfplength = 4
    reag = AllChem.GetMorganFingerprintAsBitVect(mol = reag, radius = 2, nBits = gc.fingerprint_bits*reactionfplength)
    prod = AllChem.GetMorganFingerprintAsBitVect(mol = prod, radius = 2, nBits = gc.fingerprint_bits*reactionfplength)
    reactionfp = [i - j for i, j in zip(prod, reag)]
    reaction = []
    for i in range(gc.fingerprint_bits):
        pos = 0
        for j in range(reactionfplength):
            pos += reactionfp[i+gc.fingerprint_bits*j]
        reaction.append(pos)
    input = np.array(reaction).reshape(1, gc.fingerprint_bits)  
    return input

def get_reaction_as_smiles(instance, reactions, chemicals):
    id = instance['RX_ID'][0]
    react = reactions.find_one({'_id':id})
    reaction_smiles = ''
    reactants = react['RX_RXRN']
    products = react['RX_PXRN']
    for reactant in reactants:
        #Try using newest smiles.
        try:
            chemical = chemicals.find_one({'_id':reactant})
            smiles = chemical['SMILES_new']
            mol = Chem.MolFromSmiles(smiles)
            if(smiles and (not mol)):
                new_smiles = raw_input('New smiles for reactant: {}:'.format(smiles))
                if new_smiles:
                    chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                else:
                    pass
        except TypeError:
            MyLogger.print_and_log('No reactant smiles found for reaction {} returning "NONE"'.format(id), fingerprinting_loc, 2)
            return 'NONE'
        except KeyError:
            #If not present: used older ones
            try:
                chemical = chemicals.find_one({'_id':reactant})
                smiles = chemicals.find_one({'_id':reactant})['SMILES']
                mol = Chem.MolFromSmiles(smiles)
                if(smiles and (not mol)):
                    new_smiles = raw_input('New smiles for reactant: {}:'.format(smiles))
                    if new_smiles:
                        chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                    else:
                        pass
                    
            except TypeError:
                MyLogger.print_and_log('No reactant smiles found for reaction {} returning "NONE"'.format(id), fingerprinting_loc, 2)
                return 'NONE'
        if not reaction_smiles:
            reaction_smiles = smiles
        else:
            reaction_smiles += '.'+smiles
    reaction_smiles += '>>'
    for product in products:
        try:
            chemical = chemicals.find_one({'_id':product})
            smiles = chemical['SMILES_new']
            '''
            mol = Chem.MolFromSmiles(smiles)
            if(smiles and (not mol)):
                new_smiles = raw_input('New smiles for product: {}'.format(smiles))
                if new_smiles:
                    chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                else:
                    pass
            '''
        except TypeError:
            MyLogger.print_and_log('No product smiles found for reaction {} returning "NONE"'.format(id), fingerprinting_loc, 2)
            return 'NONE'
        except KeyError:
            try:
                chemical = chemicals.find_one({'_id':product})
                smiles = chemicals.find_one({'_id':product})['SMILES']
                '''
                mol = Chem.MolFromSmiles(smiles)
                if(smiles and (not mol)):
                    new_smiles = raw_input('New smiles for product: {}'.format(smiles))
                    if new_smiles:
                        chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                    else:
                        pass
                '''
            except TypeError:
                MyLogger.print_and_log('No product smiles found for reaction {} returning "NONE"'.format(id), fingerprinting_loc, 2)
                return 'NONE'
        if reaction_smiles[len(reaction_smiles)-1] == '>':
            reaction_smiles += smiles
        else:
            reaction_smiles += '.'+smiles
    return reaction_smiles
                
def get_input_condition_as_smiles(doc,chemicals, asone = False, astwo = False, check = False, use_new = False):
    a = doc['RXD_SOLXRN']
    b = doc['RXD_RGTXRN']
    c = doc['RXD_CATXRN']
    solv = ""
    reag = ""
    cata = ""
    for sol in a:
        sa = chemicals.find_one({'_id':sol})
        if sa:
            if check:
                check_smiles(sa,chemicals)
                s = None
                try:
                    s = sa['SMILES_new']
                except:
                    s= sa['SMILES']
            elif use_new:
                try:
                    s = sa['SMILES_new']
                except:
                    s= sa['SMILES']
            else:
                s = sa['SMILES']
            '''
            mol = Chem.MolFromSmiles(s)
            if not mol:
                new_smiles = raw_input('New smiles for solvent: {},{}'.format(sa['_id'],s))
                if new_smiles:
                    chemicals.update({'_id':sa['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                    s=new_smiles
                else:
                    pass
            '''
        else:
            continue
        if solv:
            solv += '.'+s
        elif s:
            solv += s
    for rgt in b:
        sa = chemicals.find_one({'_id':rgt})
        if sa:
            if check:
                check_smiles(sa,chemicals)
                s = None
                try:
                    s = sa['SMILES_new']
                except:
                    s= sa['SMILES']
            elif use_new:
                try:
                    s = sa['SMILES_new']
                except:
                    s= sa['SMILES']
            else:
                s = sa['SMILES']
            '''
            mol = Chem.MolFromSmiles(s)
            if not mol:
                new_smiles = raw_input('New smiles for reagent: {},{}'.format(sa['_id'],s))
                if new_smiles:
                    chemicals.update({'_id':sa['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                    s=new_smiles
                else:
                    pass
            '''
        else:
            continue
        if(reag and s):
            reag += '.'+s
        elif s:
            reag += s
    for cat in c:
        sa = chemicals.find_one({'_id':cat})
        if sa:
            if check:
                check_smiles(sa,chemicals)
                s = None
                try:
                    s = sa['SMILES_new']
                except:
                    s= sa['SMILES']
            elif use_new:
                try:
                    s = sa['SMILES_new']
                except:
                    s= sa['SMILES']
            else:
                s = sa['SMILES']
            '''
            mol = Chem.MolFromSmiles(s)
            if not mol:
                new_smiles = raw_input('New smiles for catalyst: {},{}'.format(sa['_id'],s))
                if new_smiles:
                    chemicals.update({'_id':sa['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                    s=new_smiles
                else:
                    pass
            '''
        else:
            continue
        if(cata and s):
            cata += '.'+s
        elif s:
            cata += s
    
    if asone:
        ans = ""
        if solv:
            ans = solv
        if reag:
            if ans:
                ans += '.'+reag
            else:
                ans = reag
        if cata:
            if ans:
                ans += '.'+cata
            else:
                ans = cata
        return ans
    if astwo:
        ans = ""
        if reag:
            ans = reag
        if cata:
            if ans:
                ans += '.'+cata
            else:
                ans = cata
        return [solv, ans]
    else:
        return [solv, reag, cata]