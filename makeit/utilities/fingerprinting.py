import makeit.global_config as gc
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.parsing import check_smiles
fingerprinting_loc = 'fingerprinting'


def create_rxn_Morgan2FP(rxn_smiles, fpsize=gc.fingerprint_bits, useFeatures=True, useChirality=False):
    """Create a rxn Morgan (r=2) fingerprint as bit vector from a reaction SMILES string

        Modified from Schneider's code (2014)"""

    rsmi = rxn_smiles.split('>')[0].split('.')
    psmi = rxn_smiles.split('>')[2].split('.')

    rfp = None
    pfp = None
    for react in rsmi:
        mol = Chem.MolFromSmiles(react)
        fp = np.array(
            AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=fpsize, useFeatures=useFeatures, useChirality=useChirality))
        if rfp is None:
            rfp = fp
        else:
            rfp += fp
    for product in psmi:
        mol = Chem.MolFromSmiles(product)
        fp = np.array(
            AllChem.GetMorganFingerprintAsBitVect(mol=mol, radius=2, nBits=fpsize, useFeatures=useFeatures, useChirality=useChirality))
        if pfp is None:
            pfp = fp
        else:
            pfp += fp

    if pfp is not None and rfp is not None:
        pfp -= rfp
    return pfp.reshape(1, len(pfp))



def create_rxn_Morgan2FP_separately(rsmi, psmi, rxnfpsize=gc.fingerprint_bits, pfpsize=gc.fingerprint_bits, useFeatures=False, calculate_rfp=True, useChirality=False):
    # Similar as the above function but takes smiles separately and returns pfp and rfp separately

    rsmi = rsmi.encode('utf-8')
    psmi = psmi.encode('utf-8')
    try:
        mol = Chem.MolFromSmiles(rsmi)
    except Exception as e:
        print(e)
        return
    try:
        fp_bit = AllChem.GetMorganFingerprintAsBitVect(
            mol=mol, radius=2, nBits=rxnfpsize, useFeatures=useFeatures, useChirality=useChirality)
        fp = np.empty(rxnfpsize, dtype='float32')
        DataStructs.ConvertToNumpyArray(fp_bit, fp)
    except Exception as e:
        print("Cannot build reactant fp due to {}".format(e))
        return
    rfp = fp

    try:
        mol = Chem.MolFromSmiles(psmi)
    except Exception as e:
        return
    try:
        fp_bit = AllChem.GetMorganFingerprintAsBitVect(
            mol=mol, radius=2, nBits=pfpsize, useFeatures=useFeatures, useChirality=useChirality)
        fp = np.empty(pfpsize, dtype='float32')
        DataStructs.ConvertToNumpyArray(fp_bit, fp)
    except Exception as e:
        print("Cannot build product fp due to {}".format(e))
        return
    pfp = fp
    return [pfp, rfp]


def get_condition_input_from_smiles(conditions_smiles, split=False, s_fp=256, r_fp=256, c_fp=256):
    '''
    If split is used: first molecule in the conditions_smiles should be the solvent!
    '''
    if conditions_smiles == 'NONE':
        return None
    conditions_mol = []
    inputlength = 0
    # If input is a single string: immediately extract molecule
    if (type(conditions_smiles) is str) or (type(conditions_smiles) is unicode):
        try:
            conditions_mol.append(
                ('cond', Chem.MolFromSmiles(conditions_smiles)))
            inputlength += 1
        # If any unparsable: skip this iteration
        except Exception as e:
            MyLogger.print_and_log('Unparsable conditions. Leaving out this reaction: {}.'.format(
                instance['_id']), fingerprinting_loc, level=1)
    # otherwise: assume list of strings
    else:
        for (descr, smiles) in conditions_smiles:
            try:
                conditions_mol.append((descr, Chem.MolFromSmiles(smiles)))
                inputlength += 1
            # If any unparsable: skip this iteration
            except Exception as e:
                MyLogger.print_and_log('Unparsable conditions. Leaving out this reaction: {}.'.format(
                    instance['_id']), fingerprinting_loc, level=1)
                continue

    fps = []
    for (descr, mol) in conditions_mol:
        try:
            if descr == 'solv':
                p = AllChem.GetMorganFingerprintAsBitVect(
                    mol=mol, radius=2, nBits=s_fp, useFeatures=True)
            elif descr == 'reag':
                p = AllChem.GetMorganFingerprintAsBitVect(
                    mol=mol, radius=2, nBits=r_fp, useFeatures=True)
            elif descr == 'cata' or descr == 'cond':
                p = AllChem.GetMorganFingerprintAsBitVect(
                    mol=mol, radius=2, nBits=c_fp, useFeatures=True)
        except Exception as e:
            MyLogger.print_and_log('Could not generate fingerprint for {}'.format(
                conditions_smiles), fingerprinting_loc)
            if descr == 'solv':
                p = [0]*s_fp
            elif descr == 'reag':
                p = [0]*r_fp
            elif descr == 'cata':
                p = [0]*c_fp
        fps.append((descr, p))
    input = []
    if split:
        for (descr, fp) in fps:
            if descr == 'solv':
                fp = np.array(fp).reshape(1, s_fp)
            elif descr == 'reag':
                fp = np.array(fp).reshape(1, r_fp)
            elif descr == 'cata' or descr == 'cond':
                fp = np.array(fp).reshape(1, c_fp)
            input.append(fp)
    else:
        for (descr, fp) in fps:
            input += fp
        input = np.array(input).reshape(1, len(input))
    return input


def get_condition_input_from_instance(instance, chemicals, asone=False, astwo=False, use_new=False, split=False):
    conditions_smiles = get_input_condition_as_smiles(
        instance, chemicals, asone=asone, astwo=astwo, use_new=use_new)
    return get_condition_input_from_smiles(conditions_smiles, split=split)


def get_reaction_input_from_instance(instance, reactions, chemicals):

    reaction_smiles = get_reaction_as_smiles(instance, reactions, chemicals)
    return get_reaction_input_from_smiles(reaction_smiles)


def get_reaction_input_from_smiles(reaction_smiles, r_fp=1024, c_f=1):
    '''
    c_f: compression factor for the reaction fingerprint
    '''
    if reaction_smiles == 'NONE':
        return None
    react = reaction_smiles.split('>')
    reag_s = react[0]
    prod_s = react[len(react)-1]
    reag = Chem.MolFromSmiles(reag_s)
    prod = Chem.MolFromSmiles(prod_s)

    reag = AllChem.GetMorganFingerprintAsBitVect(
        mol=reag, radius=2, nBits=r_fp)
    prod = AllChem.GetMorganFingerprintAsBitVect(
        mol=prod, radius=2, nBits=r_fp)
    reactionfp = [i - j for i, j in zip(prod, reag)]
    reaction = []
    if(r_fp % c_f == 0):
        compressed = int(r_fp/c_f)
        for i in range(compressed):
            pos = 0
            for j in range(c_f):
                pos += reactionfp[i+compressed*j]
            reaction.append(pos)
        input = np.array(reaction).reshape(1, compressed)
    else:
        MyLogger.print_and_log(
            'Redefine reaction fingerprint size or reaction compression ratio. Fingerprint size should be divisible by the compression factor.', fingerprinting_loc, level=3)
    return input


def get_reaction_as_smiles(instance, reactions, chemicals):
    id = instance['RX_ID'][0]
    react = reactions.find_one({'_id': id})
    reaction_smiles = ''
    reactants = react['RX_RXRN']
    products = react['RX_PXRN']
    for reactant in reactants:
        # Try using newest smiles.
        try:
            chemical = chemicals.find_one({'_id': reactant})
            smiles = chemical['SMILES_new']
            mol = Chem.MolFromSmiles(smiles)
            if(smiles and (not mol)):
                return 'NONE'
                '''
                new_smiles = raw_input('New smiles for reactant: {}:'.format(smiles))
                if new_smiles:
                    chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                else:
                    pass
                '''

        except TypeError:
            MyLogger.print_and_log('No reactant smiles found for reaction {} returning "NONE"'.format(
                id), fingerprinting_loc, 2)
            return 'NONE'
        except KeyError:
            # If not present: used older ones
            try:
                chemical = chemicals.find_one({'_id': reactant})
                smiles = chemicals.find_one({'_id': reactant})['SMILES']
                mol = Chem.MolFromSmiles(smiles)
                if(smiles and (not mol)):
                    return 'NONE'
                    '''
                    new_smiles = raw_input('New smiles for reactant: {}:'.format(smiles))
                    if new_smiles:
                        chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                    else:
                        pass
                    '''
            except TypeError:
                MyLogger.print_and_log('No reactant smiles found for reaction {} returning "NONE"'.format(
                    id), fingerprinting_loc, 2)
                return 'NONE'
        if not reaction_smiles:
            reaction_smiles = smiles
        else:
            reaction_smiles += '.'+smiles
    reaction_smiles += '>>'
    for product in products:
        try:
            chemical = chemicals.find_one({'_id': product})
            smiles = chemical['SMILES_new']
            mol = Chem.MolFromSmiles(smiles)
            if(smiles and (not mol)):
                return 'NONE'
                '''
                new_smiles = raw_input('New smiles for product: {}'.format(smiles))
                if new_smiles:
                    chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                else:
                    pass
                '''
        except TypeError:
            MyLogger.print_and_log('No product smiles found for reaction {} returning "NONE"'.format(
                id), fingerprinting_loc, 2)
            return 'NONE'
        except KeyError:
            try:
                chemical = chemicals.find_one({'_id': product})
                smiles = chemicals.find_one({'_id': product})['SMILES']
                mol = Chem.MolFromSmiles(smiles)
                if(smiles and (not mol)):
                    return 'NONE'
                    '''
                    new_smiles = raw_input('New smiles for product: {}'.format(smiles))
                    if new_smiles:
                        chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles,'checked':True}})
                    else:
                        pass
                    '''
            except TypeError:
                MyLogger.print_and_log('No product smiles found for reaction {} returning "NONE"'.format(
                    id), fingerprinting_loc, 2)
                return 'NONE'
        if reaction_smiles[len(reaction_smiles)-1] == '>':
            reaction_smiles += smiles
        else:
            reaction_smiles += '.'+smiles
    return reaction_smiles


def get_input_condition_as_smiles(doc, chemicals, asone=False, astwo=False, check=False, use_new=False):
    a = doc['RXD_SOLXRN']
    b = doc['RXD_RGTXRN']
    c = doc['RXD_CATXRN']
    solv = ""
    reag = ""
    cata = ""
    for sol in a:
        sa = chemicals.find_one({'_id': sol})
        if sa:
            if check:
                check_smiles(sa, chemicals)
                s = None
                try:
                    s = sa['SMILES_new']
                except:
                    s = sa['SMILES']
            elif use_new:
                try:
                    s = sa['SMILES_new']
                except:
                    s = sa['SMILES']
            else:
                s = sa['SMILES']

            mol = Chem.MolFromSmiles(s)
            if not mol:
                s = 'NONE'
                '''
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
        sa = chemicals.find_one({'_id': rgt})
        if sa:
            if check:
                check_smiles(sa, chemicals)
                s = None
                try:
                    s = sa['SMILES_new']
                except:
                    s = sa['SMILES']
            elif use_new:
                try:
                    s = sa['SMILES_new']
                except:
                    s = sa['SMILES']
            else:
                s = sa['SMILES']

            mol = Chem.MolFromSmiles(s)
            if not mol:
                s = 'NONE'
                '''
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
        sa = chemicals.find_one({'_id': cat})
        if sa:
            if check:
                check_smiles(sa, chemicals)
                s = None
                try:
                    s = sa['SMILES_new']
                except:
                    s = sa['SMILES']
            elif use_new:
                try:
                    s = sa['SMILES_new']
                except:
                    s = sa['SMILES']
            else:
                s = sa['SMILES']

            mol = Chem.MolFromSmiles(s)
            if not mol:
                s = 'NONE'
                '''
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
        return ans.rstrip('.')
    if astwo:
        ans = ""
        if reag:
            ans = reag
        if cata:
            if ans:
                ans += '.'+cata
            else:
                ans = cata
        return [('solv', solv.rstrip('.').replace('..', '.').replace('..', '.')), ('reag', ans.rstrip('.').replace('..', '.').replace('..', '.'))]
    else:
        return [('solv', solv.rstrip('.').replace('..', '.').replace('..', '.')), ('reag', reag.rstrip('.').replace('..', '.').replace('..', '.')), ('cata', cata.rstrip('.').replace('..', '.').replace('..', '.'))]
