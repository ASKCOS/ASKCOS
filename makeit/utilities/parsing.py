import makeit.global_config as gc
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from makeit.utilities.io.logger import MyLogger
parsing_loc = 'parsing'

def canonicalize_smiles(smiles):
    """Returns canonicalized SMILES string.

    Args:
        smiles (str): SMILES string to canonicalize.
    """
    #rdkit defaults canonicalize to true
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

def parse_molecule_to_smiles(target):
    """Parses a molecular type (smiles, rdkit mol or mol file) to smiles format).

    Args:
        target (str or Chem.Mol): SMILES string, filename, or Chem.Mol to parse.

    Returns:
        str or None: SMILES string of target, or None if parsing fails.
    """
    try:
        mol = Chem.MolFromSmiles(target)
        if mol:
            #This in order to canonicalize the molecule
            return Chem.MolToSmiles(mol)
    except Exception as e:
        try:
            smiles = Chem.MolToSmiles(target, isomericSmiles = gc.USE_STEREOCHEMISTRY)
            if smiles:
                return smiles
        except Exception as e:
            try:
                mol = Chem.MolFromMolFile(target)
                if mol:
                    return Chem.MolToSmiles(mol, isomericSmiles = gc.USE_STEREOCHEMISTRY)
            except Exception as e:
                MyLogger.print_and_log('Unable to parse target molecule format. Parsing Only available for: Smiles, RDKIT molecule and mol files. Returning "None"', parsing_loc, level = 1)
    return None
def check_smiles(chemical, chemicals):
    """Checks recorded SMILES for a chemical against the SMILES from its name.

    Args:
        chemical (pymongo.collection): Chemical to check the SMILES of.
        chemicals (pymongo.collection): Database of recorded chemicals.
    """
    is_checked = False
    try:
        is_checked = chemical['checked']
    except:
        is_checked = False
        pass
    if(is_checked):
        return
    recorded_smiles = chemical['SMILES']
    trivial_name = None
    try:
        trivial_name = chemical['IDE_CN']
    except KeyError:
        return
    smiles_from_name = None

    #Short cut on not fully interpreted chemicals
    if trivial_name[0]=='<':
        #don't bother, skips anyway: is faster.
        #chemicals.update({'_id':chemical['_id']},{'$set':{'checked':True}})
        return
    #make spaces readable in url
    trivial_name = trivial_name.replace(' ','%20')
    try:
        smiles_from_name = urllib2.urlopen('https://cactus.nci.nih.gov/chemical/structure/{}/smiles'.format(trivial_name)).read()
    except Exception as e:
        pass

    try:
        recorded_can_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(recorded_smiles))
        can_smiles_from_name = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_from_name))
    except Exception as e:
        #don't bother, skips anyway: is faster.
        #chemicals.update({'_id':chemical['_id']},{'$set':{'checked':True}})
        return

    if(recorded_can_smiles != can_smiles_from_name):
        print('Chemical "{}" with reaxys id {}:'.format(trivial_name.replace('%20',' '), chemical['_id']))
        print('Recorded smiles: {}'.format(recorded_can_smiles))
        print('Smiles from recorded name: {}'.format(can_smiles_from_name))
        ans = None
        while not ans:
            #ans = raw_input('Should the recorded smiles be overwritten by smiles from recorded name? (Y)es or (N)o\t')
            #just write them all in new smiles.
            ans = 'Y'
            if ans == 'Y' or ans == 'y':
                chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':can_smiles_from_name,'checked':True}})
                print('Recorded smiles changed to: {}'.format(can_smiles_from_name))
            elif ans == 'N' or ans == 'n':
                change = None
                while not change:
                    change = raw_input('Do you want to manually enter a smiles string? (Y)es or (N)o\t')
                    if change == 'y' or change == 'Y':
                        new_smiles = raw_input('Please type the desired smiles:\t')
                        chemicals.update({'_id':chemical['_id']},{'$set':{'SMILES_new':new_smiles, 'checked':True}})
                        print('Smiles have been changed to: {}'.format(new_smiles))
                    elif change == 'n' or change == 'N':
                        print('Recorded smiles will be kept.')
                        chemicals.update({'_id':chemical['_id']},{'$set':{'checked':True}})
                    else:
                        print('Please answer with one of the given options.')
                        change = None

            else:
                print('Please answer with one of the given options.')
                ans = None
    chemicals.update({'_id':chemical['_id']},{'$set':{'checked':True}})

def parse_list_to_smiles(mol_list):
    """Parses reactants into SMILES string.

    Format of reactants: list of smiles, rdkit mol, or mol files; or
    single smiles, rdkit mol, or mol file.

    Args:
        mol_list (str, Chem.Mol, or list of str or Chem.Mol): Reactants to be
            parsed into SMILES strings.

    Returns:
        str or None: SMILES string of reactants or None if parsing fails.
    """
    isList = False
    if mol_list.__class__.__name__ == 'list':
        isList = True

    smiles = ''
    if isList:
        for mol in mol_list:
            # QUESTION: Won't this raise an exception if mol cannot be parsed?
            if smiles:
                smiles += '.' + parse_molecule_to_smiles(mol)
            else:
                smiles += parse_molecule_to_smiles(mol)
    else:
        smiles = parse_molecule_to_smiles(mol_list)

    return smiles
