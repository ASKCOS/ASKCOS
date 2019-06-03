import sys 
if sys.version_info[0] < 3:
    from urllib2 import urlopen 
else:
    from urllib.request import urlopen
import rdkit.Chem as Chem
name_parser_loc = 'name_parser'


def name_to_molecule(name):
    try:
        mol = Chem.MolFromSmiles(name)
        if not mol:
            raise ValueError
        return mol
    except:
        pass

    smiles = urlopen(
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/IsomericSMILES/txt'.format(name)).read()
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(
            'Could not resolve SMILES ({}) in a way parseable by RDKit, from identifier: {}'.format(smiles, name))

    return mol
