import urllib2
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

    smiles = urllib2.urlopen(
        'https://cactus.nci.nih.gov/chemical/structure/{}/smiles'.format(name)).read()
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError(
            'Could not resolve SMILES ({}) in a way parseable by RDKit, from identifier: {}'.format(smiles, name))

    return mol
