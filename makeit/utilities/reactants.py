import makeit.global_config as gc
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from makeit.utilities.io.logger import MyLogger
reactants_loc = 'util.reactants'


def clean_reactant_mapping(reactants):
    """Remaps atoms for reactants.

    Args:
        reactants (Chem.Mol): Reactants to remap.

    Returns:
        Chem.Mol: Reactants with remapped atoms.
    """
    if not reactants:
        MyLogger.print_and_log('Could not parse reactants {}'.format(reactants),reactants_loc)
        raise ValueError('Could not parse reactants')
    if gc.DEBUG: print('Number of reactant atoms: {}'.format(len(reactants.GetAtoms())))
    # Report current reactant SMILES string
    [a.ClearProp('molAtomMapNumber') for a in reactants.GetAtoms() if a.HasProp('molAtomMapNumber')]
    if gc.DEBUG: print('Reactants w/o map: {}'.format(Chem.MolToSmiles(reactants)))
    # Add new atom map numbers
    [a.SetProp('molAtomMapNumber', str(i+1)) for (i, a) in enumerate(reactants.GetAtoms())]
    # Report new reactant SMILES string
    if gc.DEBUG: print('Reactants w/ map: {}'.format(Chem.MolToSmiles(reactants)))
    return reactants
