import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir, BondStereo

from rdchiral.chiral import template_atom_could_have_been_tetra
from rdchiral.utils import vprint, PLEVEL
from rdchiral.bonds import enumerate_possible_cistrans_defs, bond_dirs_by_mapnum, \
    get_atoms_across_double_bonds

BondDirOpposite = {AllChem.BondDir.ENDUPRIGHT: AllChem.BondDir.ENDDOWNRIGHT,
                   AllChem.BondDir.ENDDOWNRIGHT: AllChem.BondDir.ENDUPRIGHT}

class rdchiralReaction():
    '''
    Class to store everything that should be pre-computed for a reaction. This
    makes library application much faster, since we can pre-do a lot of work
    instead of doing it for every mol-template pair
    '''
    def __init__(self, reaction_smarts):
        # Keep smarts, useful for reporting
        self.reaction_smarts = reaction_smarts
        
        # Initialize - assigns stereochemistry and fills in missing rct map numbers
        self.rxn = initialize_rxn_from_smarts(reaction_smarts)

        # Combine template fragments so we can play around with isotopes
        self.template_r, self.template_p = get_template_frags_from_rxn(self.rxn)

        # Define molAtomMapNumber->atom dictionary for template rct and prd
        self.atoms_rt_map = {a.GetIntProp('molAtomMapNumber'): a \
            for a in self.template_r.GetAtoms() if a.HasProp('molAtomMapNumber')}
        self.atoms_pt_map = {a.GetIntProp('molAtomMapNumber'): a \
            for a in self.template_p.GetAtoms() if a.HasProp('molAtomMapNumber')}

        # Check consistency (this should not be necessary...)
        if any(self.atoms_rt_map[i].GetAtomicNum() != self.atoms_pt_map[i].GetAtomicNum() \
                for i in self.atoms_rt_map if i in self.atoms_pt_map):
            raise ValueError('Atomic identity should not change in a reaction!')

        # Call template_atom_could_have_been_tetra to pre-assign value to atom
        [template_atom_could_have_been_tetra(a) for a in self.template_r.GetAtoms()]
        [template_atom_could_have_been_tetra(a) for a in self.template_p.GetAtoms()]

        # Pre-list chiral double bonds (for copying back into outcomes/matching)
        self.rt_bond_dirs_by_mapnum = bond_dirs_by_mapnum(self.template_r)
        self.pt_bond_dirs_by_mapnum = bond_dirs_by_mapnum(self.template_p)

        # Enumerate possible cis/trans...
        self.required_rt_bond_defs, self.required_bond_defs_coreatoms = \
            enumerate_possible_cistrans_defs(self.template_r)

class rdchiralReactants():
    '''
    Class to store everything that should be pre-computed for a reactant mol
    so that library application is faster
    '''
    def __init__(self, reactant_smiles):
        # Keep original smiles, useful for reporting
        self.reactant_smiles = reactant_smiles  

        # Initialize into RDKit mol
        self.reactants = initialize_reactants_from_smiles(reactant_smiles)

        # Set isotope->atom dictionary
        # all reactant atoms must be mapped after initialization, so this is safe
        self.atoms_r = {a.GetIsotope(): a for a in self.reactants.GetAtoms()}

        # Create copy of molecule without chiral information, used with
        # RDKit's naive runReactants
        self.reactants_achiral = initialize_reactants_from_smiles(reactant_smiles)
        [a.SetChiralTag(ChiralType.CHI_UNSPECIFIED) for a in self.reactants_achiral.GetAtoms()]
        [(b.SetStereo(BondStereo.STEREONONE), b.SetBondDir(BondDir.NONE)) \
            for b in self.reactants_achiral.GetBonds()]

        # Pre-list reactant bonds (for stitching broken products)
        self.bonds_by_isotope = [
            (b.GetBeginAtom().GetIsotope(), b.GetEndAtom().GetIsotope(), b) \
            for b in self.reactants.GetBonds()
        ]

        # Pre-list chiral double bonds (for copying back into outcomes/matching)
        self.bond_dirs_by_isotope = {}
        for (i, j, b) in self.bonds_by_isotope:
            if b.GetBondDir() != BondDir.NONE:
                self.bond_dirs_by_isotope[(i, j)] = b.GetBondDir()
                self.bond_dirs_by_isotope[(j, i)] = BondDirOpposite[b.GetBondDir()]

        # Get atoms across double bonds defined by isotope
        self.atoms_across_double_bonds = get_atoms_across_double_bonds(self.reactants)


def initialize_rxn_from_smarts(reaction_smarts):
    # Initialize reaction
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    rxn.Initialize()
    if rxn.Validate()[1] != 0:
        raise ValueError('validation failed')
    if PLEVEL >= 2: print('Validated rxn without errors')

    unmapped = 700
    for rct in rxn.GetReactants():
        rct.UpdatePropertyCache(strict=False)
        Chem.AssignStereochemistry(rct)
        # Fill in atom map numbers
        for a in rct.GetAtoms():
            if not a.HasProp('molAtomMapNumber'):
                a.SetIntProp('molAtomMapNumber', unmapped)
                unmapped += 1
    if PLEVEL >= 2: print('Added {} map nums to unmapped reactants'.format(unmapped-700))
    if unmapped > 800:
        raise ValueError('Why do you have so many unmapped atoms in the template reactants?')

    return rxn

def initialize_reactants_from_smiles(reactant_smiles):
    # Initialize reactants
    reactants = Chem.MolFromSmiles(reactant_smiles)
    Chem.AssignStereochemistry(reactants, flagPossibleStereoCenters=True)
    reactants.UpdatePropertyCache()
    # To have the product atoms match reactant atoms, we
    # need to populate the Isotope field, since this field
    # gets copied over during the reaction.
    [a.SetIsotope(i+1) for (i, a) in enumerate(reactants.GetAtoms())]
    if PLEVEL >= 2: print('Initialized reactants, assigned isotopes, stereochem, flagpossiblestereocenters')
    return reactants

def get_template_frags_from_rxn(rxn):
    # Copy reaction template so we can play around with isotopes
    for i, rct in enumerate(rxn.GetReactants()):
        if i == 0:
            template_r = rct
        else:
            template_r = AllChem.CombineMols(template_r, rct)
    for i, prd in enumerate(rxn.GetProducts()):
        if i == 0:
            template_p = prd
        else:
            template_p = AllChem.CombineMols(template_p, prd)
    return template_r, template_p