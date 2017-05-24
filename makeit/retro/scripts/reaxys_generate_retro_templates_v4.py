'''
This script (generate_reaction_templates) uses a MongoDB collection 
with atom-mapped reaction SMILES strings and parses them into a new 
collection containing the transforms.

This is intended to be used with the Reaxys database. In the database,
reagents can contribute atoms to the products. This means that those
atoms are not mapped in the RXN_SMILES field. The script currently
leaves those atoms OUT of the template.

As an example, halogenation might be performed using [Cl][Cl] as a
chlorinating agent, so the chlorine atom in the product will be 
unmapped. This script will create a retrosynthetic template that does
not included a specific precursor containing a Cl atom. Instead, an
extra field is added to the template document indicating that there 
is a necessary_reagent fragment (as a generalized SMARTS string).

Additionally, in the cases of unmapped product atoms, those atoms are
FULLY specified in the product fragment

v4 does this in parallel
'''

from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
from numpy.random import shuffle # for random selection
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
from collections import defaultdict
import rdkit.Chem.Draw as Draw
from rdkit import RDLogger
import datetime # for info files
import json # for dumping
import sys  # for commanad line
import os   # for file paths
import re 
import itertools
from makeit.retro.draw import *
from tqdm import tqdm 

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
TRANSFORM_DB = db['transforms_retro_v5']
REACTION_DB = db['reactions']
INSTANCE_DB = db['instances']
CHEMICAL_DB = db['chemicals']
reaction_smiles_field = 'RXN_SMILES'
MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS = 5

def mols_from_smiles_list(all_smiles):
    '''Given a list of smiles strings, this function creates rdkit
    molecules'''
    mols = []
    for smiles in all_smiles:
        if not smiles: continue
        mols.append(Chem.MolFromSmiles(smiles))
    return mols

def bond_to_label(bond):
    '''This function takes an RDKit bond and creates a label describing
    the most important attributes'''
    # atoms = sorted([atom_to_label(bond.GetBeginAtom().GetIdx()), \
    #               atom_to_label(bond.GetEndAtom().GetIdx())])
    a1_label = str(bond.GetBeginAtom().GetAtomicNum())
    a2_label = str(bond.GetEndAtom().GetAtomicNum())
    if bond.GetBeginAtom().HasProp('molAtomMapNumber'):
        a1_label += bond.GetBeginAtom().GetProp('molAtomMapNumber')
    if bond.GetEndAtom().HasProp('molAtomMapNumber'):
        a2_label += bond.GetEndAtom().GetProp('molAtomMapNumber')
    atoms = sorted([a1_label, a2_label])

    return '{}{}{}'.format(atoms[0], bond.GetSmarts(), atoms[1])

def get_tagged_atoms_from_mols(mols):
    '''Takes a list of RDKit molecules and returns total list of
    atoms and their tags'''
    atoms = []
    atom_tags = []
    for mol in mols:
        new_atoms, new_atom_tags = get_tagged_atoms_from_mol(mol)
        atoms += new_atoms 
        atom_tags += new_atom_tags
    return atoms, atom_tags

def get_tagged_atoms_from_mol(mol):
    '''Takes an RDKit molecule and returns list of tagged atoms and their
    corresponding numbers'''
    atoms = []
    atom_tags = []
    for atom in mol.GetAtoms():
        if ':' in atom.GetSmarts():
            atoms.append(atom)
            atom_tags.append(atom.GetSmarts().split(':')[1][:-1])
    return atoms, atom_tags

def atoms_are_different(atom1, atom2, level = 1):
    '''Compares two RDKit atoms based on basic properties'''

    if atom1.GetSmarts() != atom2.GetSmarts(): return True # should be very general
    if atom1.GetAtomicNum() != atom2.GetAtomicNum(): return True # must be true for atom mapping
    if atom1.GetTotalNumHs() != atom2.GetTotalNumHs(): return True
    if atom1.GetFormalCharge() != atom2.GetFormalCharge(): return True
    if atom1.GetDegree() != atom2.GetDegree(): return True
    #if atom1.IsInRing() != atom2.IsInRing(): return True
    if atom1.GetNumRadicalElectrons() != atom2.GetNumRadicalElectrons(): return True
    # TODO: add # pi electrons like ICSynth?

    # Check bonds and nearest neighbor identity
    if level >= 1:
        bonds1 = sorted([bond_to_label(bond) for bond in atom1.GetBonds()]) 
        bonds2 = sorted([bond_to_label(bond) for bond in atom2.GetBonds()]) 
        if bonds1 != bonds2: return True

        # # Check neighbors too (already taken care of with previous lines)
        # neighbors1 = sorted([atom.GetAtomicNum() for atom in atom1.GetNeighbors()])
        # neighbors2 = sorted([atom.GetAtomicNum() for atom in atom2.GetNeighbors()])
        # if neighbors1 != neighbors2: return True

    return False

def get_changed_atoms(reactants, products):
    '''Looks at mapped atoms in a reaction and determines which ones changed'''

    err = 0
    prod_atoms, prod_atom_tags = get_tagged_atoms_from_mols(products)

    if v: print('Products contain {} tagged atoms'.format(len(prod_atoms)))
    if v: print('Products contain {} unique atom numbers'.format(len(set(prod_atom_tags))))

    reac_atoms, reac_atom_tags = get_tagged_atoms_from_mols(reactants)
    if len(set(prod_atom_tags)) != len(set(reac_atom_tags)):
        if v: print('warning: different atom tags appear in reactants and products')
        #err = 1 # okay for Reaxys, since Reaxys creates mass
    if len(prod_atoms) != len(reac_atoms):
        if v: print('warning: total number of tagged atoms differ, stoichometry != 1?')
        #err = 1

    # Find differences 
    changed_atoms = []
    changed_atom_tags = []
    #print(reac_atom_tags)
    #print(prod_atom_tags)

    # Product atoms that are different from reactant atom equivalent
    for i, prod_tag in enumerate(prod_atom_tags):

        for j, reac_tag in enumerate(reac_atom_tags):
            if reac_tag != prod_tag: continue
            if reac_tag not in changed_atom_tags: # don't bother comparing if we know this atom changes
                # If atom changed, add
                if atoms_are_different(prod_atoms[i], reac_atoms[j]):
                    changed_atoms.append(reac_atoms[j])
                    changed_atom_tags.append(reac_tag)
                    break
                # If reac_tag appears multiple times, add (need for stoichometry > 1)
                if prod_atom_tags.count(reac_tag) > 1:
                    changed_atoms.append(reac_atoms[j])
                    changed_atom_tags.append(reac_tag)
                    break

    # Reactant atoms that do not appear in product (tagged leaving groups)
    for j, reac_tag in enumerate(reac_atom_tags):
        if reac_tag not in changed_atom_tags:
            if reac_tag not in prod_atom_tags:
                changed_atoms.append(reac_atoms[j])
                changed_atom_tags.append(reac_tag)

    if v: 
        print('{} tagged atoms in reactants change 1-atom properties'.format(len(changed_atom_tags)))
        for smarts in [atom.GetSmarts() for atom in changed_atoms]:
            print('  {}'.format(smarts))

    return changed_atoms, changed_atom_tags, err

def expand_atoms_to_use(mol, atoms_to_use, groups=[], symbol_replacements=[]):
    '''Given an RDKit molecule and a list of AtomIdX which should be included
    in the reaction, this function expands the list of AtomIdXs to include one 
    nearest neighbor with special consideration of (a) unimportant neighbors and
    (b) important functional groupings'''

    # Copy
    new_atoms_to_use = atoms_to_use[:]

    # Look for all atoms in the current list of atoms to use
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in atoms_to_use: continue
        # Ensure membership of changed atom is checked against group
        for group in groups:
            if int(atom.GetIdx()) in group[0]:
                if v: 
                    print('adding group due to match')
                    try:
                        print('Match from molAtomMapNum {}'.format(
                            atom.GetProp('molAtomMapNumber'),
                        ))
                    except KeyError:
                        pass
                for idx in group[1]:
                    if idx not in atoms_to_use:
                        new_atoms_to_use.append(idx)
                        symbol_replacements.append((idx, convert_atom_to_wildcard(mol.GetAtomWithIdx(idx))))
        # Look for all nearest neighbors of the currently-included atoms
        for neighbor in atom.GetNeighbors():
            # Evaluate nearest neighbor atom to determine what should be included
            new_atoms_to_use, symbol_replacements = \
                    expand_atoms_to_use_atom(mol, new_atoms_to_use, neighbor.GetIdx(), 
                        groups=groups, symbol_replacements=symbol_replacements)
            
    return new_atoms_to_use, symbol_replacements

def expand_atoms_to_use_atom(mol, atoms_to_use, atom_idx, groups=[], symbol_replacements=[]):
    '''Given an RDKit molecule and a list of AtomIdx which should be included
    in the reaction, this function extends the list of atoms_to_use by considering 
    a candidate atom extension, atom_idx'''

    # See if this atom belongs to any groups
    found_in_group = False
    for group in groups: # first index is atom IDs for match, second is what to include
        if int(atom_idx) in group[0]: # int correction
            if v: 
                print('adding group due to match')
                try:
                    print('Match from molAtomMapNum {}'.format(
                        mol.GetAtomWithIdx(atom_idx).GetProp('molAtomMapNumber'),
                    ))
                except KeyError:
                    pass
            # Add the whole list, redundancies don't matter 
            # *but* still call convert_atom_to_wildcard!
            for idx in group[1]:
                if idx not in atoms_to_use:
                        atoms_to_use.append(idx)
                        symbol_replacements.append((idx, convert_atom_to_wildcard(mol.GetAtomWithIdx(idx))))
            found_in_group = True
    if found_in_group:  
        return atoms_to_use, symbol_replacements
        
    # How do we add an atom that wasn't in an identified important functional group?
    # Develop special SMARTS symbol

    # Skip current candidate atom if it is already included
    if atom_idx in atoms_to_use:
        return atoms_to_use, symbol_replacements

    # Include this atom
    atoms_to_use.append(atom_idx)

    # Look for replacements
    symbol_replacements.append((atom_idx, convert_atom_to_wildcard(mol.GetAtomWithIdx(atom_idx))))

    return atoms_to_use, symbol_replacements

def convert_atom_to_wildcard(atom):
    '''This function takes an RDKit atom and turns it into a wildcard 
    using hard-coded generalization rules. This function should be used
    when candidate atoms are used to extend the reaction core for higher
    generalizability'''

    # Is this a terminal atom? We can tell if the degree is one
    if atom.GetDegree() == 1:
        symbol = '[' + atom.GetSymbol() + ';X1;H{}'.format(atom.GetTotalNumHs())
    else:

        # Initialize
        symbol = '['

        # Add atom primitive (don't use COMPLETE wildcards)
        if atom.GetAtomicNum() != 6:
            symbol += '#{};'.format(atom.GetAtomicNum())
            if atom.GetIsAromatic():
                symbol += 'a;'
        elif atom.GetIsAromatic():
            symbol += 'c;'
        else:
            symbol += 'C;'

        # Charge is important
        if atom.GetFormalCharge() != 0:
            charges = re.search('([-+]+[1-9]?)', atom.GetSmarts())
            if charges: symbol += charges.group() + ';'

        # Strip extra semicolon
        if symbol[-1] == ';': symbol = symbol[:-1]

    # Close with label or with bracket
    label = re.search('\:[0-9]+\]', atom.GetSmarts())
    if label: 
        symbol += label.group()
    else:
        symbol += ']'

    if v: 
        if symbol != atom.GetSmarts():
            print('Improved generality of atom SMARTS {} -> {}'.format(atom.GetSmarts(), symbol))

    return symbol

def get_strict_smarts_for_atom(atom):
    '''
    For an RDkit atom object, generate a SMARTS pattern that
    matches the atom as strictly as possible
    '''

    symbol = atom.GetSmarts()
    if atom.GetSymbol() == 'H':
        symbol = '[#1]'

    if '[' not in symbol:
        symbol = '[' + symbol + ']'

    if 'H' not in symbol:
        H_symbol = 'H{}'.format(atom.GetTotalNumHs())
        # Explicit number of hydrogens
        if ':' in symbol: # stick H0 before label
            symbol = symbol.replace(':', ';{}:'.format(H_symbol))
        else:
            symbol = symbol.replace(']', ';{}]'.format(H_symbol))
      
    # Explicit degree
    if ':' in symbol:
        symbol = symbol.replace(':', ';D{}:'.format(atom.GetDegree()))
    else:
        symbol = symbol.replace(']', ';D{}]'.format(atom.GetDegree()))

    # Explicit formal charge
    if '+' not in symbol and '-' not in symbol:
        charge = atom.GetFormalCharge()
        charge_symbol = '+' if (charge >= 0) else '-'
        charge_symbol += '{}'.format(abs(charge))
        if ':' in symbol: 
            symbol = symbol.replace(':', ';{}:'.format(charge_symbol))
        else:
            symbol = symbol.replace(']', ';{}]'.format(charge_symbol))

    # Explicit stereochemistry
    if USE_STEREOCHEMISTRY:
        if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            if '@' not in symbol:
                # Be explicit when there is a tetrahedral chiral tag
                if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    tag = '@'
                elif atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                    tag = '@@'
                if ':' in symbol:
                    symbol = symbol.replace(':', ';{}:'.format(tag))
                else:
                    symbol = symbol.replace(']', ';{}]'.format(tag))

    return symbol

def get_fragments_for_changed_atoms(mols, changed_atom_tags, radius=0, 
    category='reactants', expansion=[]):
    '''Given a list of RDKit mols and a list of changed atom tags, this function
    computes the SMILES string of molecular fragments using MolFragmentToSmiles 
    for all changed fragments.

    expansion: atoms added during reactant expansion that should be included and
               generalized in product fragment
    '''

    fragments = ''
    for mol in mols:
        # Initialize list of replacement symbols (updated during expansion)
        symbol_replacements = []

        # Are we looking for groups? (reactants only)
        if category == 'reactants':
            groups = get_special_groups(mol)
        else:
            groups = []

        # Build list of atoms to use
        atoms_to_use = []
        for atom in mol.GetAtoms():
            # Check self (only tagged atoms)
            if ':' in atom.GetSmarts():
                if atom.GetSmarts().split(':')[1][:-1] in changed_atom_tags:
                    atoms_to_use.append(atom.GetIdx())
                    symbol = get_strict_smarts_for_atom(atom)
                    if symbol != atom.GetSmarts():
                        symbol_replacements.append((atom.GetIdx(), symbol))
                    continue

        # Check neighbors (any atom)
        for k in range(radius):
            atoms_to_use, symbol_replacements = expand_atoms_to_use(mol, atoms_to_use, 
                groups=groups, symbol_replacements=symbol_replacements)

        if category == 'products':
            # Add extra labels to include (for products only)
            if expansion:
                for atom in mol.GetAtoms():
                    if ':' not in atom.GetSmarts(): continue
                    label = atom.GetSmarts().split(':')[1][:-1]
                    if label in expansion and label not in changed_atom_tags:
                        atoms_to_use.append(atom.GetIdx())
                        # Make the expansion a wildcard
                        symbol_replacements.append((atom.GetIdx(), convert_atom_to_wildcard(atom))) 
                        if v: print('expanded label {} to wildcard in products'.format(label))
            
            # Make sure unmapped atoms are included (from products)
            for atom in mol.GetAtoms():
                if not atom.HasProp('molAtomMapNumber'): 
                    atoms_to_use.append(atom.GetIdx())
                    symbol = get_strict_smarts_for_atom(atom)
                    symbol_replacements.append((atom.GetIdx(), symbol))

        # Define new symbols based on symbol_replacements
        symbols = [atom.GetSmarts() for atom in mol.GetAtoms()]
        for (i, symbol) in symbol_replacements:
            symbols[i] = symbol

        if not atoms_to_use: continue
        # if v:
        #   print('~~ replacement for this ' + category[:-1])
        #   print('{} -> {}'.format([mol.GetAtomWithIdx(x).GetSmarts() for (x, s) in symbol_replacements], 
        #                           [s for (x, s) in symbol_replacements]))
        # Remove molAtomMapNumber before canonicalization
        [x.ClearProp('molAtomMapNumber') for x in mol.GetAtoms()]
        
        fragments += '(' + AllChem.MolFragmentToSmiles(mol, atoms_to_use, 
        atomSymbols = symbols, allHsExplicit = True, 
        isomericSmiles = USE_STEREOCHEMISTRY, allBondsExplicit = True) + ').'
        
    return fragments[:-1]

def expand_changed_atom_tags(changed_atom_tags, reactant_fragments):
    '''Given a list of changed atom tags (numbers as strings) and a string consisting
    of the reactant_fragments to include in the reaction transform, this function 
    adds any tagged atoms found in the reactant side of the template to the 
    changed_atom_tags list so that those tagged atoms are included in the products'''

    expansion = []
    atom_tags_in_reactant_fragments = re.findall('\:([[0-9]+)\]', reactant_fragments)
    for atom_tag in atom_tags_in_reactant_fragments:
        if atom_tag not in changed_atom_tags:
            expansion.append(atom_tag)
    if v: print('after building reactant fragments, additional labels included: {}'.format(expansion))
    return expansion

def get_special_groups(mol):
    '''Given an RDKit molecule, this function returns a list of tuples, where
    each tuple contains the AtomIdx's for a special group of atoms which should 
    be included in a fragment all together. This should only be done for the 
    reactants, otherwise the products might end up with mapping mismatches

    We draw a distinction between atoms in groups that trigger that whole
    group to be included, and "unimportant" atoms in the groups that will not
    be included if another atom matches.'''

    # Define templates, based on Functional_Group_Hierarchy.txt from Greg Laandrum
    group_templates = [ 
        (range(3), '[OH0,SH0]=C[O,Cl,I,Br,F]',), # carboxylic acid / halogen
        (range(3), '[OH0,SH0]=CN',), # amide/sulfamide
        (range(4), 'S(O)(O)[Cl]',), # sulfonyl chloride
        (range(3), 'B(O)O',), # boronic acid/ester
        (range(2), 'N#C',), # nitrile
        (range(3), '[N;H0;$(N-[#6]);D2]-,=[N;D2]-,=[N;D1]',), # azide
        (range(8), 'O=C1N([Br,I,F,Cl])C(=O)CC1',), # NBS brominating agent
        (range(2), 'C=C',), # alkene
        ((2,), 'C=C-[*]',), # adj to alkene
        ((2,), 'O=C-[*]',), # adj to carbonyl
        ((3,), 'O=C([O,N])-[*]',), # adj to carboxylic acid/amide/ester
        (range(4), 'ClS(Cl)=O',), # thionyl chloride
        (range(2), '[Mg,Li][Br,Cl]',), # grinard (non-disassociated)
        (range(3), 'S(O)(O)',), # SO2 group
        (range(2), 'N~N',), # diazo
        ((1,), '[a!c]:a',), # adjacency to heteroatom in aromatic ring
        ((2,), '[a!c]:a:a',), # two-steps away from heteroatom in aromatic ring
        #((1,), 'c(-,=[*]):c([Cl,I,Br,F])',), # ortho to halogen on ring - too specific?
        #((1,), 'c(-,=[*]):c:c([Cl,I,Br,F])',), # meta to halogen on ring - too specific?
        ]

    # Build list
    groups = []
    for (add_if_match, template) in group_templates:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(template))
        for match in matches:
            add_if = []
            for pattern_idx, atom_idx in enumerate(match):
                if pattern_idx in add_if_match:
                    add_if.append(atom_idx)
            groups.append((add_if, match))
    return groups

def canonicalize_template(template):
    '''This function takes one-half of a template SMARTS string 
    (i.e., reactants or products) and re-orders them based on
    an equivalent string without atom mapping.'''

    # Strip labels to get sort orders
    template_nolabels = re.sub('\:[0-9]+\]', ']', template)

    # Split into separate molecules *WITHOUT wrapper parentheses*
    template_nolabels_mols = template_nolabels[1:-1].split(').(')
    template_mols          = template[1:-1].split(').(')

    # Split into fragments within those molecules
    for i in range(len(template_mols)):
        nolabel_mol_frags = template_nolabels_mols[i].split('.')
        mol_frags         = template_mols[i].split('.')

        # Get sort order within molecule, defined WITHOUT labels
        sortorder = [j[0] for j in sorted(enumerate(nolabel_mol_frags), key = lambda x:x[1])]

        # Apply sorting and merge list back into overall mol fragment
        template_nolabels_mols[i] = '.'.join([nolabel_mol_frags[j] for j in sortorder])
        template_mols[i]          = '.'.join([mol_frags[j] for j in sortorder])

    # Get sort order between molecules, defined WITHOUT labels
    sortorder = [j[0] for j in sorted(enumerate(template_nolabels_mols), key = lambda x:x[1])]

    # Apply sorting and merge list back into overall transform
    template = '(' + ').('.join([template_mols[i] for i in sortorder]) + ')'

    return template

def reassign_atom_mapping(transform):
    '''This function takes an atom-mapped reaction and reassigns 
    the atom-mapping labels (numbers) from left to right, once 
    that transform has been canonicalized.'''

    all_labels = re.findall('\:([0-9]+)\]', transform)

    # Define list of replacements which matches all_labels *IN ORDER*
    replacements = []
    replacement_dict = {}
    counter = 1
    for label in all_labels: # keep in order! this is important
        if label not in replacement_dict:
            replacement_dict[label] = str(counter)
            counter += 1
        replacements.append(replacement_dict[label])

    # Perform replacements in order
    transform_newmaps = re.sub('\:[0-9]+\]', 
        lambda match: (':' + replacements.pop(0) + ']'),
        transform)

    return transform_newmaps

def canonicalize_transform(transform):
    '''This function takes an atom-mapped SMARTS transform and
    converts it to a canonical form by, if nececssary, rearranging
    the order of reactant and product templates and reassigning
    atom maps.'''

    transform_reordered = '>>'.join([canonicalize_template(x) for x in transform.split('>>')])
    return reassign_atom_mapping(transform_reordered)

def process_an_example_doc(example_doc, test=False):
    '''Function for a worker to process one doc'''
    total_partialmapped = 0
    total_attempted = 0
    total_unmapped = 0


    if v: 
        print('##################################')
        print('###        RXN {}'.format(example_doc['_id']))
        print('##################################')

    if "nonmapped reaction" in example_doc['RX_SKW']:
        print('Unmapped reaction {}'.format(example_doc['_id']))
        total_unmapped += 1
        return

    # Check atom mapping
    try:
        prod_numatoms = example_doc['PROD_NUMATOMS']
        frac_prod_mapped = example_doc['FRAC_PROD_MAPPED']
        unmapped_atoms = prod_numatoms * (1 - frac_prod_mapped)
        if unmapped_atoms > MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS: 
            return
    except KeyError as e:
        return

    try:
        # Unpack
        reaction_smiles = str(example_doc[reaction_smiles_field])
        reactants, products = [mols_from_smiles_list(x) for x in 
                                    [mols.split('.') for mols in reaction_smiles.split('>>')]]
        if None in reactants: return 
        if None in products: return
        [AllChem.RemoveHs(mol) for mol in reactants + products]
        [Chem.SanitizeMol(mol) for mol in reactants + products]
    except Exception as e:
        # can't sanitize -> skip
        print(e)
        print('Could not load SMILES or sanitize')
        print('ID: {}'.format(example_doc['_id']))
        return
    
    try:
        ###
        ### Check product atom mapping to see if reagent contributes
        ###

        are_unmapped_product_atoms = False
        extra_reactant_fragment = ''
        for product in products:
            if sum([a.HasProp('molAtomMapNumber') for a in product.GetAtoms()]) < len(product.GetAtoms()):
                if v: print('Not all product atoms have atom mapping')
                if v: print('ID: {}'.format(example_doc['_id']))
                if v: print('REACTION: {}'.format(example_doc['RXN_SMILES']))
                are_unmapped_product_atoms = True
        if are_unmapped_product_atoms: # add fragment to template

            total_partialmapped += 1
            for product in products:
                # Get unmapped atoms
                unmapped_ids = [
                    a.GetIdx() for a in product.GetAtoms() if not a.HasProp('molAtomMapNumber')
                ]
                # Define new atom symbols for fragment with atom maps, generalizing fully
                atom_symbols = ['[{}]'.format(a.GetSymbol()) for a in product.GetAtoms()]
                # And bond symbols...
                bond_symbols = ['~' for b in product.GetBonds()]
                if unmapped_ids:
                    extra_reactant_fragment += \
                        AllChem.MolFragmentToSmiles(product, unmapped_ids, 
                        allHsExplicit = False, isomericSmiles = USE_STEREOCHEMISTRY, 
                        atomSymbols = atom_symbols, bondSymbols = bond_symbols) + '.'
            if extra_reactant_fragment:
                extra_reactant_fragment = extra_reactant_fragment[:-1]
                if v: print('    extra reactant fragment: {}'.format(extra_reactant_fragment))

            # Consolidate repeated fragments (stoichometry)
            extra_reactant_fragment = '.'.join(sorted(list(set(extra_reactant_fragment.split('.')))))
            #fragmatch = Chem.MolFromSmarts(extra_reactant_fragment) # no parentheses

        ###
        ### Do RX-level processing
        ###  


        if v: print(reaction_smiles)
        if None in reactants + products:
            print('Could not parse all molecules in reaction, skipping')
            print('ID: {}'.format(example_doc['_id']))
            return

        # Calculate changed atoms
        changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
        if err: 
            print('Could not get changed atoms')
            print('ID: {}'.format(example_doc['_id']))
            return
        if not changed_atom_tags:
            print('No atoms changed? (stereochemistry only?)')
            print('ID: {}'.format(example_doc['_id']))
            # print('Reaction SMILES: {}'.format(example_doc['RXN_SMILES']))
            return

        # Get fragments for reactants
        reactant_fragments = get_fragments_for_changed_atoms(reactants, changed_atom_tags, 
            radius = 1, expansion = [], category = 'reactants')
        # Get fragments for products 
        # (WITHOUT matching groups but WITH the addition of reactant fragments)
        product_fragments  = get_fragments_for_changed_atoms(products, changed_atom_tags, 
            radius = 0, expansion = expand_changed_atom_tags(changed_atom_tags, reactant_fragments),
            category = 'products')


        ###
        ### Put together and canonicalize (as best as possible)
        ###
        rxn_string = '{}>>{}'.format(reactant_fragments, product_fragments)
        rxn_canonical = canonicalize_transform(rxn_string)
        # print('Pre-resplit: {}'.format(rxn_canonical))
        # Change from inter-molecular to whatever-molecular
        rxn_canonical_split = rxn_canonical.split('>>')
        rxn_canonical = rxn_canonical_split[0][1:-1].replace(').(', '.') + \
            '>>' + rxn_canonical_split[1][1:-1].replace(').(', '.')
        
        reactants_string = rxn_canonical.split('>>')[0]
        products_string  = rxn_canonical.split('>>')[1]

        retro_canonical = products_string + '>>' + reactants_string

        # print('Original string: {}'.format(example_doc['RXN_SMILES']))
        # print('\nOverall retro transform: {}'.format(retro_canonical))
        # print('Extra necessary fragment: {}'.format(extra_reactant_fragment))

        # Load into RDKit
        rxn = AllChem.ReactionFromSmarts(retro_canonical)
        if rxn.Validate()[1] != 0: 
            print('Could not validate reaction successfully')
            print('ID: {}'.format(example_doc['_id']))
            print('retro_canonical: {}'.format(retro_canonical))
            print('original: {}'.format(example_doc['RXN_SMILES']))
            if v: raw_input('Pausing...')
            return

        ###
        ### Now look for specific RXDs.
        ###
        rxd_id_list = ['{}-{}'.format(example_doc['_id'], j) for j in range(1, int(example_doc['RX_NVAR']) + 1)]
        # All RXDs will have the same template - they'll just need different
        # references and contribute to the count differently

        # Insert - if it doesn't exist, Mongo will create the document
        # $inc, $addToSet, and $setOnInsert are necessary!
        if not test:
            TRANSFORM_DB.update_one(
                {'reaction_smarts': retro_canonical},
                {
                    '$inc': {
                        'count': len(rxd_id_list),
                    },
                    '$addToSet': {
                        'references': { '$each': rxd_id_list }
                    },
                    '$setOnInsert': {
                        'reaction_smarts': retro_canonical,
                        'necessary_reagent': extra_reactant_fragment,
                        'rxn_example': reaction_smiles,
                    }
                },
                True # upsert
            )
        if v: print('Added record for template {}'.format(retro_canonical))

        total_attempted += 1
    

    except KeyboardInterrupt:
        print('Interrupted')
        raise KeyboardInterrupt

    except Exception as e:
        print(e)
        if v: 
            print('skipping')
            #raw_input('Enter anything to continue')
        return


def main(N = 15, skip = 0, skip_id = 0):
    '''Read reactions'''
    global v

    # Define scoring variables
    total_attempted = 0 # total reactions simulated (excludes skipped)
    total_correct = 0 # actual products predicted
    total_precise = 0 # ONLY actual products predicted
    total_unmapped = 0
    total_partialmapped = 0
    total_nonreaction = 0

    # N = min([N, REACTION_DB.count({'RX_SKW': 'mapped reaction'})])

    def data_generator():
        data_generator.ctr = -1
        try: # to allow breaking
            # Look for entries
            for example_doc in tqdm(REACTION_DB.find({'_id': {'$gt': skip_id}, 'RX_SKW': 'mapped reaction'}).sort('_id', 1)[skip:]):
                data_generator.ctr += 1

                if data_generator.ctr % 1000 == 0:
                    print('count: {}'.format(data_generator.ctr))
                    print('id: {}'.format(example_doc['_id']))

                if example_doc['_id'] < skip_id: continue
                # if data_generator.ctr < skip: continue 
        
                # Are we done?
                if data_generator.ctr >= N:
                    data_generator.ctr -= 1
                    break

                yield example_doc

        except KeyboardInterrupt:
            print('Stopped early!')  

    
    generator = data_generator()
    from multiprocessing import Pool 
    pool = Pool(20)
    while True:
        print('Queueing up 1000 more examples from generator')
        res = pool.map(process_an_example_doc, itertools.islice(generator, 10000))
        if res:
            pass
        else:
            break
    pool.close()
    pool.join()

    # except Exception as e:
    #   print(e)

    total_examples = data_generator.ctr + 1
    print('...finished looking through {} reaction records'.format(min([N, ctr])))

    print('Unmapped reactions in {}/{} ({}%) cases'.format(total_unmapped,
        total_examples, total_unmapped * 100.0 / total_examples))
    print('Partially-mapped reactions in {}/{} ({}%) cases'.format(total_partialmapped,
        total_examples, total_partialmapped * 100.0 / total_examples))
    print('Non-reactions (stereochem only?) in {}/{} ({}%) cases'.format(total_nonreaction,
        total_examples, total_nonreaction * 100.0 / total_examples))
    print('Error-free parsing in {}/{} ({}%) cases'.format(total_attempted, 
        total_examples, total_attempted * 100.0 / total_examples))
    print('Correct product predictions in {}/{} ({}%) cases'.format(total_correct, 
        total_examples, total_correct * 100.0 / total_examples))
    print('Specific product predictions in {}/{} ({}%) cases'.format(total_precise, 
        total_examples, total_precise * 100.0 / total_examples))
    print('Created {} database entries'.format(TRANSFORM_DB.find().count()))

    return True


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', type = bool, default = False,
                        help = 'Verbose printing (incl. saving images); defaults to False')
    parser.add_argument('-n', '--num', type = int, default = 50,
                        help = 'Maximum number of records to examine; defaults to 50')
    parser.add_argument('--skip', type = int, default = 0,
                        help = 'Number of examples to skip, defaults 0')
    parser.add_argument('--skip_id', type = int, default = 0, 
                        help = 'IDs of reaction entries to skip to, defaults 0')
    args = parser.parse_args()

    v = args.v
    lg = RDLogger.logger()
    if not v: lg.setLevel(4)

    clear = raw_input('Do you want to clear the {} existing templates? '.format(TRANSFORM_DB.find().count()))
    if clear in ['y', 'Y', 'yes', '1', 'Yes']:
        result = TRANSFORM_DB.delete_many({})
        print('Cleared {} entries from collection'.format(result.deleted_count))
    TRANSFORM_DB.create_index([('reaction_smarts', 'hashed')])
    print('Added hashed index on reaction_smarts')

    # example_doc = REACTION_DB.find_one({'_id': 68309})
    # example_doc = REACTION_DB.find_one({'_id': 9004813})
    # process_an_example_doc(example_doc, test=True)
    # quit(1)

    main(N = args.num, skip = args.skip, skip_id = args.skip_id)
