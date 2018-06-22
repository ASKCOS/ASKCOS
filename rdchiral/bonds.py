import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir, BondStereo
from rdchiral.utils import vprint, PLEVEL

BondDirOpposite = {AllChem.BondDir.ENDUPRIGHT: AllChem.BondDir.ENDDOWNRIGHT,
                   AllChem.BondDir.ENDDOWNRIGHT: AllChem.BondDir.ENDUPRIGHT,
                   AllChem.BondDir.NONE: AllChem.BondDir.NONE}
BondDirLabel    = {AllChem.BondDir.ENDUPRIGHT: '\\',
                   AllChem.BondDir.ENDDOWNRIGHT: '/'}

def bond_dirs_by_mapnum(mol):
    bond_dirs_by_mapnum = {}
    for b in mol.GetBonds():
        i = None; j = None
        if b.GetBeginAtom().HasProp('molAtomMapNumber'):
            i = b.GetBeginAtom().GetIntProp('molAtomMapNumber')
        if b.GetEndAtom().HasProp('molAtomMapNumber'):
            j = b.GetEndAtom().GetIntProp('molAtomMapNumber')
        if i is None or j is None or b.GetBondDir() == BondDir.NONE:
            continue
        bond_dirs_by_mapnum[(i, j)] = b.GetBondDir()
        bond_dirs_by_mapnum[(j, i)] = BondDirOpposite[b.GetBondDir()]
    return bond_dirs_by_mapnum

def enumerate_possible_cistrans_defs(template_r, \
        labeling_func=lambda a: a.GetIntProp('molAtomMapNumber')):
    '''
    This function is meant to take a reactant template and fully enumerate
    all the ways in which different double-bonds can have their cis/trans
    chirality specified (based on labeling_func). This is necessary because
    double-bond chirality cannot be specified using cis/trans (global properties)
    but must be done using ENDUPRIGHT and ENDDOWNRIGHT for the attached single
    bonds (local properties). Now, the next issue is that on each side of the
    double bond, only one of the single bond directions must be specified, and
    that direction can be using either atom order. e.g.,

    A1         B1
       \      /
         C = C
       /      \
    A2         B2

    Can be specified by:
    A1-C is an ENDDOWNRIGHT, C-B1 is an ENDUPRIGHT
    A1-C is an ENDDOWNRIGHT, C-B2 is an ENDDOWNRIGHT
    A1-C is an ENDDOWNRIGHT, B1-C is an ENDDOWNRIGHT
    A1-C is an ENDDOWNRIGHT, B2-C is an ENDUPRIGHT
    ...and twelve more definitions using different A1/A2 specs.

    ALSO - we can think about horizontally reflecting this bond entirely,
    which gets us even more definitions.

    So, the point of this function is to fully enumerate *all* of the ways
    in which chirality could have been specified. That way, we can take a
    reactant atom and check if its chirality is within the list of acceptable
    definitions to determine if a match was made.

    Gross.

    The way we do this is by first defining the *local* chirality of a double
    bond, which weights side chains based purely on the unique isotope numbering.
    Once we have a local cis/trans definition for a double bond, we can enumerate
    the sixteen possible ways that a reactant could match it.

    #TODO: figure out implicit cis from ring membership
    '''

    required_bond_defs = {}
    required_bond_defs_coreatoms = set()

    if PLEVEL >= 10: print('Looking at initializing template frag')
    for b in template_r.GetBonds():
        if b.GetBondType() != BondType.DOUBLE:
            continue

        # Define begin and end atoms of the double bond
        ba = b.GetBeginAtom()
        bb = b.GetEndAtom()

        # Now check if it is even possible to specify
        if ba.GetDegree() == 1 or bb.GetDegree() == 1:
            continue

        ba_label = labeling_func(ba)
        bb_label = labeling_func(bb)
            
        if PLEVEL >= 10: print('Found a double bond with potential cis/trans (based on degree)')
        if PLEVEL >= 10: print('{} {} {}'.format(ba_label,
                               b.GetSmarts(),
                               bb_label))

        # Save core atoms so we know that cis/trans was POSSIBLE to specify
        required_bond_defs_coreatoms.add((ba_label, bb_label))
        required_bond_defs_coreatoms.add((bb_label, ba_label))
            
        # Define heaviest isotope neighbor for each atom, excluding the other side of the double bond
        ba_neighbor_labels = [labeling_func(a) for a in ba.GetNeighbors()]
        ba_neighbor_labels.remove(bb_label) # remove other side of =
        ba_neighbor_labels_max = max(ba_neighbor_labels)
        bb_neighbor_labels = [labeling_func(a) for a in bb.GetNeighbors()]
        bb_neighbor_labels.remove(ba_label) # remove other side of =
        bb_neighbor_labels_max = max(bb_neighbor_labels)

        # The direction of the bond being observed might need to be flipped,
        # based on
        #     (a) if it is the heaviest atom on this side, and 
        #     (b) if the begin/end atoms for the directional bond are 
        #         in the wrong order (i.e., if the double-bonded atom 
        #         is the begin atom)
        front_spec = None
        back_spec = None
        for bab in ba.GetBonds():
            if bab.GetBondDir() != BondDir.NONE:
                if labeling_func(bab.GetBeginAtom()) == ba_label:
                    # Bond is in wrong order - flip
                    if labeling_func(bab.GetEndAtom()) != ba_neighbor_labels_max:
                        # Defined atom is not the heaviest - flip
                        front_spec = bab.GetBondDir()
                        break
                    front_spec = BondDirOpposite[bab.GetBondDir()]
                    break
                if labeling_func(bab.GetBeginAtom()) != ba_neighbor_labels_max:
                    # Defined atom is not heaviest
                    front_spec = BondDirOpposite[bab.GetBondDir()]
                    break
                front_spec = bab.GetBondDir()
                break
        if front_spec is None:
            if PLEVEL >= 10: print('Chirality not specified at front end of the bond!')
        else:
            if PLEVEL >= 10: print('Front specification: {}'.format(front_spec))
            
            for bbb in bb.GetBonds():
                if bbb.GetBondDir() != BondDir.NONE:
                    # For the "back" specification, the double-bonded atom *should* be the BeginAtom
                    if labeling_func(bbb.GetEndAtom()) == bb_label:
                        # Bond is in wrong order - flip
                        if labeling_func(bbb.GetBeginAtom()) != bb_neighbor_labels_max:
                            # Defined atom is not the heaviest - flip
                            back_spec = bbb.GetBondDir()
                            break
                        back_spec = BondDirOpposite[bbb.GetBondDir()]
                        break
                    if labeling_func(bbb.GetEndAtom()) != bb_neighbor_labels_max:
                        # Defined atom is not heaviest - flip
                        back_spec = BondDirOpposite[bbb.GetBondDir()]
                        break
                    back_spec = bbb.GetBondDir()
                    break
        if back_spec is None:
            if PLEVEL >= 10: print('Chirality not specified at back end of the bond!')
        else:
            if PLEVEL >= 10: print('Back specification: {}'.format(back_spec))

        # Is this an overall unspecified bond? Put it in the dictionary anyway, 
        # so there is something to match
        if front_spec is None or back_spec is None:
            # Create a definition over this bond so that reactant MUST be unspecified, too
            for start_atom in ba_neighbor_labels:
                for end_atom in bb_neighbor_labels:
                    required_bond_defs[(start_atom, ba_label, bb_label, end_atom)] = (BondDir.NONE, BondDir.NONE)
                    required_bond_defs[(ba_label, start_atom, bb_label, end_atom)] = (BondDir.NONE, BondDir.NONE)
                    required_bond_defs[(start_atom, ba_label, end_atom, bb_label)] = (BondDir.NONE, BondDir.NONE)
                    required_bond_defs[(ba_label, start_atom, end_atom, bb_label)] = (BondDir.NONE, BondDir.NONE)
                    required_bond_defs[(bb_label, end_atom, start_atom, ba_label)] = (BondDir.NONE, BondDir.NONE)
                    required_bond_defs[(end_atom, bb_label, start_atom, ba_label)] = (BondDir.NONE, BondDir.NONE)
                    required_bond_defs[(bb_label, end_atom, ba_label, start_atom)] = (BondDir.NONE, BondDir.NONE)
                    required_bond_defs[(end_atom, bb_label, ba_label, start_atom)] = (BondDir.NONE, BondDir.NONE)
            continue
        
        if front_spec == back_spec:
            if PLEVEL >= 10: print('-> locally TRANS')
            b.SetProp('localChirality', 'trans')
        else:
            if PLEVEL >= 10: print('--> locally CIS')
            b.SetProp('localChirality', 'cis')

        possible_defs = {}
        for start_atom in ba_neighbor_labels:
            for end_atom in bb_neighbor_labels:
                needs_inversion = (start_atom != ba_neighbor_labels_max) != \
                    (end_atom != bb_neighbor_labels_max)
                for start_atom_dir in [BondDir.ENDUPRIGHT, BondDir.ENDDOWNRIGHT]:
                    # When locally trans, BondDir of start shold be same as end, 
                    # unless we need inversion
                    if (front_spec != back_spec) != needs_inversion: 
                        # locally cis and does not need inversion (True, False)
                        # or locally trans and does need inversion (False, True)
                        end_atom_dir = BondDirOpposite[start_atom_dir]
                    else:
                        # locally cis and does need inversion (True, True)
                        # or locally trans and does not need inversion (False, False)
                        end_atom_dir = start_atom_dir

                    # Enumerate all combinations of atom orders...
                    possible_defs[(start_atom, ba_label, bb_label, end_atom)] = (start_atom_dir, end_atom_dir)
                    possible_defs[(ba_label, start_atom, bb_label, end_atom)] = (BondDirOpposite[start_atom_dir], end_atom_dir)
                    possible_defs[(start_atom, ba_label, end_atom, bb_label)] = (start_atom_dir, BondDirOpposite[end_atom_dir])
                    possible_defs[(ba_label, start_atom, end_atom, bb_label)] = (BondDirOpposite[start_atom_dir], BondDirOpposite[end_atom_dir])

                    possible_defs[(bb_label, end_atom, start_atom, ba_label)] = (end_atom_dir, start_atom_dir)
                    possible_defs[(bb_label, end_atom, ba_label, start_atom)] = (end_atom_dir, BondDirOpposite[start_atom_dir])
                    possible_defs[(end_atom, bb_label, start_atom, ba_label)] = (BondDirOpposite[end_atom_dir], start_atom_dir)
                    possible_defs[(end_atom, bb_label, ba_label, start_atom)] = (BondDirOpposite[end_atom_dir], BondDirOpposite[start_atom_dir])

        # Save to the definition of this bond (in either direction)
        required_bond_defs.update(possible_defs)
        
    if PLEVEL >= 10: print('All bond specs for this template:' )
    if PLEVEL >= 10: print(str([(k, v) for (k, v) in required_bond_defs.items()]))
    return required_bond_defs, required_bond_defs_coreatoms

def get_atoms_across_double_bonds(mol, labeling_func=lambda a:a.GetIsotope()):
    '''
    This function takes a molecule and returns a list of cis/trans specifications
    according to the following:

    (isotopes, dirs)

    where atoms = (a1, a2, a3, a4) and dirs = (d1, d2)
    and (a1, a2) defines the ENDUPRIGHT/ENDDOWNRIGHT direction of the "front"
    of the bond using d1, and (a3, a4) defines the direction of the "back" of 
    the bond using d2.

    This is used to initialize reactants with a SINGLE definition constraining
    the chirality. Templates have their chirality fully enumerated, so we can
    match this specific definition to the full set of possible definitions
    when determining if a match should be made.

    NOTE: the atom isotopes are returned. This is so we can later use them
    to get the old_mapno property from the corresponding product atom, which is
    an outcome-specific assignment

    We also include implicit chirality here based on ring membership, but keep
    track of that separately
    '''
    atoms_across_double_bonds = []
    atomrings = None

    for b in mol.GetBonds():
        if b.GetBondType() != BondType.DOUBLE:
            continue

        # Define begin and end atoms of the double bond
        ba = b.GetBeginAtom()
        bb = b.GetEndAtom()

        # Now check if it is even possible to specify
        if ba.GetDegree() == 1 or bb.GetDegree() == 1:
            continue

        ba_label = labeling_func(ba)
        bb_label = labeling_func(bb)
            
        if PLEVEL >= 5: print('Found a double bond with potential cis/trans (based on degree)')
        if PLEVEL >= 5: print('{} {} {}'.format(ba_label,
                               b.GetSmarts(),
                               bb_label))
        
        # Try to specify front and back direction separately
        front_isotopes = None 
        front_dir = None 
        back_isotopes = None 
        back_dir = None
        is_implicit = False 
        bab = None; bbb = None;
        for bab in (z for z in ba.GetBonds() if z.GetBondType() != BondType.DOUBLE):
            if bab.GetBondDir() != BondDir.NONE:
                front_isotopes = (labeling_func(bab.GetBeginAtom()), labeling_func(bab.GetEndAtom()))
                front_dir = bab.GetBondDir()
                break 
        for bbb in (z for z in bb.GetBonds() if z.GetBondType() != BondType.DOUBLE):
            if bbb.GetBondDir() != BondDir.NONE:
                back_isotopes = (labeling_func(bbb.GetBeginAtom()), labeling_func(bbb.GetEndAtom()))
                back_dir = bbb.GetBondDir()
                break 

        # If impossible to spec, just continue
        if (bab is None or bbb is None):
            continue

        # Did we actually get a specification out?
        if (front_dir is None or back_dir is None):

            if b.IsInRing(): 
                # Implicit cis! Now to figure out right definitions...
                if atomrings is None:
                    atomrings = mol.GetRingInfo().AtomRings() # tuple of tuples of atomIdx
                for atomring in atomrings:
                    if ba.GetIdx() in atomring and bb.GetIdx() in atomring:
                        front_isotopes = (labeling_func(bab.GetOtherAtom(ba)), ba_label)
                        back_isotopes = (bb_label, labeling_func(bbb.GetOtherAtom(bb)))
                        if (bab.GetOtherAtomIdx(ba.GetIdx()) in atomring) != \
                                (bbb.GetOtherAtomIdx(bb.GetIdx()) in atomring):
                            # one of these atoms are in the ring, one is outside -> trans
                            if PLEVEL >= 10: print('Implicit trans found')
                            front_dir = BondDir.ENDUPRIGHT
                            back_dir = BondDir.ENDUPRIGHT
                        else:
                            if PLEVEL >= 10: print('Implicit cis found')
                            front_dir = BondDir.ENDUPRIGHT 
                            back_dir = BondDir.ENDDOWNRIGHT
                        is_implicit = True 
                        break

            else:
                # Okay no, actually unspecified
                # Specify direction as BondDir.NONE using whatever bab and bbb were at the end fo the loop
                # note: this is why we use "for bab in ___generator___", so that we know the current
                #       value of bab and bbb correspond to a single bond we can def. by
                front_isotopes = (labeling_func(bab.GetBeginAtom()), labeling_func(bab.GetEndAtom()))
                front_dir = BondDir.NONE
                back_isotopes = (labeling_func(bbb.GetBeginAtom()), labeling_func(bbb.GetEndAtom()))
                back_dir = BondDir.NONE

        # Save this (a1, a2, a3, a4) -> (d1, d2) spec
        atoms_across_double_bonds.append(
            (
                front_isotopes + back_isotopes,
                (front_dir, back_dir),
                is_implicit,
            )
        )

    return atoms_across_double_bonds

def restore_bond_stereo_to_sp2_atom(a, bond_dirs_by_isotope):
    '''Copy over single-bond directions (ENDUPRIGHT, ENDDOWNRIGHT) to 
    the single bonds attached to some double-bonded atom, a

    a - atom with a double bond
    bond_dirs_by_isotope - dictionary of (begin_isotope, end_isotope): bond_dir
        that defines if a bond should be ENDUPRIGHT or ENDDOWNRIGHT. The reverse
        key is also included with the reverse bond direction. If the source
        molecule did not have a specified chirality at this double bond, then
        the isotope tuples will be missing from the dict

    In some cases, like C=C/O>>C=C/Br, we should assume that stereochem was
    preserved, even though isotopes won't match. There might be some reactions
    where the chirality is inverted (like C=C/O >> C=C\Br), but let's not
    worry about those for now...

    Returns True if a bond direction was copied'''

    for bond_to_spec in a.GetBonds():
        if (bond_to_spec.GetOtherAtom(a).GetIsotope(), a.GetIsotope()) in bond_dirs_by_isotope:
            bond_to_spec.SetBondDir(
                bond_dirs_by_isotope[
                    (bond_to_spec.GetBeginAtom().GetIsotope(),
                     bond_to_spec.GetEndAtom().GetIsotope())
                ]
            )
            if PLEVEL >= 2: print('Tried to copy bond direction b/w {} and {}'.format(
                bond_to_spec.GetBeginAtom().GetIsotope(),
                bond_to_spec.GetEndAtom().GetIsotope()
            ))
            return True
    
    # Weird case, like C=C/O >> C=C/Br
    if PLEVEL >= 2: print('Bond stereo could not be restored to sp2 atom, missing the branch that was used to define before...')
    
    if a.GetDegree() == 2:
        # Either the branch used to define was replaced with H (deg 3 -> deg 2)
        # or the branch used to define was reacted (deg 2 -> deg 2)
        for bond_to_spec in a.GetBonds():
            if bond_to_spec.GetBondType() == BondType.DOUBLE:
                continue
            if not bond_to_spec.GetOtherAtom(a).HasProp('old_mapno'): 
                # new atom, deg2->deg2, assume direction preserved
                if PLEVEL >= 5: print('Only single-bond attachment to atom {} is new, try to reproduce chirality'.format(a.GetIsotope()))
                needs_inversion = False 
            else:
                # old atom, just was not used in chirality definition - set opposite
                if PLEVEL >= 5: print('Only single-bond attachment to atom {} is old, try to reproduce chirality'.format(a.GetIsotope()))
                needs_inversion = True

            for (i, j), bond_dir in bond_dirs_by_isotope.items():
                if bond_dir != BondDir.NONE:
                    if i == bond_to_spec.GetBeginAtom().GetIsotope():
                        if needs_inversion:
                            bond_to_spec.SetBondDir(BondDirOpposite[bond_dir])
                        else:
                            bond_to_spec.SetBondDir(bond_dir)
                        return True 

    if a.GetDegree() == 3:
        # If we lost the branch defining stereochem, it must have been replaced
        for bond_to_spec in a.GetBonds():
            if bond_to_spec.GetBondType() == BondType.DOUBLE:
                continue
            oa = bond_to_spec.GetOtherAtom(a)
            if oa.HasProp('old_mapno') or oa.HasProp('react_atom_idx'):
                # looking at an old atom, which should have opposite direction as removed atom
                needs_inversion = True 
            else:
                # looking at a new atom, assume same as removed atom
                needs_inversion = False 

            for (i, j), bond_dir in bond_dirs_by_isotope.items():
                if bond_dir != BondDir.NONE:
                    if i == bond_to_spec.GetBeginAtom().GetIsotope():
                        if needs_inversion:
                            bond_to_spec.SetBondDir(BondDirOpposite[bond_dir])
                        else:
                            bond_to_spec.SetBondDir(bond_dir)
                        return True 

    return False