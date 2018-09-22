from __future__ import print_function
import sys 
import os
import re

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir

from rdchiral.utils import vprint, PLEVEL
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rdchiral.chiral import template_atom_could_have_been_tetra, copy_chirality,\
    atom_chirality_matches
from rdchiral.clean import canonicalize_outcome_smiles, combine_enantiomers_into_racemic
from rdchiral.bonds import BondDirOpposite, restore_bond_stereo_to_sp2_atom

'''
This file contains the main functions for running reactions. 

An incomplete description of expected behavior is as follows:

(1) RDKit's native RunReactants is called on an achiral version of the molecule,
which has had all tetrahedral centers and bond directions stripped.

(2) For each outcome, we examine the correspondence between atoms in the
reactants and atoms in the reactant template for reasons to exclude the 
current outcome. The following conditions are checked:

    TETRAHEDRAL ATOMS
    (a) If a reactant atom is a tetrahedral center with specified chirality
        and the reactant template atom is NOT chiral but is defined in a way
        that it could have been specified, reject this outcome
    (b) If a reactant atom is a tetrahedral center with specified chirality
        and the reactant template atom is NOT chiral and is not defined in
        a way where it could have been (i.e., is generalized without spec.
        neighbors), then keep the match.
    (c) If a reactant atom is achiral but the reactant tempalte atom is chiral,
        the match is still allowed to happen. We might want to change this later
        or let it be an option.
    (d) If a reactant atom is a tetrahedral center with specified chirality
        and the reactant template also has its chirality specified, let the
        match happen if the chirality matches.


    DOUBLE BONDS
    (a) If a reactant double bond is defined with directionality specified and
        the reactant template is unspecified but COULD have been (i.e., 
        neighbors of sp2 carbons are specified), reject this outcome
    (b) If a reactant double bond is defined with directionality specified and
        the reactant template si unspecified but could NOT have been (in the
        case of generalization), allow the match to occur. This is what we
        default to when half the bond is specified, like in "C=C/O"
    note: reactants are checked for implicit bond stereo based on rings
    (c) If a reactant double bond has implicit cis due to ring membership, it is
        still allowed to match an unspecified template double bond. Might lead
        to some weird edge cases, but mostly makes sense.


(3) For each outcome, merge all products into a single molecule. During this
process, we check for bonds that are missing in the product. These are those
that were present in the reactants but were NOT matched in the reactant
template.

(4) For each outcome, examine product atoms to correct tetrahedral chirality.

(5) For each outcome, examine product double bonds to correct cis/trans-ness

'''

def rdchiralRunText(reaction_smarts, reactant_smiles, **kwargs):
    '''Run from SMARTS string and SMILES string. This is NOT recommended
    for library application, since initialization is pretty slow. You should
    separately initialize the template and molecules and call run()'''
    rxn = rdchiralReaction(reaction_smarts)
    reactants = rdchiralReactants(reactant_smiles)
    return rdchiralRun(rxn, reactants, **kwargs)

def rdchiralRun(rxn, reactants, keep_isotopes=False, combine_enantiomers=True):
    '''
    rxn = rdchiralReaction (rdkit reaction + auxilliary information)
    reactants = rdchiralReactants (rdkit mol + auxilliary information)

    note: there is a fair amount of initialization (assigning stereochem), most
    importantly assigning isotope numbers to the reactant atoms. It is 
    HIGHLY recommended to use the custom classes for initialization.
    '''

    

    ###############################################################################
    # Run naive RDKit on ACHIRAL version of molecules
    outcomes = rxn.rxn.RunReactants((reactants.reactants_achiral,))
    if PLEVEL >= (1): print('Using naive RunReactants, {} outcomes'.format(len(outcomes)))
    if not outcomes:
        return []

    ###############################################################################
    

    ###############################################################################
    # Initialize, now that there is at least one outcome

    final_outcomes = set()
    # We need to keep track of what map numbers 
    # (i.e., isotopes) correspond to which atoms
    # note: all reactant atoms must be mapped, so this is safe
    atoms_r = reactants.atoms_r

    # Copy reaction template so we can play around with isotopes
    template_r, template_p = rxn.template_r, rxn.template_p

    # Get molAtomMapNum->atom dictionary for tempalte reactants and products
    atoms_rt_map = rxn.atoms_rt_map
    atoms_pt_map = rxn.atoms_pt_map
    ###############################################################################


    for outcome in outcomes:
        ###############################################################################
        # Look for new atoms in products that were not in 
        # reactants (e.g., LGs for a retro reaction)
        if PLEVEL >= (2): print('Processing {}'.format(str([Chem.MolToSmiles(x, True) for x in outcome])))
        unmapped = 900
        for m in outcome:
            for a in m.GetAtoms():
                # Assign "map" number via isotope
                if not a.GetIsotope():
                    a.SetIsotope(unmapped)
                    unmapped += 1
        if PLEVEL >= 2: print('Added {} map numbers to product'.format(unmapped-900))
        ###############################################################################


        ###############################################################################
        # Check to see if reactants should not have been matched (based on chirality)

        # Define isotope -> reactant template atom map
        atoms_rt =  {a.GetIsotope(): atoms_rt_map[a.GetIntProp('old_mapno')] \
            for m in outcome for a in m.GetAtoms() if a.HasProp('old_mapno')}

        # Set isotopes of reactant template
        # note: this is okay to do within the loop, because ALL atoms must be matched
        # in the templates, so the isotopes will get overwritten every time
        [a.SetIsotope(i) for (i, a) in atoms_rt.items()]

        # Make sure each atom matches
        # note: this is a little weird because atom_chirality_matches takes three values,
        #       -1 (both tetra but opposite), 0 (not a match), and +1 (both tetra and match)
        #       and we only want to continue if they all equal -1 or all equal +1
        prev = None
        skip_outcome = False
        for match in (atom_chirality_matches(atoms_rt[i], atoms_r[i]) for i in atoms_rt):
            print(match)
            if match == 0: 
                if PLEVEL >= 2: print('Chirality violated! Should not have gotten this match')
                skip_outcome = True 
                break
            elif match == 2: # ambiguous case
                continue
            elif prev is None:
                prev = match
            elif match != prev:
                if PLEVEL >= 2: print('Part of the template matched reactant chirality, part is inverted! Should not match')
                skip_outcome = True 
                break
        if skip_outcome:
            continue      
        if PLEVEL >= 2: print('Chirality matches! Just checked with atom_chirality_matches')

        # Check bond chirality - iterate through reactant double bonds where
        # chirality is specified (or not). atoms defined by isotope
        skip_outcome = False
        for atoms, dirs, is_implicit in reactants.atoms_across_double_bonds:
            if all(i in atoms_rt for i in atoms):
                # All atoms definining chirality were matched to the reactant template
                # So, check if it is consistent with how the template is defined
                #...but /=/ should match \=\ since they are both trans...
                matched_atom_map_nums = tuple(atoms_rt[i].GetIntProp('molAtomMapNumber') for i in atoms)
                if matched_atom_map_nums not in rxn.required_rt_bond_defs:
                    continue # this can happen in ring openings, for example
                dirs_template = rxn.required_rt_bond_defs[matched_atom_map_nums]
                if dirs != dirs_template and \
                        (BondDirOpposite[dirs[0]], BondDirOpposite[dirs[1]]) != dirs_template and \
                        not (dirs_template == (BondDir.NONE, BondDir.NONE) and is_implicit):
                    if PLEVEL >= 5: print('Reactant bond chirality does not match template!')
                    if PLEVEL >= 5: print('Based on map numbers...')
                    if PLEVEL >= 5: print('  rct: {} -> {}'.format(matched_atom_map_nums, dirs))
                    if PLEVEL >= 5: print('  tmp: {} -> {}'.format(matched_atom_map_nums, dirs_template))
                    if PLEVEL >= 5: print('skipping this outcome, should not have matched...')
                    skip_outcome = True 
                    break
        if skip_outcome:
            continue

        ###############################################################################



        ###############################################################################
        # Convert product(s) to single product so that all 
        # reactions can be treated as pseudo-intramolecular
        # But! check for ring openings mistakenly split into multiple
        # This can be diagnosed by duplicate map numbers (i.e., SMILES)

        isotopes = [a.GetIsotope() for m in outcome for a in m.GetAtoms() if a.GetIsotope()]
        if len(isotopes) != len(set(isotopes)): # duplicate?
            if PLEVEL >= 1: print('Found duplicate isotopes in product - need to stitch')
            # need to do a fancy merge
            merged_mol = Chem.RWMol(outcome[0])
            merged_iso_to_id = {a.GetIsotope(): a.GetIdx() for a in outcome[0].GetAtoms() if a.GetIsotope()}
            for j in range(1, len(outcome)):
                new_mol = outcome[j]
                for a in new_mol.GetAtoms():
                    if a.GetIsotope() not in merged_iso_to_id:
                        merged_iso_to_id[a.GetIsotope()] = merged_mol.AddAtom(a)
                for b in new_mol.GetBonds():
                    bi = b.GetBeginAtom().GetIsotope()
                    bj = b.GetEndAtom().GetIsotope()
                    if PLEVEL >= 10: print('stitching bond between {} and {} in stich has chirality {}, {}'.format(
                        bi, bj, b.GetStereo(), b.GetBondDir()
                    ))
                    if not merged_mol.GetBondBetweenAtoms(
                            merged_iso_to_id[bi], merged_iso_to_id[bj]):
                        merged_mol.AddBond(merged_iso_to_id[bi],
                            merged_iso_to_id[bj], b.GetBondType())
                        merged_mol.GetBondBetweenAtoms(
                            merged_iso_to_id[bi], merged_iso_to_id[bj]
                        ).SetStereo(b.GetStereo())
                        merged_mol.GetBondBetweenAtoms(
                            merged_iso_to_id[bi], merged_iso_to_id[bj]
                        ).SetBondDir(b.GetBondDir())
            outcome = merged_mol.GetMol()
            if PLEVEL >= 1: print('Merged editable mol, converted back to real mol, {}'.format(Chem.MolToSmiles(outcome, True)))
        else:
            new_outcome = outcome[0]
            for j in range(1, len(outcome)):
                new_outcome = AllChem.CombineMols(new_outcome, outcome[j])
            outcome = new_outcome
        if PLEVEL >= 2: print('Converted all outcomes to single molecules')
        ###############################################################################




        ###############################################################################
        # Figure out which atoms were matched in the templates
        # atoms_rt and atoms_p will be outcome-specific.
        atoms_pt = {a.GetIsotope(): atoms_pt_map[a.GetIntProp('old_mapno')] \
            for a in outcome.GetAtoms() if a.HasProp('old_mapno')}
        atoms_p = {a.GetIsotope(): a for a in outcome.GetAtoms() if a.GetIsotope()}

        # Set isotopes of product template
        # note: this is okay to do within the loop, because ALL atoms must be matched
        # in the templates, so the isotopes will get overwritten every time
        # This makes it easier to check parity changes
        [a.SetIsotope(i) for (i, a) in atoms_pt.items()]
        ###############################################################################



        ###############################################################################
        # Check for missing bonds. These are bonds that are present in the reactants,
        # not specified in the reactant template, and not in the product. Accidental
        # fragmentation can occur for intramolecular ring openings
        missing_bonds = []
        for (i, j, b) in reactants.bonds_by_isotope:
            if i in atoms_p and j in atoms_p:
                # atoms from reactant bond show up in product
                if not outcome.GetBondBetweenAtoms(atoms_p[i].GetIdx(), atoms_p[j].GetIdx()):
                    #...but there is not a bond in the product between those atoms
                    if i not in atoms_rt or j not in atoms_rt or not template_r.GetBondBetweenAtoms(atoms_rt[i].GetIdx(), atoms_rt[j].GetIdx()):
                        # the reactant template did not specify a bond between those atoms (e.g., intentionally destroy)
                        missing_bonds.append((i, j, b))
        if missing_bonds:
            if PLEVEL >= 1: print('Product is missing non-reacted bonds that were present in reactants!')
            outcome = Chem.RWMol(outcome)
            rwmol_iso_to_id = {a.GetIsotope(): a.GetIdx() for a in outcome.GetAtoms() if a.GetIsotope()}
            for (i, j, b) in missing_bonds:
                outcome.AddBond(rwmol_iso_to_id[i], rwmol_iso_to_id[j])
                new_b = outcome.GetBondBetweenAtoms(rwmol_iso_to_id[i], rwmol_iso_to_id[j])
                new_b.SetBondType(b.GetBondType())
                new_b.SetBondDir(b.GetBondDir())
                new_b.SetIsAromatic(b.GetIsAromatic())
            outcome = outcome.GetMol()
        else:
            if PLEVEL >= 3: print('No missing bonds')
        ###############################################################################


        # Now that we've fixed any bonds, connectivity is set. This is a good time
        # to udpate the property cache, since all that is left is fixing atom/bond
        # stereochemistry.
        try:
            Chem.SanitizeMol(outcome)
            outcome.UpdatePropertyCache()
        except ValueError as e: 
            if PLEVEL >= 1: print('{}, {}'.format(Chem.MolToSmiles(outcome, True), e))
            continue


        ###############################################################################
        # Correct tetra chirality in the outcome
        tetra_copied_from_reactants = []
        for a in outcome.GetAtoms():
            # Participants in reaction core (from reactants) will have old_mapno
            # Spectators present in reactants will have react_atom_idx
            # ...so new atoms will have neither!
            if not a.HasProp('old_mapno'):
                # Not part of the reactants template
                
                if not a.HasProp('react_atom_idx'):
                    # Atoms only appear in product template - their chirality
                    # should be properly instantiated by RDKit...hopefully...
                    if PLEVEL >= 4: print('Atom {} created by product template, should have right chirality'.format(a.GetIsotope()))
                
                else:
                    if PLEVEL >= 4: print('Atom {} outside of template, copy chirality from reactants'.format(a.GetIsotope()))
                    copy_chirality(atoms_r[a.GetIsotope()], a)
                    if a.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
                        tetra_copied_from_reactants.append(a)

            else:
                # Part of reactants and reaction core
                
                if template_atom_could_have_been_tetra(atoms_rt[a.GetIsotope()]):
                    if PLEVEL >= 3: print('Atom {} was in rct template (could have been tetra)'.format(a.GetIsotope()))
                    
                    if template_atom_could_have_been_tetra(atoms_pt[a.GetIsotope()]):
                        if PLEVEL >= 3: print('Atom {} in product template could have been tetra, too'.format(a.GetIsotope()))
                        
                        # Was the product template specified?
                        
                        if atoms_pt[a.GetIsotope()].GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                            # No, leave unspecified in product
                            if PLEVEL >= 3: print('...but it is not specified in product, so destroy chirality')
                            a.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
                        
                        else:
                            # Yes
                            if PLEVEL >= 3: print('...and product is specified')
                            
                            # Was the reactant template specified?
                            
                            if atoms_rt[a.GetIsotope()].GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                                # No, so the reaction introduced chirality
                                if PLEVEL >= 3: print('...but reactant template was not, so copy from product template')
                                copy_chirality(atoms_pt[a.GetIsotope()], a)
                            
                            else:
                                # Yes, so we need to check if chirality should be preserved or inverted
                                if PLEVEL >= 3: print('...and reactant template was, too! copy from reactants')
                                copy_chirality(atoms_r[a.GetIsotope()], a)
                                if atom_chirality_matches(atoms_pt[a.GetIsotope()], atoms_rt[a.GetIsotope()]) == -1:
                                    if PLEVEL >= 3: print('but! reactant template and product template have opposite stereochem, so invert')
                                    a.InvertChirality()
                    
                    else:
                        # Reactant template chiral, product template not - the
                        # reaction is supposed to destroy chirality, so leave
                        # unspecified
                        if PLEVEL >= 3: print('If reactant template could have been ' +
                            'chiral, but the product template could not, then we dont need ' +
                            'to worry about specifying product atom chirality')

                else:
                    if PLEVEL >= 3: print('Atom {} could not have been chiral in reactant template'.format(a.GetIsotope()))
                    
                    if not template_atom_could_have_been_tetra(atoms_pt[a.GetIsotope()]):
                        if PLEVEL >= 3: print('Atom {} also could not have been chiral in product template', a.GetIsotope())
                        if PLEVEL >= 3: print('...so, copy chirality from reactant instead')
                        copy_chirality(atoms_r[a.GetIsotope()], a)
                        if a.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
                            tetra_copied_from_reactants.append(a)
                    
                    else:
                        if PLEVEL >= 3: print('Atom could/does have product template chirality!'.format(a.GetIsotope()))
                        if PLEVEL >= 3: print('...so, copy chirality from product template')
                        copy_chirality(atoms_pt[a.GetIsotope()], a)
                    
            if PLEVEL >= 3: print('New chiral tag {}'.format(a.GetChiralTag()))
        if skip_outcome:
            if PLEVEL >= 2: print('Skipping this outcome - chirality broken?')
            continue
        if PLEVEL >= 2: print('After attempting to re-introduce chirality, outcome = {}'.format(Chem.MolToSmiles(outcome, True)))
        ###############################################################################


        ###############################################################################
        # Correct bond directionality in the outcome
        for b in outcome.GetBonds():
            if b.GetBondType() != BondType.DOUBLE:
                continue

            # Ring double bonds do not need to be touched(?)
            if b.IsInRing():
                continue
            
            ba = b.GetBeginAtom()
            bb = b.GetEndAtom()

            # Is it possible at all to specify this bond?
            if ba.GetDegree() == 1 or bb.GetDegree() == 1:
                continue

            if PLEVEL >= 5: print('Looking at outcome bond {}={}'.format(ba.GetIsotope(), bb.GetIsotope()))

            if ba.HasProp('old_mapno') and bb.HasProp('old_mapno'):
                # Need to rely on templates for bond chirality, both atoms were
                # in the reactant template 
                if PLEVEL >= 5: print('Both atoms in this double bond were in the reactant template')
                if (ba.GetIntProp('old_mapno'), bb.GetIntProp('old_mapno')) in \
                        rxn.required_bond_defs_coreatoms:   
                    if PLEVEL >= 5: print('and reactant template *could* have specified the chirality!')
                    if PLEVEL >= 5: print('..product should be property instantiated')
                    continue
                if PLEVEL >= 5: print('But it was impossible to have specified chirality (e.g., aux C=C for context)')

            elif not ba.HasProp('react_atom_idx') and not bb.HasProp('react_atom_idx'):
                # The atoms were both created by the product template, so any bond
                # stereochemistry should have been instantiated by the product template
                # already...hopefully...otherwise it isn't specific enough?
                continue 

            # Need to copy from reactants, this double bond was simply carried over,
            # *although* one of the atoms could have reacted and been an auxilliary
            # atom in the reaction, e.g., C/C=C(/CO)>>C/C=C(/C[Br])
            if PLEVEL >= 5: print('Restoring cis/trans character of bond {}={} from reactants'.format(
                ba.GetIsotope(), bb.GetIsotope()))
            
            # Start with setting the BeginAtom
            begin_atom_specified = restore_bond_stereo_to_sp2_atom(ba, reactants.bond_dirs_by_isotope)
            
            if not begin_atom_specified:
                # don't bother setting other side of bond, since we won't be able to
                # fully specify this bond as cis/trans
                continue 
            
            # Look at other side of the bond now, the EndAtom
            end_atom_specified = restore_bond_stereo_to_sp2_atom(bb, reactants.bond_dirs_by_isotope)
            if not end_atom_specified:
                # note: this can happen if C=C/C-N turns into C=C/C=N 
                if PLEVEL >= 1:
                    print(reactants.bond_dirs_by_isotope)
                    print(ba.GetIsotope())
                    print(bb.GetIsotope())
                    print(Chem.MolToSmiles(reactants.reactants, True))
                    print(Chem.MolToSmiles(outcome, True))
                    print('Uh oh, looks like bond direction is only specified for half of this bond?')

        ###############################################################################

        if not keep_isotopes:
            [a.SetIsotope(0) for a in outcome.GetAtoms()]
            

        # Now, check to see if we have destroyed chirality by removing isotopes
        # this occurs when chirality was not actually possible (e.g., due to
        # symmetry) but we had assigned a tetrahedral center originating
        # from the reactants.
        #    ex: SMILES C(=O)1C[C@H](Cl)CCC1
        #        SMARTS [C:1]-[C;H0;D3;+0:2](-[C:3])=[O;H0;D1;+0]>>[C:1]-[CH2;D2;+0:2]-[C:3]
        skip_outcome = False
        if len(tetra_copied_from_reactants) > 0:
            Chem.AssignStereochemistry(outcome, cleanIt=True, force=True)
            for a in tetra_copied_from_reactants:
                if a.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                    if PLEVEL >= 2: print('Auxiliary reactant atom was chiral, now is broken -> skip outcome')
                    skip_outcome = True 
                    break 
        if skip_outcome:
            continue


        smiles = Chem.MolToSmiles(outcome, True)
        smiles_new = canonicalize_outcome_smiles(smiles)
        if smiles_new is None:
            continue

        final_outcomes.add(smiles_new)

    ###############################################################################
    # One last fix for consolidating multiple stereospecified products...
    if combine_enantiomers:
        final_outcomes = combine_enantiomers_into_racemic(final_outcomes)
    ###############################################################################

    return list(final_outcomes)


if __name__ == '__main__':
    reaction_smarts = '[C:1][OH:2]>>[C:1][O:2][C]'
    reactant_smiles = 'CC(=O)OCCCO'
    outcomes = rdchiralRunText(reaction_smarts, reactant_smiles)
    print(outcomes)

    PLEVEL = 5
    reaction_smarts = '[#7:1]-[C:2](=[O;D1;H0:3])-[CH;@;D3;+0:4]1-[CH2;D2;+0:5]-[CH;@;D3;+0:9](-[C:7](-[#8:6])=[O;D1;H0:8])-[NH;D2;+0:10]-[CH;@;D3;+0:11]-1-[c:12]>>[#7:1]-[C:2](=[O;D1;H0:3])-[CH;D2;+0:4]=[CH2;D1;+0:5].[#8:6]-[C:7](=[O;D1;H0:8])-[CH2;D2;+0:9]/[N;H0;D2;+0:10]=[CH;D2;+0:11]/[c:12]'
    reactant_smiles = 'CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1'
    outcomes = rdchiralRunText(reaction_smarts, reactant_smiles)
    print(outcomes)