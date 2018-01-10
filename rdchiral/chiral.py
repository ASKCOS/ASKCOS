from __future__ import print_function
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir

from rdchiral.utils import vprint, parity4, PLEVEL

def template_atom_could_have_been_tetra(a, strip_if_spec=False, cache=True):
    '''
    Could this atom have been a tetrahedral center?
    If yes, template atom is considered achiral and will not match a chiral rct
    If no, the tempalte atom is auxilliary and we should not use it to remove
    a matched reaction. For example, a fully-generalized terminal [C:1] 
    '''

    if a.HasProp('tetra_possible'):
        return a.GetBoolProp('tetra_possible')
    if a.GetDegree() < 3 or (a.GetDegree() == 3 and 'H' not in a.GetSmarts()):
        if cache:
            a.SetBoolProp('tetra_possible', False)
        if strip_if_spec: # Clear chiral tag in case improperly set
            a.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
        return False 
    if cache:
        a.SetBoolProp('tetra_possible', True)
    return True 



def copy_chirality(a_src, a_new):

    # Not possible to be a tetrahedral center anymore?
    if a_new.GetDegree() < 3:
        return 
    if a_new.GetDegree() == 3 and \
            any(b.GetBondType() != BondType.SINGLE for b in a_new.GetBonds()):
        return

    if PLEVEL >= 3: print('For isotope {}, copying src {} chirality tag to new'.format(
        a_src.GetIsotope(), a_src.GetChiralTag()))
    a_new.SetChiralTag(a_src.GetChiralTag())
    
    if atom_chirality_matches(a_src, a_new) == -1:
        if PLEVEL >= 3: print('For isotope {}, inverting chirality'.format(a_new.GetIsotope()))
        a_new.InvertChirality()

def atom_chirality_matches(a_tmp, a_mol):
    '''
    Checks for consistency in chirality between a template atom and a molecule atom.

    Also checks to see if chirality needs to be inverted in copy_chirality

    Returns +1 if it is a match and there is no need for inversion (or ambiguous)
    Returns -1 if it is a match but they are the opposite
    Returns 0 if an explicit NOT match
    Returns 2 if ambiguous or achiral-achiral
    '''
    if a_mol.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
        if a_tmp.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
            if PLEVEL >= 3: print('atom {} is achiral & achiral -> match'.format(a_mol.GetIsotope()))
            return 2 # achiral template, achiral molecule -> match
        # What if the template was chiral, but the reactant isn't just due to symmetry?
        if not a_mol.HasProp('_ChiralityPossible'):
            # It's okay to make a match, as long as the product is achiral (even
            # though the product template will try to impose chirality)
            if PLEVEL >= 3: print('atom {} is specified in template, but cant possibly be chiral in mol'.format(a_mol.GetIsotope()))
            return 2

        # TODO: figure out if we want this behavior - should a chiral template
        # be applied to an achiral molecule? For the retro case, if we have
        # a retro reaction that requires a specific stereochem, return False;
        # however, there will be many cases where the reaction would probably work
        if PLEVEL >= 3: print('atom {} is achiral in mol, but specified in template'.format(a_mol.GetIsotope()))
        return 0
    if a_tmp.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
        if PLEVEL >= 3: print('Reactant {} atom chiral, rtemplate achiral...'.format(a_tmp.GetIsotope()))
        if template_atom_could_have_been_tetra(a_tmp):
            if PLEVEL >= 3: print('...and that atom could have had its chirality specified! no_match')
            return 0
        if PLEVEL >= 3: print('...but the rtemplate atom could not have had chirality specified, match anyway')
        return 2

    isotopes_tmp = [a.GetIsotope() for a in a_tmp.GetNeighbors()]
    isotopes_mol = [a.GetIsotope() for a in a_mol.GetNeighbors()]

    # When there are fewer than 3 heavy neighbors, chirality is ambiguous...
    if len(isotopes_tmp) < 3 or len(isotopes_mol) < 3:
        return 2

    # Degree of 3 -> remaining atom is a hydrogen, add to list
    if len(isotopes_tmp) < 4:
        isotopes_tmp.append(-1) # H
    if len(isotopes_mol) < 4:
        isotopes_mol.append(-1) # H

    try:
        if PLEVEL >= 10: print(str(isotopes_tmp))
        if PLEVEL >= 10: print(str(isotopes_mol))
        if PLEVEL >= 10: print(str(a_tmp.GetChiralTag()))
        if PLEVEL >= 10: print(str(a_mol.GetChiralTag()))
        only_in_src = [i for i in isotopes_tmp if i not in isotopes_mol][::-1] # reverse for popping
        only_in_mol = [i for i in isotopes_mol if i not in isotopes_tmp]
        if len(only_in_src) <= 1 and len(only_in_mol) <= 1:
            tmp_parity = parity4(isotopes_tmp)
            mol_parity = parity4([i if i in isotopes_tmp else only_in_src.pop() for i in isotopes_mol])
            if PLEVEL >= 10: print(str(tmp_parity))
            if PLEVEL >= 10: print(str(mol_parity))
            parity_matches = tmp_parity == mol_parity
            tag_matches = a_tmp.GetChiralTag() == a_mol.GetChiralTag()
            chirality_matches = parity_matches == tag_matches
            if PLEVEL >= 2: print('Isotope {} chiral match? {}'.format(a_tmp.GetIsotope(), chirality_matches))
            return 1 if chirality_matches else -1
        else:
            if PLEVEL >= 2: print('Isotope {} chiral match? Based on isotope lists, ambiguous -> True'.format(a_tmp.GetIsotope()))
            return 2 # ambiguous case, just return for now
            # TODO: fix this?

    except IndexError as e:
        print(a_tmp.GetPropsAsDict())
        print(a_mol.GetPropsAsDict())
        print(a_tmp.GetChiralTag())
        print(a_mol.GetChiralTag())
        print(str(e))
        print(str(isotopes_tmp))
        print(str(isotopes_mol))
        raise KeyError('Pop from empty set - this should not happen!')