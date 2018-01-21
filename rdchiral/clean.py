from __future__ import print_function
import rdkit.Chem as Chem 
import re
from itertools import chain

from rdchiral.utils import vprint, PLEVEL


def canonicalize_outcome_smiles(smiles, ensure=True):
    # Uniquify via SMILES string - a little sloppy
    # Need a full SMILES->MOL->SMILES cycle to get a true canonical string
    # also, split by '.' and sort when outcome contains multiple molecules
    if ensure or not keep_isotopes: 
        outcome = Chem.MolFromSmiles(smiles)
        if outcome is None:
            if PLEVEL >= 1: print('~~ could not parse self?')
            if PLEVEL >= 1: print('Attempted SMILES: {}', smiles)
            return None

        smiles = Chem.MolToSmiles(outcome, True)

    return  '.'.join(sorted(smiles.split('.')))

def combine_enantiomers_into_racemic(final_outcomes):
    '''
    If two products are identical except for an inverted CW/CCW or an
    opposite cis/trans, then just strip that from the product. Return
    the achiral one instead.
    
    This is not very sophisticated, since the chirality could affect the bond
    order and thus the canonical SMILES. But, whatever. It also does not look
    to invert multiple stereocenters at once
    '''

    for smiles in list(final_outcomes)[:]:

        # Look for @@ tetrahedral center
        for match in re.finditer(r'@@', smiles):
            smiles_inv = '%s@%s' % (smiles[:match.start()], smiles[match.end():])
            if smiles_inv in final_outcomes:
                if smiles in final_outcomes:
                    final_outcomes.remove(smiles)
                final_outcomes.remove(smiles_inv)
                # Re-parse smiles so that hydrogens can become implicit
                smiles = smiles[:match.start()] + smiles[match.end():]
                outcome = Chem.MolFromSmiles(smiles)
                if outcome is None:
                    raise ValueError('Horrible mistake when fixing duplicate!')
                smiles = '.'.join(sorted(Chem.MolToSmiles(outcome, True).split('.')))
                final_outcomes.add(smiles)

        # Look for // or \\ trans bond
        # where [^=\.] is any non-double bond or period or slash
        for match in chain(re.finditer(r'(\/)([^=\.\\\/]+=[^=\.\\\/]+)(\/)', smiles), 
                re.finditer(r'(\\)([^=\.\\\/]+=[^=\.\\\/]+)(\\)', smiles)):
            # See if cis version is present in list of outcomes
            opposite = {'\\': '/', '/': '\\'}
            smiles_cis1 = '%s%s%s%s%s' % (smiles[:match.start()], 
                match.group(1), match.group(2), opposite[match.group(3)],
                smiles[match.end():]) 
            smiles_cis2 = '%s%s%s%s%s' % (smiles[:match.start()], 
                opposite[match.group(1)], match.group(2), match.group(3),
                smiles[match.end():])
            # Also look for equivalent trans
            smiles_trans2 = '%s%s%s%s%s' % (smiles[:match.start()], 
                opposite[match.group(1)], match.group(2), 
                opposite[match.group(3)], smiles[match.end():])
            # Kind of weird remove conditionals...
            remove = False
            if smiles_cis1 in final_outcomes:
                final_outcomes.remove(smiles_cis1)
                remove = True 
            if smiles_cis2 in final_outcomes:
                final_outcomes.remove(smiles_cis2)
                remove = True 
            if smiles_trans2 in final_outcomes and smiles in final_outcomes:
                final_outcomes.remove(smiles_trans2)
            if remove:
                final_outcomes.remove(smiles)
                smiles = smiles[:match.start()] + match.group(2) + smiles[match.end():]
                outcome = Chem.MolFromSmiles(smiles)
                if outcome is None:
                    raise ValueError('Horrible mistake when fixing duplicate!')
                smiles = '.'.join(sorted(Chem.MolToSmiles(outcome, True).split('.')))
                final_outcomes.add(smiles)
    return final_outcomes
