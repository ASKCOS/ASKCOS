'''
Module meant to work with DFTB+ for calculating atom-level descriptors
'''

import os
import subprocess
import rdkit.Chem as Chem 
import rdkit.Chem.AllChem as AllChem

dftb_root = os.path.join(os.path.dirname(__file__), 'dftbplus-1.3.0.x86_64-linux')
dftb = os.path.join(dftb_root, 'dftb+')
dftb_in = os.path.join(dftb_root, 'dftb_in.hsd')



template_input_file = '''
Geometry = GenFormat {
  <<< "geom.gen"
}

Driver = gDIIS {
  MovedAtoms = 1:-1
  MaxForceComponent = 1E-4
  MaxSteps = 100
  OutputPrefix = "geom.out"
}

Hamiltonian = DFTB {
  SCC = Yes
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/home/ccoley/dftb+/dftbplus-1.3.0.x86_64-linux/3ob-3-1/"
    Separator = "-"
    Suffix = ".skf"
    LowerCaseTypeName = No
  }
  Charge = %f
  MaxAngularMomentum {
%s
  }
  MaxSCCIterations = 500
  Mixer = DIIS{
    InitMixingParameter = 0.2
    Generations = 10
  }
  SCCTolerance = 1e-7
}

ParserOptions {
  ParserVersion = 4
}
'''

# Max angular momenta
L_max = {
    'Br' : 'd',    
    'Mg' : 'p',
    'C'  : 'p',    
    'N'  : 'p',
    'Ca' : 'p',    
    'Na' : 'p',
    'Cl' : 'd',    
    'O'  : 'p',
    'F'  : 'p',    
    'P'  : 'p',
    'H'  : 's',    
    'S'  : 'd',
    'I'  : 'd',    
    'Zn' : 'd',
    'K'  : 'p',
}

def make_input_files(mol, charge = 0.0):

    new_mol = AllChem.AddHs(mol)
    new_mol = replace_invalid_atoms(new_mol)
    
    symbols = list(set([a.GetSymbol() for a in new_mol.GetAtoms()]))
    num_atoms = len(new_mol.GetAtoms())

    # Do we need to generate a geometry? Only if this is the first time
    if charge == 0.0:
        AllChem.EmbedMolecule(new_mol)
        conf = new_mol.GetConformer()
        with open(os.path.join(dftb_root, 'geom.gen'), 'w') as fid:
            fid.write('{} C\n'.format(num_atoms))
            fid.write('  {}\n'.format(' '.join(symbols)))
            for (i, a) in enumerate(new_mol.GetAtoms()):
                fid.write('%i\t%i\t%f\t%f\t%f\n' % (i + 1, symbols.index(a.GetSymbol()) + 1, 
                    conf.GetAtomPosition(a.GetIdx()).x, conf.GetAtomPosition(a.GetIdx()).y, conf.GetAtomPosition(a.GetIdx()).z))

    maxangularmomentum = ''
    for symbol in symbols:
        try:
            maxangularmomentum += '    {} = "{}"\n'.format(symbol, L_max[symbol])
        except KeyError:
            raise ValueError('Cannot run DFTB+ on this molecule because symbol {} does not have a defined MaxAngularMomentum! Could it be substituted for one of the elements for which we have parameters? Check {}'.format(symbol, __file__))

    with open(dftb_in, 'w') as fid:
        if charge == 0.0:
            fid.write(template_input_file % (charge, maxangularmomentum))
        else:
            # Use optimized neutral geometry
            fid.write(template_input_file.replace('geom.gen', 'geom.out.gen').replace('MaxSteps = 200', 'MaxSteps = 0') % (charge, maxangularmomentum))

    return new_mol

def read_results(mol, new_mol, energylabel = 'energy', chargelabel = 'dftbCharge', Hchargelabel = 'dftbHCharge'):
    '''given the results file of a DFTB+ calculation, add charges to mol'''

    # Results "atom num" is 1 + AtomIdx for the RDKit object

    # Which atoms are hydrogens? Keep track of their "parent" atom
    atom_idx_map = dict()
    hydrogen_idx = set()
    for a in new_mol.GetAtoms():
        if a.GetAtomicNum() == 1:
            atom_idx_map[a.GetIdx()] = a.GetNeighbors()[0].GetIdx()
            hydrogen_idx.add(a.GetIdx())

    skip_for_energy = 0
    skip_for_charges = 0
    with open(os.path.join(dftb_root, 'detailed.out'), 'r') as fid:
        for line in fid:
            if '************' in line:
                skip_for_energy += 1
                continue
            if skip_for_energy == 1:
                skip_for_energy += 1
                continue
            if skip_for_energy == 2:
                line_split = [x for x in line.strip().split(' ') if x]
                new_mol.SetDoubleProp(energylabel, float(line_split[1]))
                mol.SetDoubleProp(energylabel, float(line_split[1]))
                skip_for_energy += 1
                continue
            if skip_for_charges == 0:
                if 'Net atomic charges' not in line: continue

            # Skip "Net atomic charge" line and header line
            if skip_for_charges < 2:
                skip_for_charges += 1
                continue 
            line_split = [x for x in line.strip().split(' ') if x]

            # Done with charge table?
            if not line_split: break

            # Parse
            idx = int(line_split[0]) - 1
            charge = float(line_split[1])

            # Set values in explicit-H molecule
            new_mol.GetAtomWithIdx(idx).SetDoubleProp(chargelabel, charge)

            # Set values in implicit-H molecule
            if idx in hydrogen_idx:
                a = mol.GetAtomWithIdx(atom_idx_map[idx])
                prev_H_pcharge = a.GetDoubleProp(Hchargelabel) if a.HasProp(Hchargelabel) else 0.0
                a.SetDoubleProp(Hchargelabel, prev_H_pcharge + charge)
            else:
                mol.GetAtomWithIdx(idx).SetDoubleProp(chargelabel, charge)

    # Set H_pcharge for remaining atoms
    for a in mol.GetAtoms():
        if not a.HasProp(Hchargelabel):
            a.SetDoubleProp(Hchargelabel, 0.0)


def replace_invalid_atoms(new_mol):
    '''
    Because DFTB+ only has parameters for some elements, we need to replace others
    to trick the solver into working'''

    for a in new_mol.GetAtoms():
        if a.GetSymbol() in L_max: continue 
        if a.GetAtomicNum() in range(21, 31) or a.GetAtomicNum() in range(39, 49): # transition metals
            a.SetAtomicNum(30) # zinc
        elif a.GetAtomicNum() == 3: # lithium
            a.SetAtomicNum(11) # sodium
        elif a.GetSymbol() == 'Si':
            a.SetAtomicNum(6) # carbon
        else:
            raise ValueError('Atom {} not in list of possible DFTB elements, nor in our list of custom replacements'.format(a.GetSymbol()))
    return new_mol


def atom_dftb(mol, v = False, to_file = False):
    '''Given a molecule, run a DFTB calculation and get descriptors'''

    attributes = [[] for i in range(mol.GetNumAtoms())]

    settings = [
        (0.0, ''),
        (-1.0, 'neg_'),
        (+1.0, 'pos_'),
    ]

    for (charge, label) in settings:

        energylabel  = label + 'energy'
        chargelabel  = label + 'dftbCharge'
        Hchargelabel = label + 'dftbHCharge'

        new_mol = make_input_files(mol, charge = charge)
        os.chdir(os.path.dirname(dftb_root))
        with open(os.path.join(os.path.dirname(dftb_root), 'log.txt'), 'w') as fid:
            subprocess.call(dftb, shell = True, stdout = fid)
        read_results(mol, new_mol, energylabel = energylabel, chargelabel = chargelabel, Hchargelabel = Hchargelabel)

        #print('Gastieger:')
        #import rdkit.Chem.rdPartialCharges as rdPartialCharges
        #rdPartialCharges.ComputeGasteigerCharges(mol)
        #print([(a.GetIdx(), a.GetSymbol(), a.GetProp('_GasteigerCharge'), a.GetProp('_GasteigerHCharge')) for a in mol.GetAtoms()])

        if v:
            print('DFTB+:')
            print([(a.GetIdx(), a.GetSymbol(), a.GetProp(chargelabel), a.GetProp(Hchargelabel)) for a in mol.GetAtoms()])

            print('Energy:')
            print(mol.GetProp(energylabel))


    # Now perform molecule-level calculations
    IP = mol.GetDoubleProp('pos_energy') - mol.GetDoubleProp('energy')
    EA = mol.GetDoubleProp('neg_energy') - mol.GetDoubleProp('energy')
    mu = - (IP + EA) / 2.0 # negative chem potential
    eta = IP - EA # hardness
    S = 1.0 / eta # softness
    omega = mu ** 2.0 / (2.0 * eta) # philicity


    # Apply atom-level calculations
    coefficients = [
        (S, 'S_'),
        (omega, 'W_'),
    ]

    for i, a in enumerate(mol.GetAtoms()):
        a.SetDoubleProp('Fukui_nuc', - a.GetDoubleProp('neg_dftbCharge') + a.GetDoubleProp('dftbCharge'))
        a.SetDoubleProp('Fukui_elec', - a.GetDoubleProp('dftbCharge') + a.GetDoubleProp('pos_dftbCharge'))
        a.SetDoubleProp('Fukui_rad', - (a.GetDoubleProp('neg_dftbCharge') - a.GetDoubleProp('pos_dftbCharge')) / 2.0)
        a.SetDoubleProp('Fukui', a.GetDoubleProp('Fukui_nuc') - a.GetDoubleProp('Fukui_elec'))

        for (coeff, label) in coefficients:
            a.SetDoubleProp(label + 'Fukui_nuc', a.GetDoubleProp('Fukui_nuc') * coeff)
            a.SetDoubleProp(label + 'Fukui_elec', a.GetDoubleProp('Fukui_elec') * coeff)
            a.SetDoubleProp(label + 'Fukui_rad', a.GetDoubleProp('Fukui_rad') * coeff)
            a.SetDoubleProp(label + 'Fukui', a.GetDoubleProp('Fukui') * coeff)

    
    cols = [
        #'_GasteigerCharge',
        #'_GasteigerHCharge',
        'dftbCharge',
        'dftbHCharge',
        'pos_dftbCharge',
        'pos_dftbHCharge',
        'neg_dftbCharge',
        'neg_dftbHCharge',
        'Fukui_nuc',
        'Fukui_elec',
        'Fukui_rad',
        'Fukui',
        'S_Fukui_nuc',
        'S_Fukui_elec',
        'S_Fukui_rad',
        'S_Fukui',
        'W_Fukui_nuc',
        'W_Fukui_elec',
        'W_Fukui_rad',
        'W_Fukui',
    ]

    if to_file:
        with open(os.path.join(os.path.dirname(dftb_root), '{}.dat'.format(smiles)), 'w') as fid: 
            fid.write('SMILES\t{}\n'.format(smiles))
            fid.write('Energy\t{}\n'.format(mol.GetDoubleProp('energy')))
            fid.write('IP\t{}\n'.format(IP))
            fid.write('EA\t{}\n'.format(EA))
            fid.write('Mu (chem pot)\t{}\n'.format(mu))
            fid.write('Eta (hardness)\t{}\n'.format(eta))
            fid.write('S (softness)\t{}\n'.format(S))
            fid.write('Omega (philicity)\t{}\n'.format(omega))

            fid.write('Idx\tSymbol\t' + '\t'.join(cols) + '\n')
            for a in mol.GetAtoms():
                fid.write('{}\t{}\t'.format(a.GetIdx(), a.GetSymbol()) + '\t'.join([a.GetProp(label) for label in cols]) + '\n')

    for i, a in enumerate(mol.GetAtoms()):
        attributes[i].extend(
            [a.GetDoubleProp(label) for label in cols]
        )

    return attributes



if __name__ == '__main__':

    while True:
        smiles = raw_input('Enter smiles: ')
        mol = Chem.MolFromSmiles(smiles)
        if not mol: continue 


        attributes = atom_dftb(mol, v = False, to_file = False)
        print(attributes)