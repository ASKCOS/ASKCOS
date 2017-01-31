'''
Module meant to work with DFTB+ for calculating atom-level descriptors
'''

import os
import subprocess
import rdkit.Chem as Chem 
import rdkit.Chem.AllChem as AllChem

dftb_root = '/home/ccoley/dftb+/dftbplus-1.3.0.x86_64-linux/'
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
    InitMixingParameter = 0.1
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
	AllChem.EmbedMolecule(new_mol)
	conf = new_mol.GetConformer()

	symbols = list(set([a.GetSymbol() for a in new_mol.GetAtoms()]))
	num_atoms = len(new_mol.GetAtoms())

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
			raise ValueError('Cannot run DFTB+ on this molecule because symbol {} does not have a defined MaxAngularMomentum!'.format(symbol))

	with open(dftb_in, 'w') as fid:
		fid.write(template_input_file % (charge, maxangularmomentum))

	return new_mol

def read_results(mol, new_mol):
	'''given the results file of a DFTB+ calculation, add charges to mol'''

	# Results "atom num" is 1 + AtomIdx for the RDKit object

	# Which atoms are hydrogens? Keep track of their "parent" atom
	atom_idx_map = dict()
	hydrogen_idx = set()
	for a in new_mol.GetAtoms():
		if a.GetAtomicNum() == 1:
			atom_idx_map[a.GetIdx()] = a.GetNeighbors()[0].GetIdx()
			hydrogen_idx.add(a.GetIdx())

	skipped_lines = 0
	with open(os.path.join(dftb_root, 'detailed.out'), 'r') as fid:
		for line in fid:
			if skipped_lines == 0:
				if 'Net atomic charges' not in line: continue

			# Skip "Net atomic charge" line and header line
			if skipped_lines < 2:
				skipped_lines += 1
				continue 
			line_split = [x for x in line.strip().split(' ') if x]

			# Done with charge table?
			if not line_split: break

			# Parse
			idx = int(line_split[0]) - 1
			charge = float(line_split[1])

			# Set values in explicit-H molecule
			new_mol.GetAtomWithIdx(idx).SetDoubleProp('DFTB_pcharge', charge)

			# Set values in implicit-H molecule
			if idx in hydrogen_idx:
				a = mol.GetAtomWithIdx(atom_idx_map[idx])
				prev_H_pcharge = a.GetDoubleProp('DFTB_H_pcharge') if a.HasProp('DFTB_H_pcharge') else 0.0
				a.SetDoubleProp('DFTB_H_pcharge', prev_H_pcharge + charge)
			else:
				mol.GetAtomWithIdx(idx).SetDoubleProp('DFTB_pcharge', charge)

	# Set H_pcharge for remaining atoms
	for a in mol.GetAtoms():
		if not a.HasProp('DFTB_H_pcharge'):
			a.SetDoubleProp('DFTB_H_pcharge', 0.0)

if __name__ == '__main__':

	while True:
		smiles = raw_input('Enter smiles: ')
		mol = Chem.MolFromSmiles(smiles)
		if not mol: continue 


		charges = [0.0, -1.0, +1.0]

		for charge in charges:

			new_mol = make_input_files(mol, charge = charge)
			os.chdir(os.path.dirname(dftb_root))
			subprocess.call(dftb, shell = True)
			read_results(mol, new_mol)

			print('Gastieger:')
			import rdkit.Chem.rdPartialCharges as rdPartialCharges
			rdPartialCharges.ComputeGasteigerCharges(mol)
			print([(a.GetIdx(), a.GetSymbol(), a.GetProp('_GasteigerCharge'), a.GetProp('_GasteigerHCharge')) for a in mol.GetAtoms()])

			print('DFTB+:')
			print([(a.GetIdx(), a.GetSymbol(), a.GetProp('DFTB_pcharge'), a.GetProp('DFTB_H_pcharge')) for a in mol.GetAtoms()])

			if charge:
				fname = os.path.join(dftb_root, '{} charge{}.dat'.format(smiles, charge))
			else:
				fname = os.path.join(dftb_root, '{}.dat'.format(smiles))
			with open(fname, 'w') as fid:
				fid.write('smiles: {}, DFTB charge: {}\n'.format(smiles, charge))
				fid.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
					'Idx', 'Symbol', 'GasteigerCharge', 'GasteigerHCharge', 'DFTB_pcharge', 'DFTB_H_pcharge'
				))
				for a in mol.GetAtoms():
					fid.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
						a.GetIdx(), a.GetSymbol(), a.GetProp('_GasteigerCharge'), a.GetProp('_GasteigerHCharge'), a.GetProp('DFTB_pcharge'), a.GetProp('DFTB_H_pcharge')
					))