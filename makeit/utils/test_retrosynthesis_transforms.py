from __future__ import print_function
import argparse
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem

def load_transforms(tform_file):
	'''This function loads retrosynthetic transforms from a text file, 
	where each line consists of '# \t SMARTS_string \t score'. For now,
	the score is the frequency of that transform's occurrence in the 
	original database.'''
	
	rxns = []
	scores = []
	with open(tform_file, 'r') as in_fid:
		for i, line in enumerate(in_fid):
			line_split = line.split('\t')

			# Get transform as RDKit reaction
			rxn = AllChem.ReactionFromSmarts(line_split[1])

			# Check
			if rxn.Validate() != (0, 0):
				if v: print('Couldn\'t load transform {}'.format(line_split[1]))
				continue
			if not rxn: continue
			rxns.append(rxn)

			# Score?
			if len(line_split) >= 3:
				scores.append(float(line_split[2]))

	return rxns, scores



def main(tform_file, target_smiles, out_file):
	'''Given a SMILES string representation of a molecule and a list of 
	one-step retrosynthesis transforms, this function generates a list 
	of unique reactant sets that could feasibly make the target molecule 
	in a single step. Optionally, scores for each transform are provided;
	if multiple transforms suggest the same starting materials, then the
	scores are simply added.'''

	# Load data
	rxns, scores = load_transforms(tform_file)
	if not scores: scores = [1] * len(rxns)

	# Get molecule
	mol = Chem.MolFromSmiles(target_smiles)

	precursors = {}
	for i, rxn in enumerate(rxns):
		outcomes = rxn.RunReactants([mol])
		if not outcomes: continue
		for j, outcome in enumerate(outcomes):
			#[Chem.SanitizeMol(x) for x in outcome]
			precursor = ' + '.join(sorted([Chem.MolToSmiles(x, isomericSmiles = True) for x in outcome]))
			if v: print('Candidate precursors ({}) found with score {}'.format(precursor, scores[i]))
			if precursor not in precursors:
				precursors[precursor] = scores[i]
			else:
				precursors[precursor] += scores[i]

	# Write as 'num \t precursors \t score'
	# (sorted by score)
	with open(out_file, 'w') as out_fid:
		out_fid.write('For the synthesis of {}:\n'.format(target_smiles))
		for i, precursor in enumerate(sorted(precursors, key = precursors.get, reverse = True)):
			out_fid.write('{}\t{}\t{}\n'.format(i, precursor, precursors[precursor]))


if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('tform_file', type = str, 
		 				help = 'File where each line is a mapped SMARTS retrosynthesis')
	parser.add_argument('target_smiles', type = str, 
						help = 'SMILES string of target molecule')
	parser.add_argument('out_file', type = str,
						help = 'File to output candidate starting materials to')
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	args = parser.parse_args()

	v = args.v
	main(args.tform_file, args.target_smiles, args.out_file)