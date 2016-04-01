from __future__ import print_function
import argparse
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
import sys  # for commanad line
import os   # for file paths
import re 

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

def main(tform_in, tform_out):
	'''This function takes a list of atom-mapped SMARTS string reaction
	transforms in a tab delimited file (# \t SMARTS) and consolidates 
	them into a second transform file where duplicates are removed'''

	transform_dict = {}
	with open(tform_in, 'r') as in_fid:
		for i, line in enumerate(in_fid):

			# Print header
			if v: 
				print('\n-------------')
				print('Transform {}'.format(i))
				print('-------------\n')

			# Extract transform
			transform = line.split('\t')[1].strip()
			if v: print('pre:  {}'.format(transform))

			# Convert using auxiliary functions
			transform_canonical = canonicalize_transform(transform)
			if v: print('post: {}'.format(transform_canonical))

			# Add to dictionary
			if transform_canonical not in transform_dict:
				transform_dict[transform_canonical] = 1
				if v: print('(+) New transform!')
			else:
				transform_dict[transform_canonical] += 1
				if v: print('(-) Already seen this transform')

			if v:
				# Pause
				raw_input('Enter anything to continue...')

	# Write as 'counter' \t 'canonical transform' \t 'number of records'
	# (sorted by frequency of occurrence)
	with open(tform_out, 'w') as out_fid:
		for i, tform in enumerate(sorted(transform_dict, key = transform_dict.get, reverse = True)):
			out_fid.write('{}\t{}\t{}\n'.format(i, tform, transform_dict[tform]))

	# Report
	print('Saved {} unique transforms'.format(len(transform_dict)))

	return 

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('data_file', type = str, 
		 				help = 'File where each line is a mapped SMARTS reaction template')
	parser.add_argument('out_file', type = str,
						help = 'File to output consolidated SMARTS transforms to')
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	args = parser.parse_args()

	v = args.v
	main(args.data_file, args.out_file)