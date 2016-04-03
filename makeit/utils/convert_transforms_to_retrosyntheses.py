from __future__ import print_function
import argparse
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
import sys  # for commanad line
import os   # for file paths
import re 

def convert_to_retro(transform):
	'''This function takes a forward synthesis and converts it to a
	retrosynthesis. Only transforms with a single product are kept, since 
	retrosyntheses should have a single reactant (and split it up accordingly).'''
	
	reactants = transform.split('>>')[0]
	products  = transform.split('>>')[1]
	if ').(' in products: return None # Multiple product mols
	return '>>'.join([products, reactants])

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

			# Convert to retrosynthesis
			retro_transform = convert_to_retro(transform)
			
			# Add to dictionary
			if retro_transform:
				if v: print('post: {}'.format(retro_transform))
				transform_dict[retro_transform] = int(line.split('\t')[2].strip())
			else:
				if v: print('unsuitable as retrosynthesis step!')

			if v:
				# Pause
				raw_input('Enter anything to continue...')

	# Write as 'counter' \t 'canonical transform' \t 'number of records'
	# (sorted by frequency of occurrence)
	with open(tform_out, 'w') as out_fid:
		for i, tform in enumerate(sorted(transform_dict, key = transform_dict.get, reverse = True)):
			out_fid.write('{}\t{}\t{}\n'.format(i, tform, transform_dict[tform]))

	# Report
	print('Saved {} unique retrosynthesis transforms'.format(len(transform_dict)))

	return 

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('data_file', type = str, 
		 				help = 'File where each line is a mapped SMARTS reaction template')
	parser.add_argument('out_file', type = str,
						help = 'File to output retrosyntheses to')
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	args = parser.parse_args()

	v = args.v
	main(args.data_file, args.out_file)