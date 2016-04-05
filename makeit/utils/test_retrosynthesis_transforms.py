from __future__ import print_function
import argparse
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
from makeit.utils.draw import *
import os

def load_transforms(tform_file):
	'''This function loads retrosynthetic transforms from a text file, 
	where each line consists of '# \t SMARTS_string \t score'. For now,
	the score is the frequency of that transform's occurrence in the 
	original database.'''
	
	rxns = []
	scores = []
	tforms = []
	with open(tform_file, 'r') as in_fid:
		for i, line in enumerate(in_fid):
			line_split = line.split('\t')
			smarts = line_split[1]

			for force_intramolecular in [False, True]:
				# Force grouped products?
				if force_intramolecular: smarts = smarts.replace('>>', '>>(') + ')'

				# Get transform as RDKit reaction
				rxn = AllChem.ReactionFromSmarts(smarts)

				# Check
				if rxn.Validate() != (0, 0):
					if v: print('Couldn\'t load transform {}'.format(smarts))
					continue
				
				# Save
				if rxn:
					rxns.append(rxn)
					tforms.append(smarts)

					# Score?
					if len(line_split) >= 3:
						scores.append(float(line_split[2]))
	
	# Default scores to all ones
	if not scores: scores = [1] * len(rxns)
	print('Successfully loaded {} out of {} transforms'.format(len(rxns), i ))
	return rxns, scores, tforms

def save_retrosynthesis(precursors, folder, target_smiles, n = 50, label = '', ext = 'jpg'):
	'''This function generates an html page for the single step'''

	# Prepare folders
	if not os.path.exists(folder):
		os.makedirs(folder)
	media_folder = os.path.join(folder, 'media')
	if not os.path.exists(media_folder):
		os.makedirs(media_folder)

	# Write text file as 'num \t precursors \t score'
	# (sorted by score)
	out_file = os.path.join(folder, 'full_retro_list.txt')
	with open(out_file, 'w') as out_fid:
		if label:
			out_fid.write('For the synthesis of {} ({}):\n'.format(target_smiles, label))
		else:
			out_fid.write('For the synthesis of {}:\n'.format(target_smiles))
		for i, precursor in enumerate(sorted(precursors, key = precursors.get, reverse = True)):
			out_fid.write('{}\t{}\t{}\n'.format(i, precursor, precursors[precursor]))

	# Write webpage
	with open(os.path.join(folder, 'index.html'), 'w') as out_fid:
		# Begin webpage
		out_fid.write('<html>\n')
		out_fid.write('\t<body>\n')
		if label:
			out_fid.write('\t\t<h1>List of one-step precursors to {} ({})</h1>\n'.format(target_smiles, label))
		else:
			out_fid.write('\t\t<h1>List of one-step precursors to {}</h1>\n'.format(target_smiles))
		img = MolsSmilesToImage(target_smiles)
		if ext != 'png':
			img = StripAlphaFromImage(img)
		img.save(os.path.join(media_folder, 'target.{}'.format(ext)))
		out_fid.write('\t\t<img src="media/target.{}"><br>\n'.format(ext))
		out_fid.write('\t\t<table border="1" align="center">\n')

		# Headers
		out_fid.write('\t\t\t<tr>\n')
		out_fid.write('\t\t\t\t<th align="center">Rank</th>\n')
		out_fid.write('\t\t\t\t<th align="center">\n')
		out_fid.write('\t\t\t\t\t<font face="Courier New, Courier monospace">Precursor(s)</font><br>\n')
		out_fid.write('\t\t\t\t</th>\n')
		out_fid.write('\t\t\t\t<th align="center">Score</th>\n')
		out_fid.write('\t\t\t</tr>\n')

		# Add precursors
		for i, precursor in enumerate(sorted(precursors, key = precursors.get, reverse = True)):
			if n and i == n: break
			img = MolsSmilesToImage(precursor.replace(' + ', '.'))
			if ext != 'png':
				img = StripAlphaFromImage(img)
			img.save(os.path.join(media_folder, 'precursor{}_score{}.{}'.format(i, precursors[precursor], ext)))

			# Add to table
			out_fid.write('\t\t\t<tr align="center">\n')
			out_fid.write('\t\t\t\t<td align="center">{}</td>\n'.format(i))
			out_fid.write('\t\t\t\t<td align="center">\n')
			out_fid.write('\t\t\t\t\t<font face="Courier New, Courier monospace">{}</font><br>\n'.format(precursor))
			out_fid.write('\t\t\t\t\t<img src="media/precursor{}_score{}.{}">\n'.format(i, precursors[precursor], ext))
			out_fid.write('\t\t\t\t</td>\n')
			out_fid.write('\t\t\t\t<td align="center">{}</td>\n'.format(precursors[precursor]))
			out_fid.write('\t\t\t</tr>\n')

		# Finish webpage
		out_fid.write('\t\t</table>\n')
		out_fid.write('\t</body>\n')
		out_fid.write('</html>\n')


def main(tform_file, target_smiles, folder, label = '', ext = 'jpg'):
	'''Given a SMILES string representation of a molecule and a list of 
	one-step retrosynthesis transforms, this function generates a list 
	of unique reactant sets that could feasibly make the target molecule 
	in a single step. Optionally, scores for each transform are provided;
	if multiple transforms suggest the same starting materials, then the
	scores are simply added.'''

	# Load data
	rxns, scores, tforms = load_transforms(tform_file)

	# Get molecule and re-write SMILES
	mol = Chem.MolFromSmiles(target_smiles)
	if not mol:
		print('error: could not parse target_smiles: {}'.format(target_smiles))
		quit(1)
	Chem.SanitizeMol(mol)
	target_smiles = Chem.MolToSmiles(mol, isomericSmiles = True)

	precursors = {}
	for i, rxn in enumerate(rxns):
		outcomes = rxn.RunReactants([mol])
		if not outcomes: continue
		for j, outcome in enumerate(outcomes):
			try:
				[x.UpdatePropertyCache() for x in outcome]
				[Chem.SanitizeMol(x) for x in outcome]
			except Exception as e:
				if v:
					print('warning: could not sanitize products from transform {}'.format(tforms[i]))
					[print('    {}'.format(Chem.MolToSmiles(x, isomericSmiles = True))) for x in outcome]
				continue
			precursor = ' + '.join(sorted([Chem.MolToSmiles(x, isomericSmiles = True) for x in outcome]))
			if v: 
				print('Candidate precursors ({}) found with score {}'.format(precursor, scores[i]))
				print('using transform: {}'.format(tforms[i]))
			if precursor not in precursors:
				precursors[precursor] = scores[i]
			else:
				precursors[precursor] += scores[i]

	save_retrosynthesis(precursors, folder, target_smiles, label = label, ext = ext)


if __name__ == '__main__':

	parser = argparse.ArgumentParser()

	parser.add_argument('tform_file', type = str, 
		 				help = 'File where each line is a mapped SMARTS retrosynthesis')
	parser.add_argument('-s', '--smiles', type = str, 
						help = 'SMILES string of target molecule')
	parser.add_argument('-n', '--name', type = str, default = '',
						help = 'Label for target molecule; defaults to None')
	parser.add_argument('-o', '--out', type = str,
						help = 'Folder to output candidate starting materials to')
	parser.add_argument('-e', '--ext', type = str, default = 'jpg',
		  				help = 'Extension of images; defaults to jpg')
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing; defaults to False')
	args = parser.parse_args()

	v = args.v
	main(args.tform_file, args.smiles, args.out, label = args.name, ext = args.ext)