import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.Draw as Draw
from PIL import Image, ImageOps
from collections import defaultdict
from rdkit.Chem.Draw.cairoCanvas import Canvas 
import os
import numpy as np
import re

'''
Many of these functions are taken from RDKit.
'''

def mols_from_smiles_list(all_smiles):
	'''Given a list of smiles strings, this function creates rdkit
	molecules'''
	mols = []
	for smiles in all_smiles:
		if not smiles: continue
		mols.append(Chem.MolFromSmiles(smiles))
	return mols

def defaultDrawOptions():
	'''This function returns an RDKit drawing options object with 
	default drawing options.'''

	opts = Draw.DrawingOptions()
	#opts.elemDict = defaultdict(lambda: (0,0,0)) # all atoms are black
	opts.noCarbonSymbols = True
	opts.selectColor = (1, 0, 0)
	opts.wedgeBonds = True
	return opts

def StripAlphaFromImage(img):
	'''This function takes an RGBA PIL image and returns an RGB image'''

	if len(img.split()) == 3: return img
	return Image.merge('RGB', img.split()[:3])

def MolToImage(mol, max_size = (1000, 1000), kekulize = True, options = None,
	canvas = None, **kwargs):
	'''Wrapper for RDKit's MolToImage. If mol == None, an arrow is drawn'''

	if not options:   options = defaultDrawOptions()
	if mol:
		return Draw.MolToImage(mol, size = max_size, kekulize = kekulize, options = options, 
			canvas = canvas, **kwargs)
	else: # mol == None means this is an arrow
		subImgSize = (250, 250)
		img,canvas = Draw._createCanvas(subImgSize) 
		p0 = (10,subImgSize[1]//2) 
		p1 = (subImgSize[0]-10,subImgSize[1]//2) 
		p3 = (subImgSize[0]-20,subImgSize[1]//2-10) 
		p4 = (subImgSize[0]-20,subImgSize[1]//2+10) 
		canvas.addCanvasLine(p0,p1,lineWidth=2,color=(0,0,0)) 
		canvas.addCanvasLine(p3,p1,lineWidth=2,color=(0,0,0)) 
		canvas.addCanvasLine(p4,p1,lineWidth=2,color=(0,0,0)) 
		if hasattr(canvas,'flush'): 
			canvas.flush() 
		else: 
			canvas.save()
		img.save('test_arrow.png')
		return img

def TrimImgByWhite(img, padding = 0):
	'''This function takes a PIL image, img, and crops it to the minimum rectangle 
	based on its whiteness/transparency. 5 pixel padding used automatically.'''

	# Convert to array
	as_array = np.array(img) # N x N x (r,g,b,a)

	# Set previously-transparent pixels to white
	as_array[as_array[:, :, 3] == 0] = [255, 255, 255, 255]

	# Content defined as non-white and non-transparent pixel
	has_content = np.sum(as_array, axis = 2, dtype = np.uint32) != 255 * 4
	xs, ys = np.nonzero(has_content)

	# Crop down
	x_range = max([min(xs) - 5, 0]), min([max(xs) + 5, as_array.shape[0]])
	y_range = max([min(ys) - 5, 0]), min([max(ys) + 5, as_array.shape[1]])
	as_array_cropped = as_array[x_range[0]:x_range[1], y_range[0]:y_range[1], 0:3]

	img = Image.fromarray(as_array_cropped, mode = 'RGB')

	return ImageOps.expand(img, border = padding, fill = (255, 255, 255, 0))

def StitchPILsHorizontally(imgs):
	'''This function takes a list of PIL images and concatenates
	them onto a new image horizontally, with each one
	vertically centered.'''

	# Create blank image (def: transparent white)
	heights = [img.size[1] for img in imgs]
	height = max(heights)
	widths = [img.size[0] for img in imgs]
	width = sum(widths)
	res = Image.new('RGBA', (width, height), (255, 255, 255, 0))

	# Add in sub-images
	for i, img in enumerate(imgs):
		offset_x = sum(widths[:i]) # left to right
		offset_y = (height - heights[i]) / 2
		res.paste(img, (offset_x, offset_y))

	return res

def CheckAtomForGeneralization(atom):
	'''Given an RDKit atom, this function determines if that atom's SMART 
	representation was likely a result of generalization. This assumes that
	the transform string was generated using explicit Hs with aliphatic 
	carbons as C, aromatic carbons as c, and non-carbons as #N where N is the 
	atomic number of the generalized species.'''

	smarts = atom.GetSmarts()
	# Check if this was a result of generalization
	# non-carbon atom, generlized
	if '#' in smarts: 
		atomSymbol = atom.GetSymbol()
		atom.SetAtomicNum(0)
		atom.SetProp('dummyLabel', '[{}]'.format(atomSymbol))
		atom.UpdatePropertyCache()
	# aliphatic carbon, generalized (all non-generalized use explicit Hs)
	elif '[C' in smarts and 'H' not in smarts:
		atom.SetAtomicNum(0)
		atom.SetProp('dummyLabel', 'C[al]')
		atom.UpdatePropertyCache()
	elif '[c' in smarts and 'H' not in smarts:
		atom.SetAtomicNum(0)
		atom.SetProp('dummyLabel', 'C[ar]')
		atom.UpdatePropertyCache()


def ReactionToImage(rxn, dummyAtoms = False, options = None, **kwargs):
	'''Modification of RDKit's ReactionToImage to allow for each molecule 
	to have a different drawn size. rxn is an RDKit reaction object

	warning: this function adds hydrogens as it sees fit'''

	# Extract mols from reaction
	mols = []
	for i in range(rxn.GetNumReactantTemplates()):
		mol = rxn.GetReactantTemplate(i)
		mol.UpdatePropertyCache(False)
		mols.append(mol)
		if dummyAtoms: [CheckAtomForGeneralization(atom) for atom in mol.GetAtoms()]


	mols.append(None) # placeholder for arrow
	for j in range(rxn.GetNumProductTemplates()):
		mol = rxn.GetProductTemplate(j)
		mol.UpdatePropertyCache(False)
		mols.append(mol)
		if dummyAtoms: [CheckAtomForGeneralization(atom) for atom in mol.GetAtoms()]

	# Generate images for all molecules/arrow
	imgs = [TrimImgByWhite(MolToImage(mol, kekulize = False, options = options), padding = 15) for mol in mols]

	# Combine
	return StitchPILsHorizontally(imgs)

def ReactionStringToImage(rxn_string, strip = True, update = True, options = None, **kwargs):
	'''This function takes a SMILES rxn_string as input, not an 
	RDKit reaction object, and draws it.'''

	reactants, agents, products = [mols_from_smiles_list(x) for x in 
										[mols.split('.') for mols in rxn_string.split('>')]]
	if None in reactants + agents + products:
		raise ValueError('Could not parse entirety of reaction: {}'.format(rxn_string))

	# Stich together mols (ignore agents)
	mols = reactants + [None] + products
	if update: [mol.UpdatePropertyCache(False) for mol in mols if mol]
	if strip: mols = [Chem.MolFromInchi(Chem.MolToInchi(mol)) if mol else None for mol in mols]

	# Generate images
	imgs = [TrimImgByWhite(MolToImage(mol, kekulize = False, options = options), padding = 15) for mol in mols]

	# Combine
	return StitchPILsHorizontally(imgs)

def TransformStringToImage(transform, **kwargs):
	'''Wrapper function meant to take a SMARTS transform and return a PIL image
	of that transform.

	TODO: Need to improve generalization visually! Right now it still shows'''

	options = defaultDrawOptions()
	options.dotsPerAngstrom = 50
	rxn = AllChem.ReactionFromSmarts(transform)
	return ReactionToImage(rxn, dummyAtoms = True, options = options, **kwargs)

def MolsSmilesToImage(smiles, options = None, **kwargs):
	'''This function takes a SMILES string of one or more molecules
	and generates a combined image for that molecule set.'''

	# Generate mols
	mols = mols_from_smiles_list(smiles.split('.'))
	# Generate images
	imgs = [TrimImgByWhite(MolToImage(mol, kekulize = False, options = options), padding = 15) for mol in mols]
	# Combine
	return StitchPILsHorizontally(imgs)


if __name__ == '__main__':

	# Simple test cases
	rxn_string = '[Na+].[CH3:2][C:3](=[O:5])[O-].[CH3:6][c:7]1[cH:12][cH:11][cH:10][cH:9][cH:8]1>>CN3[C@H]1CC[C@@H]3C[C@@H](C1)OC(=O)C(CO)c2ccccc2.[c:7]1([CH3:6])[c:12]([C:3]([c:2]2[cH:11][cH:12][cH:7][cH:8][c:9]2[CH3:10])=[O:5])[cH:11][cH:10][cH:9][cH:8]1'
	rxn = AllChem.ReactionFromSmarts(rxn_string)
	rxn_image = ReactionToImage(rxn)
	rxn_image.save('test_rxn.png')
	rxn_image_string = ReactionStringToImage(rxn_string, strip = True)
	rxn_image_string.save('test_rxn_string.png')
	tform = '([O;H0:3]=[C;H0:4](-[C:5])-[NH:2]-[C:1])>>([C:1]-[NH2:2]).([OH:3]-[C;H0:4](=O)-[C:5])'
	img = TransformStringToImage(tform)
	img.save('transform.png')