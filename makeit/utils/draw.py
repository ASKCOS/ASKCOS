import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.Draw as Draw
from PIL import Image, ImageOps
from collections import defaultdict
from rdkit.Chem.Draw.cairoCanvas import Canvas 
import os
import numpy as np

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
	opts.elemDict = defaultdict(lambda: (0,0,0)) # all atoms are black
	opts.noCarbonSymbols = True
	opts.selectColor = (1, 0, 0)
	opts.wedgeBonds = True
	return opts

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
	heights = [img.size[0] for img in imgs]
	height = max(heights)
	widths = [img.size[1] for img in imgs]
	width = sum(widths)
	res = Image.new('RGBA', (height, width), (255, 255, 255, 0))

	# Add in sub-images
	for i, img in enumerate(imgs):
		offset_x = sum(widths[:i]) # left to right
		offset_y = (height - heights[i]) / 2
		res.paste(img, (offset_y, offset_x))

	return res

def ReactionToImage(rxn, **kwargs):
	'''Modification of RDKit's ReactionToImage to allow for each molecule 
	to have a different drawn size. rxn is an RDKit reaction object'''

	# Extract mols from reaction
	mols = []
	for i in range(rxn.GetNumReactantTemplates()):
		mol = rxn.GetReactantTemplate(i)
		mol.UpdatePropertyCache(False)
		mols.append(mol)
	mols.append(None) # placeholder for arrow
	for j in range(rxn.GetNumProductTemplates()):
		mol = rxn.GetProductTemplate(j)
		mol.UpdatePropertyCache(False)
		mols.append(mol)

	# Generate images for all molecules/arrow
	imgs = [TrimImgByWhite(MolToImage(mol), padding = 15) for mol in mols]
	for i, img in enumerate(imgs):
		img.save('{}.png'.format(i))

	# Combine
	return StitchPILsHorizontally(imgs)

def ReactionStringToImage(rxn_string, strip = False, **kwargs):
	'''This function takes a SMILES rxn_string as input, not an 
	RDKit reaction object, and draws it.'''

	reactants, agents, products = [mols_from_smiles_list(x) for x in 
										[mols.split('.') for mols in rxn_string.split('>')]]
	if None in reactants + agents + products:
		raise ValueError('Could not parse entirety of reaction: {}'.format(rxn_string))

	# Stich together mols (ignore agents)
	mols = reactants + [None] + products
	if strip: mols = [Chem.MolFromInchi(Chem.MolToInchi(mol)) for mol in mols]

	# Generate images
	imgs = [TrimImgByWhite(MolToImage(mol), padding = 15) for mol in mols]

	# Combine
	return StitchPILsHorizontally(imgs)


mol = Chem.MolFromSmiles('[C:1][C:4]OCCN')
fpath = 'test_mol.png'
img = MolToImage(mol)
img.save('test.png')
img_cropped = TrimImgByWhite(img, padding = 15)
img_cropped.save('test_crop.png')
rxn_string = '[C:1]/[C:2]=[C:3]/[C:5][C:4]>>[C:1][C:2]\[C:3]=[C:4]\[C:5]'
rxn = AllChem.ReactionFromSmarts(rxn_string)
rxn_image = ReactionToImage(rxn)
rxn_image.save('test_rxn.png')

rxn_image_from_string = ReactionStringToImage(rxn_string)
rxn_image_from_string.save('test_rxn_string.png')