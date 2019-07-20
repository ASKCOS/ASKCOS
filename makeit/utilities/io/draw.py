from __future__ import absolute_import
import numpy as np
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.Draw as Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry
from PIL import Image, ImageOps
from collections import defaultdict
# from rdkit.Chem.Draw.cairoCanvas import Canvas
import os
import re


def get_scaled_drawer(mol):
    #Draw all molecules with same proportions
    dpa = 26
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    point_min = Geometry.rdGeometry.Point2D()
    point_max = Geometry.rdGeometry.Point2D()
    point_min.x = min(xs) - 1
    point_min.y = min(ys) - 1
    point_max.x = max(xs) + 1
    point_max.y = max(ys) + 1
    w = int(dpa * (point_max.x - point_min.x))
    h = int(dpa * (point_max.y - point_min.y))
    drawer = rdMolDraw2D.MolDraw2DCairo(w,h)
    drawer.SetScale(w, h, point_min, point_max)
    return drawer

def MolsSmilesToImageHighlight(smiles, options=None, **kwargs):
    '''This function takes a SMILES string of one or more molecules
    and generates a combined image for that molecule set.'''
    mol = Chem.MolFromSmiles(str(smiles))
    d2 = get_scaled_drawer(mol)
    dopts = d2.drawOptions()
    reacting_atoms = kwargs.get('reacting_atoms', [])
    bonds = kwargs.get('bonds', False)
    if len(reacting_atoms) != 0:
        #highlightAtoms = list(range(Chem.MolFromSmiles(smiles).GetNumAtoms()))
        highlightAtoms = set([x for tup in ra for x in tup])
        highlightAtomColors = {x:(0,1,0) for x in highlightAtoms}
        if bonds:
            #TODO some edits have multiple atoms changing in one tuple
            highlightBonds = [mol.GetBondBetweenAtoms(x[0],x[1]).GetIdx() for x in reacting_atoms]
        else:
            highlightBonds = []
    else:
        atomScores = []
        highlightAtoms = []
        highlightBonds = []
        highlightAtomColors = []
    
    m2=Draw.PrepareMolForDrawing(mol)
    d2.DrawMolecule(m2,highlightAtoms=highlightAtoms, \
        highlightBonds=highlightBonds, highlightAtomColors=highlightAtomColors)
    d2.FinishDrawing()
    txt = d2.GetDrawingText()
    
    return txt
 

'''
Many of these functions are taken from RDKit.
'''


def mols_from_smiles_list(all_smiles):
    '''Given a list of smiles strings, this function creates rdkit
    molecules'''
    mols = []
    for smiles in all_smiles:
        if not smiles:
            continue
        mols.append(Chem.MolFromSmiles(smiles))
    return mols


def defaultDrawOptions():
    '''This function returns an RDKit drawing options object with 
    default drawing options.'''

    opts = Draw.DrawingOptions()
    # opts.elemDict = defaultdict(lambda: (0,0,0)) # all atoms are black
    opts.noCarbonSymbols = True
    opts.selectColor = (1, 0, 0)
    opts.wedgeBonds = True
    
    opts.elemDict = defaultdict(lambda: (0, 0, 0))
    opts.dotsPerAngstrom = 20
    opts.bondLineWidth = 1.5
    atomLabelFontFace = 'arial'

    return opts


def StripAlphaFromImage(img):
    '''This function takes an RGBA PIL image and returns an RGB image'''

    if len(img.split()) == 3:
        return img
    return Image.merge('RGB', img.split()[:3])


def MolToImage(mol, max_size=(1000, 1000), kekulize=True, options=None,
               canvas=None, **kwargs):
    '''Wrapper for RDKit's MolToImage. If mol == None, an arrow is drawn'''

    if not options:
        options = defaultDrawOptions()
    if mol == '->':
        subImgSize = (100, 100)
        img, canvas = Draw._createCanvas(subImgSize)
        p0 = (10, subImgSize[1]//2)
        p1 = (subImgSize[0]-10, subImgSize[1]//2)
        p3 = (subImgSize[0]-20, subImgSize[1]//2-10)
        p4 = (subImgSize[0]-20, subImgSize[1]//2+10)
        canvas.addCanvasLine(p0, p1, lineWidth=2, color=(0, 0, 0))
        canvas.addCanvasLine(p3, p1, lineWidth=2, color=(0, 0, 0))
        canvas.addCanvasLine(p4, p1, lineWidth=2, color=(0, 0, 0))
        if hasattr(canvas, 'flush'):
            canvas.flush()
        else:
            canvas.save()
        return img        
    elif mol == '<-':  # retro arrow or error
        subImgSize = (100, 100)
        (a, b) = subImgSize
        img, canvas = Draw._createCanvas(subImgSize)
        canvas.addCanvasLine((10, b//2-7), (a-17, b//2-7),
                             lineWidth=2, color=(0, 0, 0))
        canvas.addCanvasLine((10, b//2+7), (a-17, b//2+7),
                             lineWidth=2, color=(0, 0, 0))
        canvas.addCanvasLine((a-24, b//2-14), (a-10, b//2),
                             lineWidth=2, color=(0, 0, 0))
        canvas.addCanvasLine((a-24, b//2+14), (a-10, b//2),
                             lineWidth=2, color=(0, 0, 0))
        if hasattr(canvas, 'flush'):
            canvas.flush()
        else:
            canvas.save()
        return img
    elif mol is not None:
        return Draw.MolToImage(mol, size=max_size, kekulize=kekulize, options=options,
                               canvas=canvas, **kwargs)


def TrimImgByWhite(img, padding=0):
    '''This function takes a PIL image, img, and crops it to the minimum rectangle 
    based on its whiteness/transparency. 5 pixel padding used automatically.'''

    # Convert to array
    as_array = np.array(img)  # N x N x (r,g,b,a)

    # Set previously-transparent pixels to white
    if as_array.shape[2] == 4:
        as_array[as_array[:, :, 3] == 0] = [255, 255, 255, 0]

    as_array = as_array[:, :, :3]

    # Content defined as non-white and non-transparent pixel
    has_content = np.sum(as_array, axis=2, dtype=np.uint32) != 255 * 3
    xs, ys = np.nonzero(has_content)

    # Crop down
    margin = 5
    x_range = max([min(xs) - margin, 0]), min([max(xs) + margin, as_array.shape[0]])
    y_range = max([min(ys) - margin, 0]), min([max(ys) + margin, as_array.shape[1]])
    as_array_cropped = as_array[
        x_range[0]:x_range[1], y_range[0]:y_range[1], 0:3]

    img = Image.fromarray(as_array_cropped, mode='RGB')

    return ImageOps.expand(img, border=padding, fill=(255, 255, 255))


def StitchPILsHorizontally(imgs):
    '''This function takes a list of PIL images and concatenates
    them onto a new image horizontally, with each one
    vertically centered.'''

    # Create blank image (def: transparent white)
    heights = [img.size[1] for img in imgs]
    height = max(heights)
    widths = [img.size[0] for img in imgs]
    width = sum(widths)
    res = Image.new('RGB', (width, height), (255, 255, 255))

    # Add in sub-images
    for i, img in enumerate(imgs):
        offset_x = sum(widths[:i])  # left to right
        offset_y = (height - heights[i]) // 2
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
    elif '[C:' in smarts and 'H' not in smarts:
        atom.SetAtomicNum(0)
        atom.SetProp('dummyLabel', 'C[al]')
        atom.UpdatePropertyCache()
    elif '[c:' in smarts and 'H' not in smarts:
        atom.SetAtomicNum(0)
        atom.SetProp('dummyLabel', 'C[ar]')
        atom.UpdatePropertyCache()

    # Clear atom map number of 0 -> this is a dummy assignment!
    if ':0]' in smarts:
        atom.ClearProp('molAtomMapNumber')


def ReactionToImage(rxn, dummyAtoms=False, kekulize=True, options=None, **kwargs):
    '''Modification of RDKit's ReactionToImage to allow for each molecule 
    to have a different drawn size. rxn is an RDKit reaction object

    warning: this function adds hydrogens as it sees fit'''
    # Extract mols from reaction
    mols = []
    for i in range(rxn.GetNumReactantTemplates()):
        mol = rxn.GetReactantTemplate(i)
        mol.UpdatePropertyCache(False)
        mols.append(mol)
        if dummyAtoms:
            [CheckAtomForGeneralization(atom) for atom in mol.GetAtoms()]

    if kwargs.pop('retro', True):
        mols.append('<-')  # placeholder for arrow
    else:
        mols.append('->')

    for j in range(rxn.GetNumProductTemplates()):
        mol = rxn.GetProductTemplate(j)
        mol.UpdatePropertyCache(False)
        mols.append(mol)
        if dummyAtoms:
            [CheckAtomForGeneralization(atom) for atom in mol.GetAtoms()]

    # Generate images for all molecules/arrow
    imgs = [TrimImgByWhite(MolToImage(
        mol, kekulize=kekulize, options=options), padding=10) for mol in mols]

    # Combine
    return StitchPILsHorizontally(imgs)


def ReactionStringToImage(rxn_string, strip=True, update=True, options=None,
        retro=False, **kwargs):
    '''This function takes a SMILES rxn_string as input, not an 
    RDKit reaction object, and draws it.'''

    reactants, agents, products = [mols_from_smiles_list(x) for x in
                                   [mols.split('.') for mols in rxn_string.split('>')]]
    if None in reactants + products:
        raise ValueError(
            'Could not parse entirety of reaction: {}'.format(rxn_string))

    # Stich together mols (ignore agents)
    if retro:
        mols = reactants + ['<-'] + products
    else:
        mols = reactants + ['->'] + products
    if update:
        [mol.UpdatePropertyCache(False) for mol in mols if mol is not None and type(mol) != str]
    if strip:
        for mol in mols:
            if mol is not None and type(mol) != str:
                [a.ClearProp('molAtomMapNumber') for a in mol.GetAtoms()]

    # Generate images
    imgs = [TrimImgByWhite(MolToImage(
        mol, kekulize=True, options=options), padding=10) for mol in mols]

    # Combine
    return StitchPILsHorizontally(imgs)


def TransformStringToImage(transform, retro=True, **kwargs):
    '''Wrapper function meant to take a SMARTS transform and return a PIL image
    of that transform.

    TODO: Need to improve generalization visually! Right now it still shows'''

    options = defaultDrawOptions()
    options.dotsPerAngstrom = 40

    # To generalize un-mapped atoms in transform, need to identify square brackets
    # without colon in the middle (e.g., [C]) and replace with dummy label [C:0] so
    # generalization display works
    old_tags = re.findall('\[[^:]+\]', transform)
    for old_tag in old_tags:
        new_tag = old_tag.replace(']', ':0]')
        transform = transform.replace(old_tag, new_tag)
    rxn = AllChem.ReactionFromSmarts(transform)
    return ReactionToImage(rxn, dummyAtoms=True, options=options, retro=retro, **kwargs)


def MolsSmilesToImage(smiles, options=None, **kwargs):
    '''This function takes a SMILES string of one or more molecules
    and generates a combined image for that molecule set.'''

    # Generate mols
    mols = mols_from_smiles_list(smiles.split('.'))
    # Generate images
    imgs = [TrimImgByWhite(MolToImage(
        mol, kekulize=True, options=options), padding=10) for mol in mols]
    # Combine
    return StitchPILsHorizontally(imgs)


def main():
    # Simple test cases
    rxn_string = 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1.c1nc[nH]n1.Cl.O=C([O-])O.[Na+]>>OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F'
    # rxn = AllChem.ReactionFromSmarts(rxn_string)
    # rxn_image = ReactionToImage(rxn)
    # rxn_image.save('test_rxn.png')
    rxn_image_string = ReactionStringToImage(rxn_string, strip=True)
    rxn_image_string.save('draw_test_rxn_string.png')
    rxn_image_string = ReactionStringToImage(rxn_string, strip=True, retro=True)
    rxn_image_string.save('draw_retro_test_rxn_string.png')

    tform = '([O;H0:3]=[C;H0:4](-[C:5])-[NH:2]-[C:1])>>([C:1]-[NH2:2]).([OH:3]-[C;H0:4](=O)-[C:5])'
    img = TransformStringToImage(tform)
    img.save('draw_transform.png')

if __name__ == '__main__':
    main()
