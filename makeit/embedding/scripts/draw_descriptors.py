from rdkit.Chem.Draw import SimilarityMaps
import rdkit.Chem as Chem
from makeit.embedding.descriptors import *

mol = Chem.MolFromSmiles('Fc1ccccc1')
(labels, attributes) = atom_level_descriptors(mol)
selected_attribute = [x[-2] for x in attributes]

# for i, atom in enumerate(mol.GetAtoms()):
# 	atom.SetProp('molAtomMapNumber', str(int(selected_attribute[i] * 1000)))

fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, selected_attribute, 
	size = (400, 400), colorMap='jet', contourLines=10)
fig.savefig('test.png', bbox_inches='tight')