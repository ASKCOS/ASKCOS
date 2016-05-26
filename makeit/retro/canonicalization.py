import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from global_config import USE_STEREOCHEMISTRY

class SmilesFixer():
	'''
	This class stores RDKit reactions which help turn molecules with
	weird representations into ones with common ones (so they will match
	a template)
	'''
	def __init__(self):
		self.rxns = [
			# Double bonds on aromatic rings
			AllChem.ReactionFromSmarts('[NH0:1]=[c:2][nH:3]>>[NH:1][c:2][nH0:3]'),
			AllChem.ReactionFromSmarts('[NH1:1]=[c:2][nH:3]>>[NH2:1][c:2][nH0:3]'),
			# Iminol / amide 
			AllChem.ReactionFromSmarts('[NH0:1]=[C:2]-[OH:3]>>[NH1:1]-[C:2]=[OH0:3]'),
			AllChem.ReactionFromSmarts('[NH1:1]=[C:2]-[OH:3]>>[NH2:1]-[C:2]=[OH0:3]')
		]

	def fix_smiles(self, old_smiles):
		'''
		For a given SMILES string, this function "fixes" common mistakes
		found in the Lowe parsed database:
		- N=c[nH] structures are turned into the normal [NH]-c[n] forms
		- iminols are turned into amides/carbamates

		It applies the reactions in self.rxns until the SMILES string doesn't change
		'''
		mol = Chem.MolFromSmiles(old_smiles)
		if not mol: 
			return old_smiles 

		new_smiles = old_smiles
		old_smiles = ''
		while new_smiles != old_smiles:
			old_smiles = new_smiles
			for rxn in self.rxns:
				outcomes = rxn.RunReactants((mol,))
				if not outcomes: 
					continue
				else:
					mol = outcomes[0][0]
					Chem.SanitizeMol(mol)
					new_smiles = Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY)

		return new_smiles