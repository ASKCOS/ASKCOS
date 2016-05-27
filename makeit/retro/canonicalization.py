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
			# Double bonds on aromatic rings (dist 1)
			AllChem.ReactionFromSmarts('[NH0:1]=[a:2][nH:3]>>[NH:1][a:2][nH0:3]'),
			AllChem.ReactionFromSmarts('[NH1:1]=[a:2][nH:3]>>[NH2:1][a:2][nH0:3]'),
			AllChem.ReactionFromSmarts('[OH0:1]=[a:2][nH:3]>>[OH:1][a:2][nH0:3]'),
			# Double bonds on aromatic rings (dist 2)
			AllChem.ReactionFromSmarts('[NH0:1]=[a:2][a:4][nH:3]>>[NH:1][a:2][a:4][nH0:3]'),
			AllChem.ReactionFromSmarts('[NH1:1]=[a:2][a:4][nH:3]>>[NH2:1][a:2][a:4][nH0:3]'),
			AllChem.ReactionFromSmarts('[OH0:1]=[a:2][a:4][nH:3]>>[OH:1][a:2][a:4][nH0:3]'),
			# Double bonds on aromatic rings (dist 2)
			AllChem.ReactionFromSmarts('[NH0:1]=[a:2][a:4][a:5][nH:3]>>[NH:1][a:2][a:4][a:5][nH0:3]'),
			AllChem.ReactionFromSmarts('[NH1:1]=[a:2][a:4][a:5][nH:3]>>[NH2:1][a:2][a:4][a:5][nH0:3]'),
			AllChem.ReactionFromSmarts('[OH0:1]=[a:2][a:4][a:5][nH:3]>>[OH:1][a:2][a:4][a:5][nH0:3]'),
			# Iminol / amide 
			AllChem.ReactionFromSmarts('[NH0:1]=[C:2]-[OH:3]>>[NH1:1]-[C:2]=[OH0:3]'),
			AllChem.ReactionFromSmarts('[NH1:1]=[C:2]-[OH:3]>>[NH2:1]-[C:2]=[OH0:3]'),
			# Thiourea
			AllChem.ReactionFromSmarts('[NH0:1]=[C:2]-[SH:3]>>[NH1:1]-[C:2]=[SH0:3]'),
			AllChem.ReactionFromSmarts('[NH1:1]=[C:2]-[SH:3]>>[NH2:1]-[C:2]=[SH0:3]'),
			# Azide
			AllChem.ReactionFromSmarts('[NH0:1][NH0:2]=[NH0;-:3]>>[NH0;-:1]=[NH0;+:2]=[NH0;-:3]'),
			# Cyanide salts
			AllChem.ReactionFromSmarts('([K,Na;H1:1].[C;X1;H0:2]#[N:3])>>[*;H0:1][*:2]#[N:3]'),
			AllChem.ReactionFromSmarts('([Cu:1].[C;X1;H0:2]#[N:3])>>[*:1][*:2]#[N:3]'),
			# Grinards
			AllChem.ReactionFromSmarts('([MgH+:1].[C;v3:2][*:3])>>[Mg+:1][*:2][*:3]')
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