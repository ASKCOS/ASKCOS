from django import forms
from rdkit.Chem import MolFromSmiles

class SmilesInputForm(forms.Form):
	smiles = forms.CharField(label = 'Enter SMILES string',
		max_length = 150)

def is_valid_smiles(smiles):
	return (MolFromSmiles(smiles) != None)
