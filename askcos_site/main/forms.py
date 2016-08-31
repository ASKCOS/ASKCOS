from django import forms
from rdkit.Chem import MolFromSmiles


def is_valid_smiles(smiles):
	return (MolFromSmiles(smiles) != None)


class SmilesInputForm(forms.Form):
	smiles = forms.CharField(label = 'Enter SMILES string',
		max_length = 150)

class DrawingInputForm(forms.Form):
	text = forms.CharField(label = 'Enter text to draw', 
		max_length = 800)
	# style = forms.MultipleChoiceField(
	# 	required = True,
	# 	widget = forms.RadioSelect,
	# 	choices = (('mol', 'Molecule SMILES'), ('rxn', 'Reaction SMILES'), ('tform', 'Template SMARTS'))
	# )
