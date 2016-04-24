from django import forms

class SmilesInputForm(forms.Form):
	smiles = forms.CharField(label = 'Enter SMILES string',
		max_length = 150)