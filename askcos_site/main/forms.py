from django import forms
from rdkit.Chem import MolFromSmiles
from django.forms import modelform_factory
from models import ConditionPredictorSetup, SeparationInput

def is_valid_smiles(smiles):
    return (MolFromSmiles(smiles) != None)


class SmilesInputForm(forms.Form):
    smiles = forms.CharField(label = 'Enter SMILES string',
        max_length = 400)

class DrawingInputForm(forms.Form):
    text = forms.CharField(label = 'Enter text to draw',
        max_length = 1600)
    # style = forms.MultipleChoiceField(
    # 	required = True,
    # 	widget = forms.RadioSelect,
    # 	choices = (('mol', 'Molecule SMILES'), ('rxn', 'Reaction SMILES'), ('tform', 'Template SMARTS'))
    # )

nnSetup = modelform_factory(ConditionPredictorSetup, exclude=('N_conditions',))
sepInput = modelform_factory(SeparationInput, exclude=('PHASE_CHOICES', 'OPT_TGT_CHOICES'))