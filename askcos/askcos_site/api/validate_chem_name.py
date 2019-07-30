from django.http import JsonResponse
from rdkit import Chem

def validate_chem_name(request):
    '''Check the syntax and validity of a SMILES string
    Input:
    smiles:             string,             SMILES string

    Output:
    correct_syntax:     bool,               correctness of the SMILES syntax
    valid_chem_name:    bool,               validity of the chemical

    RDKit ref: https://github.com/rdkit/rdkit/issues/2430
    '''
    resp = {}
    resp['request'] = dict(**request.GET)
    smiles = request.GET.get('smiles', None)

    # results
    correct_syntax = None
    valid_chem_name = None

    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        correct_syntax = False
    else:
        correct_syntax = True
        try:
            Chem.SanitizeMol(m)
            valid_chem_name = True
        except:
            valid_chem_name = False

    resp['correct_syntax'] = correct_syntax
    resp['valid_chem_name'] = valid_chem_name
    return JsonResponse(resp)