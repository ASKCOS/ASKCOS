def chem_dict(_id, smiles, ppg, children = []):
    '''Chemical object as expected by website'''
    return {
        'id': _id,
        'is_chemical': True,
        'smiles' : smiles,
        'ppg' : ppg,
        'children': children,}
    
def rxn_dict(_id, info, template_score =1, necessary_reagent = '', num_examples = 0, children = [], smiles = ''):
    '''Reaction object as expected by website'''
    return {
        'id': _id,
        'is_reaction': True,
        'info': info,
        'necessary_reagent': necessary_reagent,
        'num_examples': num_examples,
        'children': children,
        'smiles': smiles,
        'template_score':template_score
        }
        