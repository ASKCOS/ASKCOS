def chem_dict(_id, children=[], **kwargs):
    '''Chemical object as expected by website. Removes rct_of, prod_of, and 
    depth information from the tree_dict entry'''
    kwargs.pop('rct_of', None)
    kwargs.pop('prod_of', None)
    kwargs.pop('depth', None)
    kwargs['id'] = _id 
    kwargs['is_chemical'] = True 
    kwargs['children'] = children 
    return kwargs 


def rxn_dict(_id, smiles, children=[], **kwargs):
    '''Reaction object as expected by website. Removes depth and rcts, prod, depth'''
    kwargs.pop('rcts', None)
    kwargs.pop('prod', None)
    kwargs.pop('depth', None)
    kwargs['id'] = _id
    kwargs['is_reaction'] = True 
    kwargs['children'] = children
    kwargs['smiles'] = smiles
    return kwargs