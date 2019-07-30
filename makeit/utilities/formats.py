def chem_dict(_id, children=[], **kwargs):
    """Returns chemical dictionary in the format required by the website.

    Chemical object as expected by website. Removes ``rct_of``, ``prod_of``, and
    ``depth`` information from the ``tree_dict`` entry.

    Args:
        _id (int): Chemical ID.
        children (list, optional): Children of the node.
        **kwargs: The ``tree_dict`` to be modified.
    """
    kwargs.pop('rct_of', None)
    kwargs.pop('prod_of', None)
    kwargs.pop('depth', None)
    kwargs['id'] = _id
    kwargs['is_chemical'] = True
    kwargs['children'] = children
    return kwargs


def rxn_dict(_id, smiles, children=[], **kwargs):
    """Returns reaction dictionary in the format required by the website.

    Reaction object as expected by website. Removes ``rct``, ``prod``, and
    ``depth`` information from the ``tree_dict`` entry.

    Args:
        _id (int): Chemical ID.
        smiles (str): SMILES string of reaction.
        children (list, optional): Children of the node.
        **kwargs: The ``tree_dict`` to be modified.
    """
    kwargs.pop('rcts', None)
    kwargs.pop('prod', None)
    kwargs.pop('depth', None)
    kwargs['id'] = _id
    kwargs['is_reaction'] = True
    kwargs['children'] = children
    kwargs['smiles'] = smiles
    return kwargs
