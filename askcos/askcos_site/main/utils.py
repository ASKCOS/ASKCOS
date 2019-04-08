import os
from django.http import JsonResponse
from .globals import CHEMICAL_DB
import rdkit.Chem as Chem
import sys
if sys.version_info[0] < 3:
    from urllib2 import urlopen
    from urllib2 import HTTPError
else:
    from urllib.request import urlopen
    from urllib.error import HTTPError

from askcos_site.celery import app

def ajax_error_wrapper(ajax_func):
    def ajax_func_call(*args, **kwargs):
        try:
            return ajax_func(*args, **kwargs)
        except Exception as e:
            print(e)
            return JsonResponse({'err':True, 'message': str(e)})

    return ajax_func_call

def fancyjoin(lst, nonemessage='(none)'):
    if not lst:
        return nonemessage
    if len(lst) == 1:
        return lst[0]
    if len(lst) == 2:
        return '%s and %s' % (lst[0], lst[1])
    return ', '.join(lst[:-1]) + ', and %s' % lst[-1]

def xrn_lst_to_name_lst(xrn_lst):
    lst = []
    xrn_to_smiles = {}
    for xrn in xrn_lst:
        if xrn not in xrn_to_smiles:
            chem_doc = CHEMICAL_DB.find_one({'_id': xrn})
            if chem_doc is None:
                xrn_to_smiles[xrn] = 'Chem-%i' % xrn
            elif 'IDE_CN' not in chem_doc:
                if 'SMILES' not in chem_doc:
                    xrn_to_smiles[xrn] = 'Chem-%i' % xrn
                else:
                    xrn_to_smiles[xrn] = chem_doc['SMILES']
            else:
                xrn_to_smiles[xrn] = chem_doc['IDE_CN']
        lst.append(xrn_to_smiles[xrn])
    return lst

def resolve_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        # Try to resolve using NIH
        new_smiles = []
        for smiles in smiles.split(' and '):
            try:
                smiles = urlopen('https://cactus.nci.nih.gov/chemical/structure/{}/smiles'.format(smiles)).read()
            except HTTPError:
                return None
            mol = Chem.MolFromSmiles(smiles)
            if not mol: return None
            new_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        return '.'.join(new_smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def get_name_from_smiles(smiles):
    try:
        names = urlopen('https://cactus.nci.nih.gov/chemical/structure/{}/names'.format(smiles)).read()
        if '&lt;!DOCTYPE html&gt;' in names:
            return smiles
        return names.split('\n')[0]
    except HTTPError:
        return smiles

def fix_rgt_cat_slvt(rgt1, cat1, slvt1):
    # Merge cat and reagent for forward predictor
    if rgt1 and cat1:
        rgt1 = rgt1 + '.' + cat1
    elif cat1:
        rgt1 = cat1
    # Reduce solvent to single one
    if '.' in slvt1:
        slvt1 = slvt1.split('.')[0]
    return (rgt1, cat1, slvt1)

def trim_trailing_period(txt):
    if txt:
        if txt[-1] == '.':
            return txt[:-1]
    return txt

def restrict(view_function):
    def authorize(request, *args, **kwargs):
        RESTRICT_API = os.getenv('RESTRICT_API', False)
        request_ip = request.META['REMOTE_ADDR']
        first_block = request_ip.split('.')[0]
        if not RESTRICT_API or first_block.startswith('18'):
            return view_function(request, *args, **kwargs)
        return JsonResponse({'error': 'You are not authorized to access this endpoint. Contact admin for more information'})
    return authorize
