import os, sys, json, urllib
import rdkit.Chem as Chem 
from makeit.utilities.io.name_parser import name_to_molecule, urlopen

## Syntax of banned list:
# json-dumped dictionary of "name": [isomeric_smiles, flat_smiles]
# where smiles could be None

# Where is the banned chemical list stored?
banned_names_fpath = os.path.join(os.path.dirname(__file__), 'banned_names.txt')
banned_list_fpath = os.path.join(os.path.dirname(__file__), 'banned_list.json')

# Default
banned = {}

# Open list
try:
    with open(banned_list_fpath, 'r') as fid:
        banned = json.load(fid)
except:
    print('Warning: did not find existing banned_list.json file')
    pass

# Check that all names are in the banned_list
try:
    with open(banned_names_fpath, 'r') as fid:
        names = [line.split('#')[0].strip() for line in fid.readlines()]
    for name in names:
        if not name:
            continue
        if name not in banned:
            banned[name] = (None, None)
except:
    print('Warning: did not find banned_names.txt file')
    pass

# Try to fill in SMILES that are missing automatically
for (name, smis) in banned.items():
    
    if smis[0] is None: # use NIH resolver
        escaped_name = urllib.parse.quote(name)
        try:
            mol = name_to_molecule(escaped_name)
            if mol is None:
                continue

            banned[name] = (
                Chem.MolToSmiles(mol, True),
                Chem.MolToSmiles(mol, False),
            )
            print('Parsed {} --> {}'.format(name, banned[name]))
        except:
            pass
    
    else: # try to canonicalize
        mol = Chem.MolFromSmiles(smis[0])
        banned[name] = (
            Chem.MolToSmiles(mol, True),
            Chem.MolToSmiles(mol, False),
        )


# Resave list
with open(banned_list_fpath, 'w') as fid:
    json.dump(banned, fid, indent=4, sort_keys=True)
