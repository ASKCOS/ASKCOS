'''
This script moves templates from one database to another. While doing so,
it removes wildcards that are (hopefully?) unnecessary. This will make matching
go faster. Hopefully by a factor of ~2
'''

import re
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from tqdm import tqdm 

replacements = [
    (r'^\[\*\:[0-9]+\][-=~#]?', ''),        # leading wildcard
    (r'[-=~#]?\[\*\:[0-9]+\]$', ''),        # ending wildcard
    (r'\([-=~#]?\[\*\:[0-9]+\]\)', ''),    # alone in a side branch, w/ or w/o bond specified
    (r'[-=~#]?\[\*\:[0-9]+\]\)', ')'),     # last in a side branch, w/ or w/o bond specified
]

def clean_fragment(fragment):
    for pattern, sub in replacements:
        fragment = re.sub(pattern, sub, fragment)
    return fragment

if __name__ == '__main__':
    # DATABASE
    from pymongo import MongoClient    # mongodb plugin
    client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
    db = client['reaxys']
    original_templates = db['transforms_forward_v1']
    new_templates = db['transforms_forward_v1_nowildcards']

    docs_to_insert = []
    for template_doc in tqdm(original_templates.find()):
        reaction_smarts = template_doc['reaction_smarts']

        reactants, products = reaction_smarts.split('>>')
        reactants = '.'.join([clean_fragment(x) for x in reactants.split('.')])
        products = '.'.join([clean_fragment(x) for x in products.split('.')])
        reaction_smarts_new = reactants + '>>' + products

        #print('Original SMARTS: {}'.format(reaction_smarts))
        #print('Proposed SMARTS: {}'.format(reaction_smarts_new))
        
        template_doc['reaction_smarts'] = reaction_smarts_new
        docs_to_insert.append(template_doc)

        if len(docs_to_insert) == 1000:
            new_templates.insert(docs_to_insert)
            del docs_to_insert 
            docs_to_insert = []
