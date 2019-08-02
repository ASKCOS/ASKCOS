import gzip
import json
import copy

from rdkit import Chem

with gzip.open('buyables.all.json.gz', 'rb') as f:
    content = f.read().decode('utf-8')
    buyables_all = json.loads(content)

buyables = {}
for doc in buyables_all:
    for smiles in doc['smiles'].split('.'):
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
        if buyables.get(smiles):
            buyables[smiles]['ppg'] = min(doc['ppg'], buyables[smiles]['ppg'])
        else:
            new_doc = copy.deepcopy(doc)
            new_doc.pop('_id', '')
            new_doc['smiles'] = smiles
            buyables[smiles] = new_doc

with gzip.open('buyables.json.gz', 'wt') as f:
    json.dump(list(buyables.values()), f)