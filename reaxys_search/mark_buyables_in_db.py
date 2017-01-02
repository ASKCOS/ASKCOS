from pymongo import MongoClient
from bson.objectid import ObjectId
from tqdm import tqdm 

'''This script is meant to take the 'buyables' database and match products against
the 'chemicals' database to mark them.'''


DEBUG = False

client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
CHEMICAL_DB = db['chemicals']
BUYABLES_DB = db['buyables']

print('Found {} invalid entries'.format(CHEMICAL_DB.count({'SMILES': '', 'buyable_id': {'$exists': True}})))
CHEMICAL_DB.update(
	{'SMILES': '', 'buyable_id': {'$exists': True}},
	{'$unset': {'buyable_id': 1}},
	multi = True 
	)
print('Found {} invalid entries'.format(CHEMICAL_DB.count({'SMILES': '', 'buyable_id': {'$exists': True}})))


for i, buyable_doc in enumerate(BUYABLES_DB.find(no_cursor_timeout = True)):
	if 'smiles' not in buyable_doc: continue
	smiles = buyable_doc['smiles']
	if not smiles: 
		print('### empty smiles, skip')
		continue 
	print('{}: {}'.format(i, smiles))

	chemical_doc = None
	for chemical_doc in CHEMICAL_DB.find({'SMILES': smiles}):
		CHEMICAL_DB.update_one(
			{'_id': chemical_doc['_id']},
			{'$set': 
				{'buyable_id': buyable_doc['_id']}
			}
		)
		if DEBUG: print('Updated chemical ID {} to match buyable ID {}'.format(chemical_doc['_id'], buyable_doc['_id']))

	if DEBUG:
		if chemical_doc == None:
			print('Could not find a match for #{} - {}'.format(buyable_doc['_id'], smiles))
		else:
			print('Found a match for #{} - {}'.format(buyable_doc['_id'], smiles))
		raw_input('Pause...')