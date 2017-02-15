from pymongo import MongoClient

client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['askcos_transforms']
SMARTS_DB = db['EFGs']
SMARTS_DB.remove({})

from EFGs import *

for i, fgroup in enumerate(fgroups):
	fgroups[i]['_id'] = i 

SMARTS_DB.insert(fgroups)