# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
INSTANCE_DB = db['instances']


print('Num that have yield, T, at least one solvent, and time: {}'.format(
	INSTANCE_DB.count({
			'RXD_NYD': {'$ne': -1}, 
			'RXD_T': {'$ne': -1}, 
			'RXD_SOLXRN': {'$ne': []},
			'RXD_TIM': {'$ne': -1},
		}, no_cursor_timeout = True)
	))

INSTANCE_DB.update({
	'RXD_NYD': {'$ne': -1}, 
	'RXD_T': {'$ne': -1}, 
	'RXD_SOLXRN': {'$ne': []},
	'RXD_TIM': {'$ne': -1},
}, {'$set': {'complete': True}}, multi = True)

print('{} instances are marked as complete'.format(INSTANCE_DB.count({'complete': True})))