import pymongo
from pymongo import MongoClient

global db_client
try:
	db_client
except NameError:
	db_client = MongoClient('mongodb://guest:guest@rmg.mit.edu:27017/admin', 
		serverSelectionTimeoutMS = 2000,
		connect=True)
	print('Pymongo version {}'.format(pymongo.version))
	print('Initialized Mongo connection, DB names follow:')
	print(db_client.database_names())
