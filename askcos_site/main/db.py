from pymongo import MongoClient

global db_client
try:
	db_client
except NameError:
	import makeit.global_config as gc
	db_client = MongoClient(gc.MONGO['path'], gc.MONGO[ 'id'], connect=gc.MONGO['connect'])