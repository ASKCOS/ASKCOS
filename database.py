import pymongo
from pymongo import MongoClient

global db_client
try:
    try:
        db_client
    except NameError:
        db_client = MongoClient('mongodb://guest:guest@askcos2.mit.edu:27017/admin', 
            serverSelectionTimeoutMS = 2000,
            connect=True)
except Exception as e:
    db_client = None
    print('##### DB COULD NOT CONNECT ####')
    # print('Pymongo version {}'.format(pymongo.version))
    # print('Initialized Mongo connection, DB names follow:')
    # print(db_client.database_names())