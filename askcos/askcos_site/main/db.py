import gzip
import json
import pymongo
import makeit.global_config as gc

db_client = pymongo.MongoClient(gc.MONGO['path'], gc.MONGO[ 'id'], connect=gc.MONGO['connect'])

def seed_mongo_db(buyables=True, chemicals=True, reactions=True, retro_templates=True, forward_templates=True):
    if buyables:
        seed_mongo_documents('buyables', gc.BUYABLES['file_name'])
        db_client[gc.BUYABLES['database']][gc.BUYABLES['collection']].create_index([('smiles', pymongo.TEXT)])
    if chemicals:
        seed_mongo_documents('chemicals', gc.CHEMICALS['file_name'])
        db_client[gc.CHEMICALS['database']][gc.CHEMICALS['collection']].create_index([('smiles', pymongo.HASHED)])
    if reactions:
        seed_mongo_documents('reactions', gc.REACTIONS['file_name'])
    if retro_templates:
        seed_mongo_documents('retro_templates', gc.RETRO_TEMPLATES['file_name'])
        db_client[gc.RETRO_TEMPLATES['database']][gc.RETRO_TEMPLATES['collection']].create_index([('index', pymongo.ASCENDING)])
    if forward_templates:
        seed_mongo_documents('forward_templates', gc.FORWARD_TEMPLATES['file_name'])
    
def seed_mongo_documents(collection_name, file_path, db_name='askcos', clear=True, bypass_document_validation=True):
    coll = db_client[db_name][collection_name]
    if clear:
        res = coll.delete_many({})
    
    with gzip.open(file_path, 'rb') as f:
        docs = json.loads(f.read().decode('utf-8'))
    
    res = coll.insert_many(docs, ordered=False, bypass_document_validation=bypass_document_validation)
    print('Inserted {} documents into {}.{}'.format(len(res.inserted_ids), db_name, collection_name))
