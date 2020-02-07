import gzip
import json
import pymongo
import makeit.global_config as gc

db_client = pymongo.MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect=gc.MONGO['connect'])


def seed_mongo_db(buyables=True, chemicals=True, reactions=True, retro_templates=True, forward_templates=True):
    """Seed mongo database with default data files as defined in the makeit/global_config.py file."""
    if buyables:
        seed_buyables(gc.BUYABLES['file_name'])
    if chemicals:
        seed_chemicals(gc.CHEMICALS['file_name'])
    if reactions:
        seed_reactions(gc.REACTIONS['file_name'])
    if retro_templates:
        seed_retro_templates(gc.RETRO_TEMPLATES['file_name'])
    if forward_templates:
        seed_forward_templates(gc.FORWARD_TEMPLATES['file_name'])


def seed_buyables(filename):
    """Use the specified file to seed the buyables database"""
    seed_mongo_documents('buyables', filename)
    db_client[gc.BUYABLES['database']][gc.BUYABLES['collection']].create_index([('smiles', pymongo.TEXT)])


def seed_chemicals(filename):
    """Use the specified file to seed the chemicals database"""
    seed_mongo_documents('chemicals', filename)
    db_client[gc.CHEMICALS['database']][gc.CHEMICALS['collection']].create_index([('smiles', pymongo.HASHED)])


def seed_reactions(filename):
    """Use the specified file to seed the reactions database"""
    seed_mongo_documents('reactions', filename)


def seed_retro_templates(filename):
    """Use the specified file to seed the retrosynthetic template database"""
    seed_mongo_documents('retro_templates', filename)
    db_client[gc.RETRO_TEMPLATES['database']][gc.RETRO_TEMPLATES['collection']].create_index([('index', pymongo.ASCENDING)])


def seed_forward_templates(filename):
    """Use the specified file to seed the forward template database"""
    seed_mongo_documents('forward_templates', filename)


def seed_mongo_documents(collection_name, file_path, db_name='askcos', clear=True, bypass_document_validation=True):
    """
    Seed the specified data file into a mongo collection.

    Args:
        collection_name (str): name of the collection to seed into
        file_path (str): path to the data file to seed
        db_name (str, optional): name of the mongo database to seed into
        clear (bool, optional): clear existing data in the specified collection
        bypass_document_validation (bool, optional): skip document validation by mongo db
    """
    coll = db_client[db_name][collection_name]
    if clear:
        res = coll.delete_many({})
    
    with gzip.open(file_path, 'rb') as f:
        docs = json.loads(f.read().decode('utf-8'))
    
    res = coll.insert_many(docs, ordered=False, bypass_document_validation=bypass_document_validation)
    print('Inserted {} documents into {}.{}'.format(len(res.inserted_ids), db_name, collection_name))
