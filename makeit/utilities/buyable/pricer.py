import gzip
import json
import makeit.global_config as gc
import rdkit.Chem as Chem
from collections import defaultdict
from tqdm import tqdm
from makeit.utilities.io.logger import MyLogger
import makeit.utilities.io.pickle as pickle
from pymongo import MongoClient, errors
from multiprocessing import Manager
import os
pricer_loc = 'pricer'


class Pricer:
    '''
    The Pricer class is used to look up the ppg of chemicals if they
    are buyable.
    '''

    def __init__(self, use_db=False, BUYABLES_DB=None):

        self.BUYABLES_DB = BUYABLES_DB
        self.use_db = use_db
        self.prices = defaultdict(float)  # default 0 ppg means not buyable

    def load(self, file_name=gc.BUYABLES['file_name']):
        '''
        Load pricer information. Either create connection to MongoDB or load from local file.
        If connection to MongoDB cannot be made, fallback and try to load from local file.
        '''
        if self.use_db:
            self.load_databases()
            return

        if not os.path.isfile(file_name):
            MyLogger.print_and_log('Buyables file does not exist file: {}'.format(file_name), pricer_loc)
            return

        self.load_from_file(file_name)

    def load_databases(self):
        '''
        Load the pricing data from the online database
        '''
        db_client = MongoClient(
            gc.MONGO['path'], 
            gc.MONGO['id'],
            connect=gc.MONGO['connect'], 
            serverSelectionTimeoutMS=1000
        )

        try:
            db_client.server_info()
        except errors.ServerSelectionTimeoutError:
            MyLogger.print_and_log('Cannot connect to mongodb to load prices', pricer_loc)
            self.use_db = False
            self.load()
            return
        
        db = db_client[gc.BUYABLES['database']]
        self.BUYABLES_DB = db[gc.BUYABLES['collection']]

    def dump_to_file(self, file_path):
        '''
        Write prices to a local file
        '''
        prices = []
        for k, v in self.prices.items():
            tmp = v.copy()
            tmp['smiles'] = k
            prices.append(tmp)

        with gzip.open(file_path, 'wb') as f:
            json.dump(prices, f)

    def load_from_file(self, file_name):
        '''
        Load buyables information from local file
        '''
        with gzip.open(file_name, 'rb') as f:
            prices = json.loads(f.read().decode('utf-8'))

        for p in prices:
            smiles = p.pop('smiles', '')
            if smiles:
                self.prices[smiles] = p.pop('ppg')
        MyLogger.print_and_log('Loaded prices from flat file', pricer_loc)

    def lookup_smiles(self, smiles, alreadyCanonical=False, isomericSmiles=True):
        '''
        Looks up a price by SMILES. Canonicalize smiles string unless 
        the user specifies that the smiles string is definitely already 
        canonical. If the DB connection does not exist, look up from 
        prices dictionary attribute, otherwise lookup from DB.
        If multiple entries exist in the DB, return the lowest price.
        '''
        if not alreadyCanonical:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0.
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)

        if self.use_db:
            cursor = self.BUYABLES_DB.find({
                'smiles': smiles,
                'source': {'$ne': 'LN'}
            })
            if cursor.count():
                return min([doc['ppg'] for doc in cursor])
            else:
                return 0.
        else:
            return self.prices[smiles]


if __name__ == '__main__':
    pricer = Pricer()
    pricer.load()
    print(pricer.lookup_smiles('CCCCCO'))
    print(pricer.lookup_smiles('CCCCXCCO'))
