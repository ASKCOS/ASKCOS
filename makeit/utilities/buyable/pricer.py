import makeit.global_config as gc
import rdkit.Chem as Chem
from collections import defaultdict
from tqdm import tqdm
from makeit.utilities.io.logger import MyLogger
import makeit.utilities.io.pickle as pickle
from pymongo import MongoClient
from multiprocessing import Manager
import os
pricer_loc = 'pricer'


class Pricer:
    '''
    The Pricer class is used to look up the ppg of chemicals if they
    are buyable.
    '''

    def __init__(self, by_xrn=False, done=None, BUYABLES=None, CHEMICALS=None):

        self.CHEMICAL_DB = CHEMICALS
        self.BUYABLE_DB = BUYABLES
        self.by_xrn = False
        self.done = done
        self.prices = defaultdict(float)  # default 0 ppg means not buyable
        # default 0 ppg means not buyable
        self.prices_flat = defaultdict(float)
        self.prices_by_xrn = defaultdict(float)

    def load_databases(self):
        '''
        Load the pricing data from the online database
        '''
        db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                'id'], connect=gc.MONGO['connect'])
        db = db_client[gc.BUYABLES['database']]
        self.BUYABLE_DB = db[gc.BUYABLES['collection']]
        db = db_client[gc.CHEMICALS['database']]
        self.CHEMICAL_DB = db[gc.CHEMICALS['collection']]

    def dump_to_file(self, file_path):
        '''
        Write the data from the online datbases to a local file
        '''
        with open(file_path, 'wb') as file:
            pickle.dump(self.prices, file)
            pickle.dump(self.prices_flat, file)
            pickle.dump(self.prices_by_xrn, file)

    def load(self):
        '''
        Load the data for the pricer from a locally stored file instead of from the online database.
        '''
        from makeit.utilities.io.files import get_pricer_path
        file_path = get_pricer_path(
            gc.CHEMICALS['database'],
            gc.CHEMICALS['collection'],
            gc.BUYABLES['database'],
            gc.BUYABLES['collection'],
        )
        if os.path.isfile(file_path):
            with open(file_path, 'rb') as file:
                self.prices = defaultdict(float, pickle.load(file))
                self.prices_flat = defaultdict(float, pickle.load(file))
                self.prices_by_xrn = defaultdict(float, pickle.load(file))
        else:
            self.load_databases()
            # self.load_from_database()
            # self.dump_to_file(file_path)

    def load_from_database(self, max_ppg=1e10):
        '''
        Loads the object from a MongoDB collection containing transform
        template records.
        '''
        MyLogger.print_and_log('Loading pricer with buyable limit of ${} per gram.'.format(max_ppg), pricer_loc)

        # Save collection source use online option to load either from local
        # file or from online database.
        self.prices = defaultdict(float)  # default 0 ppg means not buyable
        # default 0 ppg means not buyable
        self.prices_flat = defaultdict(float)
        self.prices_by_xrn = defaultdict(float)

        buyable_dict = {}

        # First pull buyables source (smaller)
        for buyable_doc in self.BUYABLE_DB.find({'source': {'$ne': 'LN'}},
                        ['ppg', 'smiles', 'smiles_flat'],
                        no_cursor_timeout=True):

            if buyable_doc['ppg'] > max_ppg:
                continue

            # Deal with multi-species ones (pretty common):
            for smiles in buyable_doc['smiles'].split('.'):

                if self.prices[smiles]:  # already in dict as non-zero, so keep cheaper
                    self.prices[smiles] = min(
                        buyable_doc['ppg'], self.prices[smiles])
                else:
                    self.prices[smiles] = buyable_doc['ppg']

            for smiles_flat in buyable_doc['smiles_flat'].split('.'):
                if self.prices_flat[smiles_flat]:
                    self.prices_flat[smiles_flat] = min(
                        buyable_doc['ppg'], self.prices_flat[smiles_flat])
                else:
                    self.prices_flat[smiles_flat] = buyable_doc['ppg']

        if self.by_xrn:
            # Then pull chemicals source for XRNs (larger)
            for chemical_doc in tqdm(self.CHEMICAL_DB.find({'buyable_id': {'$gt': -1}},
                    ['buyable_id'], no_cursor_timeout=True)):
                if 'buyable_id' not in chemical_doc:
                    continue
                self.prices_by_xrn[chemical_doc['_id']] = buyable_dict[
                    chemical_doc['buyable_id']]

        MyLogger.print_and_log('Pricer has been loaded.', pricer_loc)

        # multiprocessing notify done
        if self.done == None:
            pass
        else:
            self.done.value = 1

    def lookup_smiles(self, smiles, alreadyCanonical=False, isomericSmiles=True):
        if not alreadyCanonical:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0.
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)

        entry = self.BUYABLE_DB.find_one({'smiles': smiles})

        if entry:
            return entry['ppg']
        else:
            return 0.

    def lookup_smiles_old(self, smiles, alreadyCanonical=False, isomericSmiles=True):
        '''
        Looks up a price by SMILES. Tries it as-entered and then
        re-canonicalizes it in RDKit unl ess the user specifies that
        the string is definitely already canonical.
        '''
        ppg = self.prices_flat[smiles]
        if ppg:
            return ppg

        ppg = self.prices[smiles]
        if ppg:
            return ppg

        if not alreadyCanonical:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0.
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)

            ppg = self.prices_flat[smiles]
            if ppg:
                return ppg

            ppg = self.prices[smiles]
            if ppg:
                return ppg

        return ppg

    def lookup_xrn(self, xrn):
        '''
        Looks up a price by Reaxys XRN.
        '''
        if not self.by_rxn:
            raise ValueError('Not initialized to look up prices by XRN!')
        return self.prices_by_xrn[xrn]


if __name__ == '__main__':
    pricer = Pricer()
    pricer.load()
    print(pricer.lookup_smiles('CCCCCO'))
    print(pricer.lookup_smiles('CCCCXCCO'))
