import makeit.global_config as gc
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(4)
import rdkit.Chem as Chem
from collections import defaultdict
from tqdm import tqdm
from makeit.utilities.io.logging import MyLogger
import cPickle as pickle
from pymongo import MongoClient
from multiprocessing import Manager
import time
import os
import hashlib
historian_loc = 'chemhistorian'


def tup_to_dict(info=[0, 0, [], []], refs=False):
    '''Takes a historical entry as a tuple and returns a labeled dict'''
    return {
        'as_reactant': info[0],
        'as_product': info[1],
        'as_reactant_refs': info[2] if refs else [],
        'as_product_refs': info[3] if refs else [],
    }


class ChemHistorian:
    '''
    The Historian class is used to look up chemicals to see how often they have
    been seen in Reaxys
    '''

    def __init__(self, by_xrn=False, REACTIONS=None, CHEMICALS=None):

        self.CHEMICAL_DB = CHEMICALS
        self.REACTION_DB = REACTIONS

        self.occurrences = defaultdict(lambda: [0, 0, [], []])
        self._loaded = False
        self._compressed = False

    def load_databases(self):
        '''
        Load the pricing data from the online database
        '''
        db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                'id'], connect=gc.MONGO['connect'])
        db = db_client[gc.REACTIONS['database']]
        self.REACTION_DB = db[gc.REACTIONS['collection']]
        db = db_client[gc.CHEMICALS['database']]
        self.CHEMICAL_DB = db[gc.CHEMICALS['collection']]

    def dump_to_file(self, file_path=gc.historian_data, refs=False, compressed=False):
        '''
        Write the data from the online datbases to a local file
        '''

        if not refs:
            file_path += '_no_refs'
            for k in self.occurrences.keys():
                self.occurrences[k] = tuple(self.occurrences[k][0:2])

        if compressed:
            file_path += '_compressed'

        with open(file_path, 'wb') as file:
            pickle.dump(dict(self.occurrences), file, gc.protocol)
        MyLogger.print_and_log(
                "Saved to {}".format(file_path), historian_loc, level=1)

    def load_from_file(self, file_path=gc.historian_data, refs=False, compressed=False):
        '''
        Load the data for the pricer from a locally stored file instead of from the online database.
        '''

        if not refs:
            file_path += '_no_refs'
        if compressed:
            file_path += '_compressed'

        if os.path.isfile(file_path):
            with open(file_path, 'rb') as file:
                self.occurrences = pickle.load(file)
                self._loaded = True
                if compressed:
                    self._compressed = True
        else:
            raise ValueError('File does not exist!')

    def upload_to_db(self):
        db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                'id'], connect=gc.MONGO['connect'])
        self.CHEMICAL_HISTORY_DB = db_client[gc.CHEMICAL_HISTORY['database']][gc.CHEMICAL_HISTORY['collection']]

        docs = []
        failed_docs = []
        ctr = 0
        for (i, (smi, info)) in tqdm(enumerate(self.occurrences.iteritems())):
            doc = tup_to_dict(info, refs=True)
            doc['smiles'] = smi 
            docs.append(doc)
            ctr += 1
            if ctr >= 500:
                try:
                    self.CHEMICAL_HISTORY_DB.insert_many(docs, ordered=False)
                except KeyboardInterrupt:
                    raise KeyboardInterrupt
                except Exception as e:
                    print(e)
                    time.sleep(5)
                    for doc in docs:
                        try:
                            self.CHEMICAL_HISTORY_DB.insert_one(doc)
                        except Exception as e:
                            print('#### {}'.format(e))
                            print(doc)
                            failed_docs.append(doc)
                docs = []
                ctr = 0
        try:
            self.CHEMICAL_HISTORY_DB.insert_many(docs)
        except Exception as e:
                    print(e)
                    for doc in docs:
                        try:
                            self.CHEMICAL_HISTORY_DB.insert_one(doc)
                        except Exception as e:
                            print('## {}'.format(e))
                            failed_docs.append(doc)

        with open('failed_docs.pickle', 'wb') as fid:
            pickle.dump(failed_docs, fid, -1)
        # Note: going back, failed insertions were due to duplicatekeyerrors
        # (i.e., should not have been an error) or one SMILES "Cl" had too many
        # references and exceeded the maximum document size. Manually truncating
        # the reference list solved that problem

    def load(self, online=True, CHEMICAL_DB=None, REACTION_DB=None, refs=True):
        '''
        Loads the object from a MongoDB collection
        '''
        MyLogger.print_and_log(
            'Loading historian.', historian_loc)

        # Save collection source use online option to load either from local
        # file or from online database.
        self.occurrences = defaultdict(lambda: [0, 0, [], []])

        if self.REACTION_DB and self.CHEMICAL_DB:
            pass
        elif not (CHEMICAL_DB and REACTION_DB):
            if online:
                self.load_databases()
            else:
                self.load_from_file()
        else:
            self.CHEMICAL_DB = CHEMICAL_DB
            self.REACTION_DB = REACTION_DB

        # First get xrn to smiles dict
        xrn_to_smiles = defaultdict(lambda: [])
        for chem_doc in tqdm(self.CHEMICAL_DB.find({'SMILES': {'$ne': ''}},
                                                   ['_id', 'SMILES'], no_cursor_timeout=True)):
            m = Chem.MolFromSmiles(str(chem_doc['SMILES']))
            if not m:
                continue
            try:
                xrn_to_smiles[chem_doc['_id']] = set(Chem.MolToSmiles(m, True).split('.'))
            except Exception as e:
                continue

        MyLogger.print_and_log(
            'Created xrn to smiles dict ({} entries)'.format(len(xrn_to_smiles)), historian_loc)

        # Then get reactions
        for reaction_doc in tqdm(self.REACTION_DB.find({}, ['_id', 'RX_PXRN', 'RX_RXRN', 'RX_NVAR'],
                                                       no_cursor_timeout=True)):
            for reactant_xrn in reaction_doc['RX_RXRN']:
                for smi in xrn_to_smiles[reactant_xrn]:
                    self.occurrences[smi][0] += reaction_doc['RX_NVAR']
                    if refs:
                        self.occurrences[smi][2].append(reaction_doc['_id'])

            for product_xrn in reaction_doc['RX_PXRN']:
                for smi in xrn_to_smiles[product_xrn]:
                    self.occurrences[smi][1] += reaction_doc['RX_NVAR']
                    if refs:
                        self.occurrences[smi][3].append(reaction_doc['_id'])

        self._loaded = True
        MyLogger.print_and_log('Historian is fully loaded.', historian_loc)

    def lookup_smiles(self, smiles, alreadyCanonical=False, isomericSmiles=True, refs=False):
        '''
        Looks up a price by SMILES. Tries it as-entered and then 
        re-canonicalizes it in RDKit unless the user specifies that
        the string is definitely already canonical.
        '''

        if not isomericSmiles:
            raise ValueError('Not intended to be used for non-isomeric!')

        if not alreadyCanonical:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0.
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)

        try:
            if self._compressed:
                info = self.occurrences[int(hashlib.md5(smiles).hexdigest(), 16)]
            else:
                info = self.occurrences[smiles]
        except KeyError:
            info = [0, 0, [], []]

        return tup_to_dict(info, refs=refs)

    def compress_keys(self):
        '''Convert keys to hashed values to save space'''
        new_occurrences = {}
        for k in self.occurrences.keys():
            k_compressed = int(hashlib.md5(k).hexdigest(), 16)
            new_occurrences[k_compressed] = self.occurrences[k]
        del self.occurrences
        self.occurrences = new_occurrences
        self._compressed = True


if __name__ == '__main__':
    import time
    time.sleep(1)

    # Load and dump chemhistorian
    chemhistorian = ChemHistorian()
    # print('loading refs')
    # chemhistorian.load_from_file(refs=True)
    # # print('uploading to db')
    # # chemhistorian.upload_to_db()
    # print('dumping counts only')
    # chemhistorian.dump_to_file(refs=False)


    #### Compress keys
    # print('loading no refs')
    # chemhistorian.load_from_file(refs=False)
    # print('compressing')
    # chemhistorian.compress_keys()
    # print('dumping compressed')
    # chemhistorian.dump_to_file(refs=False, compressed=True)


    print('loading no refs, compressed')
    chemhistorian.load_from_file(refs=False, compressed=True)
    print(chemhistorian.lookup_smiles('CCCCO'))

    time.sleep(10)