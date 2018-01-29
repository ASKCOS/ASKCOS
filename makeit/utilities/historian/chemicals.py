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
import os
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
        self.occurrences_flat = defaultdict(lambda: [0, 0, [], []])

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

    def dump_to_file(self, file_path=gc.historian_data):
        '''
        Write the data from the online datbases to a local file
        '''
        if not (self.REACTION_DB and self.CHEMICAL_DB):
            MyLogger.print_and_log(
                "No database information to output to file.", historian_loc, level=3)
            return

        with open(file_path, 'wb') as file:
            pickle.dump(dict(self.occurrences), file, gc.protocol)
            pickle.dump(dict(self.occurrences_flat), file, gc.protocol)
        MyLogger.print_and_log(
                "Saved to {}".format(file_path), historian_loc, level=1)

    def load_from_file(self, file_path=gc.historian_data):
        '''
        Load the data for the pricer from a locally stored file instead of from the online database.
        '''
        if os.path.isfile(file_path):
            with open(file_path, 'rb') as file:
                self.occurrences = defaultdict(lambda: [0, 0, [], []], pickle.load(file))
                self.occurrences_flat = defaultdict(lambda: [0, 0, [], []], pickle.load(file))
        else:
            self.load_databases()
            self.load(refs=True)
            self.dump_to_file()

    def load(self, online=True, CHEMICAL_DB=None, REACTION_DB=None, refs=False):
        '''
        Loads the object from a MongoDB collection
        '''
        MyLogger.print_and_log(
            'Loading historian.', historian_loc)

        # Save collection source use online option to load either from local
        # file or from online database.
        self.occurrences = defaultdict(lambda: [0, 0, [], []])
        self.occurrences_flat = defaultdict(lambda: [0, 0, [], []])

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
        xrn_to_smiles = {}
        xrn_to_smiles_flat = {}
        for chem_doc in tqdm(self.CHEMICAL_DB.find({'SMILES': {'$ne': ''}},
                                                   ['_id', 'SMILES'], no_cursor_timeout=True)):
            m = Chem.MolFromSmiles(str(chem_doc['SMILES']))
            if not m:
                continue
            try:
                xrn_to_smiles[chem_doc['_id']] = Chem.MolToSmiles(m, True)
                xrn_to_smiles_flat[chem_doc['_id']] = Chem.MolToSmiles(m, False)
            except Exception as e:
                continue

        MyLogger.print_and_log(
            'Created xrn to smiles dict ({} entries)'.format(len(xrn_to_smiles)), historian_loc)

        # Then get reactions
        for reaction_doc in tqdm(self.REACTION_DB.find({}, ['_id', 'RX_PXRN', 'RX_RXRN', 'RX_NVAR'],
                                                       no_cursor_timeout=True)):
            for reactant_xrn in reaction_doc['RX_RXRN']:
                if reactant_xrn in xrn_to_smiles:
                    self.occurrences[xrn_to_smiles[reactant_xrn]][0] += reaction_doc['RX_NVAR']
                    self.occurrences[xrn_to_smiles[reactant_xrn]][2].append(reaction_doc['_id'])
                if reactant_xrn in xrn_to_smiles_flat:
                    self.occurrences_flat[
                        xrn_to_smiles_flat[reactant_xrn]][0] += reaction_doc['RX_NVAR']
                    self.occurrences_flat[xrn_to_smiles_flat[reactant_xrn]][2].append(reaction_doc['_id'])

            for product_xrn in reaction_doc['RX_PXRN']:
                if product_xrn in xrn_to_smiles:
                    self.occurrences[xrn_to_smiles[product_xrn]][1] += reaction_doc['RX_NVAR']
                    self.occurrences[xrn_to_smiles[product_xrn]][3].append(reaction_doc['_id'])
                if product_xrn in xrn_to_smiles_flat:
                    self.occurrences_flat[
                        xrn_to_smiles_flat[product_xrn]][1] += reaction_doc['RX_NVAR']
                    self.occurrences_flat[xrn_to_smiles_flat[product_xrn]][3].append(reaction_doc['_id'])

        MyLogger.print_and_log('Historian is fully loaded.', historian_loc)

    def lookup_smiles(self, smiles, alreadyCanonical=False, isomericSmiles=True, refs=False):
        '''
        Looks up a price by SMILES. Tries it as-entered and then 
        re-canonicalizes it in RDKit unl ess the user specifies that
        the string is definitely already canonical.
        '''
        info = self.occurrences_flat[smiles]
        if info != [0, 0, [], []]:
            return tup_to_dict(info, refs=refs)

        info = self.occurrences[smiles]
        if info != [0, 0, [], []]:
            return tup_to_dict(info, refs=refs)

        if not alreadyCanonical:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0.
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)

            info = self.occurrences_flat[smiles]
            if info != [0, 0, [], []]:
                return tup_to_dict(info, refs=refs)

            info = self.occurrences[smiles]
            if info != [0, 0, [], []]:
                return tup_to_dict(info, refs=refs)

        return tup_to_dict(info, refs=refs)


if __name__ == '__main__':
    # Load and dump chemhistorian
    chemhistorian = ChemHistorian()
    chemhistorian.load_from_file()

    print(chemhistorian.lookup_smiles('CCCCO'))