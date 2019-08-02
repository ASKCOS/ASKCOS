import makeit.global_config as gc
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(4)
import rdkit.Chem as Chem
from collections import defaultdict
from tqdm import tqdm
from makeit.utilities.io.logger import MyLogger
import makeit.utilities.io.pickle as pickle
from pymongo import MongoClient
from multiprocessing import Manager
import time
import os
import hashlib
historian_loc = 'chemhistorian'


def tup_to_dict(info=[0, 0, [], []], refs=False):
    """Returns a labeled dict from a given historical entry as a tuple.

    Args:
        info ((int, int, list, list), optional): Entry to convert into dict.
            (default: {[0, 0, [], []]})
        refs (bool, optional): Whether to include references or only counts.
            (default: {False})
    """
    return {
        'as_reactant': info[0],
        'as_product': info[1],
        'as_reactant_refs': info[2] if refs else [],
        'as_product_refs': info[3] if refs else [],
    }


class ChemHistorian:
    """Looks up chemicals to see how often they have been seen in Reaxys.

    Attributes:
        CHEMICAL_DB (MongoDB): Database of chemicals.
        REACTION_DB (MongoDB): Database of reactions.
        CHEMICAL_HISTORY_DB (MongoDB): Database of chemical history.
        occurrences (defaultdict):
    """

    def __init__(self, by_xrn=False, REACTIONS=None, CHEMICALS=None):
        """Initializes ChemHistorian.

        Args:
            by_xrn (bool, optional): ?? (default: {False})
            REACTIONS (None or MongoDB): ?? (default: {None})
            CHEMICALS (None or MongoDB): ?? (default: {None})
        """
        self.CHEMICAL_DB = CHEMICALS
        self.REACTION_DB = REACTIONS

        self.occurrences = defaultdict(lambda: [0, 0, [], []])
        self._loaded = False
        self._compressed = False

    def load_databases(self):
        """Loads the history data from the online database."""
        db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                'id'], connect=gc.MONGO['connect'])
        db = db_client[gc.REACTIONS['database']]
        self.REACTION_DB = db[gc.REACTIONS['collection']]
        db = db_client[gc.CHEMICALS['database']]
        self.CHEMICAL_DB = db[gc.CHEMICALS['collection']]

    def dump_to_file(self, file_path=gc.historian_data, refs=False, compressed=False):
        """Writes the data from the online datbases to a local file.

        Args:
            file_path (str, optional): Path to the output file.
                (default: {gc.historian_data})
            refs (bool, optional): Whether to include the references or just
                the counts. (default: {False})
            compressed (bool, optional): Whether the data is compressed.
                (default: {False})
        """

        if not refs:
            file_path += '_no_refs'
            for k in self.occurrences.keys():
                self.occurrences[k] = tuple(self.occurrences[k][0:2])

        if compressed:
            file_path += '_compressed'

        with open(file_path, 'wb') as file:
            pickle.dump(dict(self.occurrences), file)
        MyLogger.print_and_log(
                "Saved to {}".format(file_path), historian_loc, level=1)

    def load_from_file(self, file_path=gc.historian_data, refs=False, compressed=False):
        """Loads the data for the pricer from a locally stored file.

        Args:
            file_path (str, optional): Path to the input file.
                (default: {gc.historian_data})
            refs (bool, optional): Whether to include the references or just
                the counts. (default: {False})
            compressed (bool, optional): Whether the data is compressed.
                (default: {False})

        Raises:
            ValueError: If file does not exist.
        """

        MyLogger.print_and_log('Loading chemhistorian from file...', historian_loc)

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
        """Uploads the data to the online database."""
        db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                'id'], connect=gc.MONGO['connect'])
        self.CHEMICAL_HISTORY_DB = db_client[gc.CHEMICAL_HISTORY['database']][gc.CHEMICAL_HISTORY['collection']]

        docs = []
        failed_docs = []
        ctr = 0
        for (i, (smi, info)) in tqdm(enumerate(self.occurrences.items())):
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
        """Loads the object from a MongoDB collection.

        Args:
            online (bool, optional): Whether to load from the online database or
                a local file. (default: {True})
            CHEMICAL_DB (None or ??, optional): ??
            REACTION_DB (None or ??, optional): ??
            refs (bool, optional): Whether to include references.
                (default: {True})
        """
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
        """Looks up number of occurances by SMILES.

        Tries it as-entered and then re-canonicalizes it in RDKit unless the
        user specifies that the string is definitely already canonical.

        Args:
            smiles (str): SMILES string of molecule to look up.
            alreadyCanonical (bool, optional): Whether to SMILES string has
                already been canonicalized. (default: {False})
            isomericSmiles (bool, optional): Whether the SMILES string is
                isomeric. Must be true, since function is only intended for
                isomeric. (default: {True})
            refs (bool, optional): Whether to include references.
                (default: {False})

        Returns:
            dict: Labeled dictionaty from tup_to_dict().

        Raises:
            ValueError: If isomericSmiles is False.
        """

        if not isomericSmiles:
            raise ValueError('Not intended to be used for non-isomeric!')

        if not alreadyCanonical:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0.
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)

        try:
            if self._compressed:
                info = self.occurrences[int(hashlib.md5(smiles.encode('utf-8')).hexdigest(), 16)]
            else:
                info = self.occurrences[smiles]
        except KeyError:
            info = [0, 0, [], []]

        return tup_to_dict(info, refs=refs)

    def compress_keys(self):
        """Converts keys to hashed values to save space."""
        new_occurrences = {}
        for k in self.occurrences.keys():
            k_compressed = int(hashlib.md5(k.encode('utf-8')).hexdigest(), 16)
            new_occurrences[k_compressed] = self.occurrences[k]
        del self.occurrences
        self.occurrences = new_occurrences
        self._compressed = True


if __name__ == '__main__':

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
