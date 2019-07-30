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
import os
historian_loc = 'reactionhistorian'


def tup_to_dict(info=[0, []], refs=False):
    """Returns a labeled dict from a given historical entry as a tuple.

    Args:
        info ((int, list), optional): Entry to convert into dict.
            (default: {[0, []]})
        refs (bool, optional): Whether to include references or only counts.
            (default: {False})
    """
    return {
        'count': info[0],
        'refs': info[1] if refs else [],
    }


class ReactionHistorian:
    """Looks up chemicals to see how often they have been seen in Reaxys.

    Attributes:
        REACTION_DB (MongoDB): Database of reactions.
        CHEMICAL_HISTORY_DB (MongoDB): Database of chemical history.
        occurrences (defaultdict):
        occurrences_flat (defaultdict):
    """

    def __init__(self, by_xrn=False, REACTIONS=None):
        """Initializes ReactionHistorian.

        Args:
            by_xrn (bool, optional): ?? (default: {False})
            REACTIONS (None or MongoDB, optional): ?? (default: {None})
        """
        self.REACTION_DB = REACTIONS

        self.occurrences = defaultdict(lambda: [0, []])
        self.occurrences_flat = defaultdict(lambda: [0, []])

    def load_databases(self):
        """Loads the history data from the online database."""
        db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                'id'], connect=gc.MONGO['connect'])
        db = db_client[gc.REACTIONS['database']]
        self.REACTION_DB = db[gc.REACTIONS['collection']]

    def dump_to_file(self, file_path=gc.reactionhistorian_data):
        """Writes the data from the online datbases to a local file.

        Args:
            file_path (str, optional): Path to the output file.
                (default: {gc.reactionhistorian_data})
        """
        if not self.REACTION_DB:
            MyLogger.print_and_log(
                "No database information to output to file.", historian_loc, level=3)
            return

        with open(file_path, 'wb') as file:
            pickle.dump(dict(self.occurrences), file)
            pickle.dump(dict(self.occurrences_flat), file)

        MyLogger.print_and_log(
                "Saved to {}".format(file_path), historian_loc, level=1)

    def load_from_file(self, file_path=gc.reactionhistorian_data, testing=False):
        """Loads the data for the pricer from a locally stored file.

        Args:
            file_path (str, optional): Path to the input file.
                (default: {gc.reactionhistorian_data})
            testing (bool, optional): Whether to only run a test.
                (default: {False})
        """

        if testing:
            self.occurrences = defaultdict(lambda: [0, []], {
                    'CCO>>CCBr': [2, ['rxn1', 'rxn2']],
                    'CCCCC>>CCC=CC': [1, ['rxn3']],
                })
            self.occurrences_flat = defaultdict(lambda: [0, []], {
                    'CCO>>CCBr': [2, ['rxn1', 'rxn2']],
                    'CCCCC>>CCC=CC': [1, ['rxn3']],
                })
            return


        if os.path.isfile(file_path):
            with open(file_path, 'rb') as file:
                self.occurrences = defaultdict(lambda: [0, []], pickle.load(file))
                self.occurrences_flat = defaultdict(lambda: [0, []], pickle.load(file))
        else:
            self.load_databases()
            self.load()
            self.dump_to_file()

    def load(self, online=True, REACTION_DB=None):
        """Loads the object from a MongoDB collection.

        Args:
            online (bool, optional): Whether to load from the online database or
                a local file. (default: {True})
            REACTION_DB (None or ??, optional): ??
        """
        MyLogger.print_and_log(
            'Loading reaction historian.', historian_loc)

        # Save collection source use online option to load either from local
        # file or from online database.
        self.occurrences = defaultdict(lambda: [0, []])
        self.occurrences_flat = defaultdict(lambda: [0, []])

        if self.REACTION_DB:
            pass
        elif not REACTION_DB:
            if online:
                self.load_databases()
            else:
                self.load_from_file()
        else:
            self.REACTION_DB = REACTION_DB

        # Then get reactions
        for reaction_doc in tqdm(self.REACTION_DB.find({'RXN_SMILES': {'$exists': True}},
                    ['_id', 'RXN_SMILES', 'RX_NVAR'], no_cursor_timeout=True)):

            # Parse and canonicalize
            reactants = Chem.MolFromSmiles(str(reaction_doc['RXN_SMILES'].split('>')[0]))
            products = Chem.MolFromSmiles(str(reaction_doc['RXN_SMILES'].split('>')[2]))

            if not reactants or not products:
                continue

            [a.ClearProp('molAtomMapNumber') for a in reactants.GetAtoms()]
            [a.ClearProp('molAtomMapNumber') for a in products.GetAtoms()]

            try:
                reaction_smiles = '{}>>{}'.format(
                    '.'.join(sorted(Chem.MolToSmiles(reactants, True).split('.'))),
                    '.'.join(sorted(Chem.MolToSmiles(products, True).split('.'))),
                )
                reaction_smiles_flat = '{}>>{}'.format(
                    '.'.join(sorted(Chem.MolToSmiles(reactants, False).split('.'))),
                    '.'.join(sorted(Chem.MolToSmiles(products, False).split('.'))),
                )
            except Exception as e:
                continue

            self.occurrences[reaction_smiles][0] += reaction_doc['RX_NVAR']
            self.occurrences[reaction_smiles][1].extend(
                ['{}-{}'.format(reaction_doc['_id'], i+1) for i in range(reaction_doc['RX_NVAR'])]
            )

            self.occurrences_flat[reaction_smiles_flat][0] += reaction_doc['RX_NVAR']
            self.occurrences_flat[reaction_smiles_flat][1].extend(
                ['{}-{}'.format(reaction_doc['_id'], i+1) for i in range(reaction_doc['RX_NVAR'])]
            )

        MyLogger.print_and_log('Reaction historian is fully loaded.', historian_loc)


    def lookup_smiles(self, smiles, alreadyCanonical=False, isomericSmiles=True, refs=False):
        """Looks up number of occurances by SMILES.

        Tries it as-entered and then re-canonicalizes it in RDKit unless the
        user specifies that the string is definitely already canonical.

        Args:
            smiles (str): SMILES string of molecule to look up.
            alreadyCanonical (bool, optional): Whether to SMILES string has
                already been canonicalized. (default: {False})
            isomericSmiles (bool, optional): Whether the SMILES string is
                isomeric. (default: {True})
            refs (bool, optional): Whether to include references.
                (default: {False})

        Returns:
            dict: Labeled dictionaty from tup_to_dict().
        """
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
    reactionhistorian = ReactionHistorian()
    reactionhistorian.load_from_file()

    print(reactionhistorian.lookup_smiles('CCCCO>>CCCCBr'))
