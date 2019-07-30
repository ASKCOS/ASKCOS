from makeit.prioritization.prioritizer import Prioritizer
import makeit.global_config as gc
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.buyable.pricer import Pricer
from makeit.utilities.io.logger import MyLogger
import math
import sys
import os
import time
import os
import makeit.utilities.io.pickle as pickle
from numpy import inf
scscore_prioritizer_loc = 'scscoreprioritizer'

class SCScorePrecursorPrioritizer(Prioritizer):
    """Standalone, importable SCScorecorer model.

    It does not have tensorflow as a dependency and is a more attractive option
    for deployment. The calculations are fast enough that there is no real
    reason to use GPUs (via tf) instead of CPUs (via np).

    Attributes:
        vars (list of np.ndarry of np.ndarray of np.float32): Weights and bias
            of given model.
        FP_rad (int): Fingerprint radius.
        FP_len (int): Fingerprint length.
        score_scale (float): Upper-bound of scale for scoring.
        pricer (Pricer or None): Pricer instance to lookup chemical costs.
    """

    def __init__(self, score_scale=5.0):
        """Initializes SCScorePrecursorPrioritizer.

        Args:
            score_scale (float, optional): Upper-bound of scale for scoring.
                (default: {5.0})
        """
        self.vars = []
        self.FP_rad = 2
        self.score_scale = score_scale
        self._restored = False
        self.pricer = None
        self._loaded = False

    def load_model(self, FP_len=1024, model_tag='1024bool'):
        """Loads model from given tag.

        Args:
            FP_len (int, optional): Fingerprint length. (default: {1024})
            model_tag (str, optional): Tag of model to load.
                (default: {'1024bool'})
        """
        self.FP_len = FP_len
        if model_tag != '1024bool' and model_tag != '1024uint8' and model_tag != '2048bool':
            MyLogger.print_and_log(
                'Non-existent SCScore model requested: {}. Using "1024bool" model'.format(model_tag), scscore_prioritizer_loc, level=2)
            model_tag = '1024bool'
        filename = 'trained_model_path_'+model_tag
        with open(gc.SCScore_Prioritiaztion[filename], 'rb') as fid:
            self.vars = pickle.load(fid)
        if gc.DEBUG:
            MyLogger.print_and_log('Loaded synthetic complexity score prioritization model from {}'.format(
            gc.SCScore_Prioritiaztion[filename]), scscore_prioritizer_loc)

        if 'uint8' in gc.SCScore_Prioritiaztion[filename]:
            def mol_to_fp(mol):
                """Returns fingerprint of molecule for uint8 model.

                Args:
                    mol (Chem.rdchem.Mol or None): Molecule to get fingerprint
                        of.

                Returns:
                    np.ndarray of np.uint8: Fingerprint of given molecule.
                """
                if mol is None:
                    return np.array((self.FP_len,), dtype=np.uint8)
                fp = AllChem.GetMorganFingerprint(
                    mol, self.FP_rad, useChirality=True)  # uitnsparsevect
                fp_folded = np.zeros((self.FP_len,), dtype=np.uint8)
                for k, v in fp.GetNonzeroElements().items():
                    fp_folded[k % self.FP_len] += v
                return np.array(fp_folded)
        else:
            def mol_to_fp(mol):
                """Returns fingerprint of molecule for bool model.

                Args:
                    mol (Chem.rdchem.Mol or None): Molecule to get fingerprint
                        of.

                Returns:
                    np.ndarray of np.bool or np.float32: Fingerprint of given
                        molecule.
                """
                if mol is None:
                    return np.zeros((self.FP_len,), dtype=np.float32)
                return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, self.FP_rad, nBits=self.FP_len,
                                                                      useChirality=True), dtype=np.bool)
        self.mol_to_fp = mol_to_fp

        self.pricer = Pricer()
        self.pricer.load()
        self._restored = True
        self._loaded = True

    def smi_to_fp(self, smi):
        """Returns fingerprint of molecule from given SMILES string.

        Args:
            smi (str): SMILES string of given molecule.
        """
        if not smi:
            return np.zeros((self.FP_len,), dtype=np.float32)
        return self.mol_to_fp(Chem.MolFromSmiles(str(smi)))

    def apply(self, x):
        """Applies model to a fingerprint to calculate score.

        Args:
            x (np.ndarray): Fingerprint of molecule to apply model to.

        Returns:
            float: Score of molecule.
        """
        if not self._restored:
            raise ValueError('Must restore model weights!')
        # Each pair of vars is a weight and bias term
        for i in range(0, len(self.vars), 2):
            last_layer = (i == (len(self.vars)-2))
            W = self.vars[i]
            b = self.vars[i+1]
            x = np.dot(W.T, x) + b
            if not last_layer:
                x = x * (x > 0)  # ReLU
        x = 1 + (self.score_scale - 1) * sigmoid(x)
        return x

    def get_priority(self, retroProduct, **kwargs):
        """Returns priority of given product.

        Args:
            retroProduct (str or RetroPrecursor): Product to calculate score
                for.
            **kwargs: Additional optional arguments. Used for mode.

        Returns:
            float: Priority of given product.
        """
        mode = kwargs.get('mode', gc.max)
        if not self._loaded:
            self.load_model()

        if not isinstance(retroProduct, str):
            scores = []
            for smiles in retroProduct.smiles_list:
                scores.append(self.get_score_from_smiles(smiles))
            return -self.merge_scores(scores, mode=mode)
        else:
            return -self.get_score_from_smiles(retroProduct)
        if not retroProduct:
            return -inf

    def merge_scores(self, list_of_scores, mode=gc.max):
        """Merges list of scores into a single score based on a given mode.

        Args:
            list_of_scores (list of floats): Scores to be merged.
            mode (str, optional): Function to merge by. (default: {gc.max})

        Returns:
            float: Merged scores.
        """
        if mode == gc.mean:
            return np.mean(list_of_scores)
        elif mode == gc.geometric:
            return np.power(np.prod(list_of_scores), 1.0/len(list_of_scores))
        elif mode == gc.pow8:
            pow8 = []
            for score in list_of_scores:
                pow8.append(8**score)
            return np.sum(pow8)
        else:
            return np.max(list_of_scores)

    def get_score_from_smiles(self, smiles, noprice=False):
        """Returns score of molecule from given SMILES string.

        Args:
            smiles (str): SMILES string of molecule.
            noprice (bool, optional): Whether to not use the molecules price as
                its score, if available. (default: {False})
        """
        # Check buyable
        if not noprice:
            ppg = self.pricer.lookup_smiles(smiles, alreadyCanonical=True)
            if ppg:
                return ppg / 100.

        fp = np.array((self.smi_to_fp(smiles)), dtype=np.float32)
        if sum(fp) == 0:
            cur_score = 0.
        else:
            # Run
            cur_score = self.apply(fp)
        return cur_score


def sigmoid(x):
    """Returns sigmoid of x.

    Args:
        x (float): Input value.
    """
    if x < -10:
        return 0
    if x > 10:
        return 1
    return 1. / (1 + math.exp(-x))

if __name__ == '__main__':
    model = SCScorePrecursorPrioritizer()
    model.load_model(model_tag='1024bool')
    smis = ['CCCOCCC', 'CCCC']
    for smi in smis:
        sco = model.get_score_from_smiles(smi, noprice=True)
        print('{} <--- {}'.format(sco, smi))

    # model = SCScorer()
    # model.load_model(model_tag='2048bool', FP_len=2048)
    # smis = ['CCCOCCC', 'CCCNc1ccccc1']
    # for smi in smis:
    #     sco = model.get_priority(smi)
    #     print('{} <--- {}'.format(sco, smi))

    # model = SCScorer()
    # model.load_model(model_tag='1024uint8')
    # smis = ['CCCOCCC', 'CCCNc1ccccc1']
    # for smi in smis:
    #     sco = model.get_priority(smi)
    #     print('{} <--- {}'.format(sco, smi))
