import makeit.global_config as gc
from makeit.prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.io.logger import MyLogger
import math
import sys
import random
import time
import os
import makeit.utilities.io.pickle as pickle
import tensorflow as tf
import math
from bson.objectid import ObjectId
from rdchiral.initialization import rdchiralReaction, rdchiralReactants
from pymongo import MongoClient



relevance_template_prioritizer_loc = 'relevance_template_prioritizer'

def top1(y_true, y_pred, k=1):
    return tf.keras.metrics.sparse_top_k_categorical_accuracy(y_true, y_pred, k=k)
def top10(y_true, y_pred, k=10):
    return tf.keras.metrics.sparse_top_k_categorical_accuracy(y_true, y_pred, k=k)
def top50(y_true, y_pred, k=50):
    return tf.keras.metrics.sparse_top_k_categorical_accuracy(y_true, y_pred, k=k)
def top100(y_true, y_pred, k=100):
    return tf.keras.metrics.sparse_top_k_categorical_accuracy(y_true, y_pred, k=k)

def loss(labels, logits):
    return tf.keras.losses.sparse_categorical_crossentropy(labels, logits, from_logits=True)

def doc_to_template(document, chiral):
    """Returns a template given a document from the database or file.

    Args:
        document (dict): Document of template from database or file.
        chiral (bool): Whether to pay attention to chirality.

    Returns:
        dict: Retrosynthetic template.
    """
    if 'reaction_smarts' not in document:
        return
    reaction_smarts = str(document['reaction_smarts'])
    if not reaction_smarts:
        return

    # different thresholds for chiral and non chiral reactions
    chiral_rxn = False
    for c in reaction_smarts:
        if c in ('@', '/', '\\'):
            chiral_rxn = True
            break

    # Define dictionary
    template = {
        'name':                 document['name'] if 'name' in document else '',
        'reaction_smarts':      reaction_smarts,
        'incompatible_groups':  document['incompatible_groups'] if 'incompatible_groups' in document else [],
        'reference':            document['reference'] if 'reference' in document else '',
        'references':           document['references'] if 'references' in document else [],
        'rxn_example':          document['rxn_example'] if 'rxn_example' in document else '',
        'explicit_H':           document['explicit_H'] if 'explicit_H' in document else False,
        '_id':                  document['_id'] if '_id' in document else -1,
        'product_smiles':       document['product_smiles'] if 'product_smiles' in document else [],
        'necessary_reagent':    document['necessary_reagent'] if 'necessary_reagent' in document else '',
        'efgs':                 document['efgs'] if 'efgs' in document else None,
        'intra_only':           document['intra_only'] if 'intra_only' in document else False,
        'dimer_only':           document['dimer_only'] if 'dimer_only' in document else False,
    }
    template['chiral'] = chiral_rxn

    # Frequency/popularity score
    template['count'] = document.get('count', 1)

    # Define reaction in RDKit and validate
    try:
        # Force reactants and products to be one pseudo-molecule (bookkeeping)
        reaction_smarts_one = '(' + reaction_smarts.replace('>>', ')>>(') + ')'

        if chiral:
            rxn = rdchiralReaction(str(reaction_smarts_one))
            template['rxn'] = rxn
        else:
            rxn = AllChem.ReactionFromSmarts(
                str(reaction_smarts_one))
            if rxn.Validate()[1] == 0:
                template['rxn'] = rxn
            else:
                template['rxn'] = None

    except Exception as e:
        if gc.DEBUG:
            MyLogger.print_and_log('Couldnt load : {}: {}'.format(
                reaction_smarts_one, e), relevance_template_prioritizer_loc, level=1)
        template['rxn'] = None
        template['rxn_f'] = None
    return template

class RelevanceTemplatePrioritizer(Prioritizer):
    """A template Prioritizer based on template relevance.

    Attributes:
        retro (bool):
        FP_rad (int): Fingerprint radius.
        FP_len (int): Fingerprint length.
        vars (list of np.ndarry of np.ndarray of np.float32): Weights and bias
            of model.
        template_count (int): Maximum number of templates to prioritize.
        max_cum_prob (float):
        batch_size ():
        NK (int):
        session (tensorflow.python.client.session.Session):
        input_mol (tensorflow.python.framework.ops.Tensor):
        mol_hiddens (tensorflow.python.framework.ops.Tensor):
        score (tensorflow.python.framework.ops.Tensor):
        topk (tensorflow.python.framework.ops.Tensor):
        coord (tensorflow.python.training.coordinator.Coordinator):
    """

    def __init__(self, retro=True):
        self.retro = retro
        self.FP_len = 2048
        self.FP_rad = 2
        self.vars = []
        self.template_count = 100
        self.max_cum_prob = 1
        self.batch_size = 1
        self.NK = 100

    def load_model(self, model_path=gc.Relevance_Prioritization['trained_model_path']):
        self.model = tf.keras.models.load_model(model_path)

    def get_topk_from_mol(self, mol, k=100):
        fp = self.mol_to_fp(mol).astype(np.float32)
        cur_scores = self.model.predict(fp.reshape(1, -1)).reshape(-1)
        indices = np.argsort(-cur_scores)[:k]
        scores = softmax(cur_scores[indices])
        return scores.tolist(), indices.tolist()

    def mol_to_fp(self, mol):
        """Returns fingerprint of molecule.

        Args:
            mol (Chem.rdchem.Mol or None): Molecule to get fingerprint
                of.

        Returns:
            np.ndarray of np.float32: Fingerprint of given molecule.
        """
        if mol is None:
            return np.zeros((self.FP_len,), dtype=np.float32)
        return np.array(
            AllChem.GetMorganFingerprintAsBitVect(
                mol, self.FP_rad, nBits=self.FP_len, useChirality=True
            ), dtype=np.float32
        )

    def smi_to_fp(self, smi):
        """Returns fingerprint of molecule from given SMILES string.

        Args:
            smi (str): SMILES string of given molecule.
        """
        if not smi:
            return np.zeros((self.FP_len,), dtype=np.float32)
        return self.mol_to_fp(Chem.MolFromSmiles(smi))

    def get_priority(self, input_tuple, **kwargs):
        """Returns list of templates ordered by relevance.

        Args:
            input_tuple (2-tuple of (list of ??, ??)): Templates to get the
                priority of.
            **kwargs: Additional optional parameters. Used for template_count
                and max_cum_prob.
        """
        (templates, target) = input_tuple
        template_count = kwargs.get('template_count', 100)
        max_cum_prob = kwargs.get('max_cum_prob', 0.995)
        # Templates should be sorted by popularity for indices to be correct!
        probs, top_ids = self.get_topk_from_smi(smi=target, k = min(template_count, len(templates)))
        top_templates = []
        cum_score = 0
        mincount = kwargs.get('mincount', 25)
        mincount_chiral = kwargs.get('mincount_chiral', 10)
        chiral = kwargs.get('chiral', True)
        use_db = kwargs.get('use_db', True)
        load_all = kwargs.get('load_all', gc.PRELOAD_TEMPLATES)
        template_cache = kwargs.get('template_cache', None)
        if not load_all and use_db:
            db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                    'id'], connect=gc.MONGO['connect'])

            db_name = gc.RETRO_TRANSFORMS_CHIRAL['database']
            collection = gc.RETRO_TRANSFORMS_CHIRAL['collection']
            TEMPLATE_DB = db_client[db_name][collection]
        for id, prob in zip(top_ids, probs):
            if load_all:
                template = templates[id]
            elif use_db:
                if template_cache is not None:
                    template = template_cache[ObjectId(templates[id][0])]
                else:
                    document = TEMPLATE_DB.find_one({'_id': ObjectId(templates[id][0])})
                    template = doc_to_template(document, chiral)
            else:
                template = doc_to_template(templates[id], chiral)
            template['score'] = prob
            top_templates.append(template)
            cum_score += prob
            #End loop if max cumulative score is exceeded
            if cum_score >= max_cum_prob:
                break
        return top_templates

    def apply(self, x):
        # Each pair of vars is a weight and bias term
        # (only used for numpy)
        for i in range(0, len(self.vars), 2):
            last_layer = (i == len(self.vars)-2)
            W = self.vars[i]
            b = self.vars[i+1]
            x = np.matmul(x, W) + b
            if not last_layer:
                x = x * (x > 0)  # ReLU
        return x

    def get_topk_from_smi(self, smi='', k=100):
        if not smi:
            return []
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            return []
        return self.get_topk_from_mol(mol, k=k)

    def sigmoid(x):
        """Returns sigmoid of x.

        Args:
            x (float): Input value.
        """
        return 1 / (1 + math.exp(-x))

def softmax(x):
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()

if __name__ == '__main__':
    model = RelevanceTemplatePrioritizer()
    model.load_model()
    smis = ['CCCOCCC', 'CCCNc1ccccc1']
    for smi in smis:
        lst = model.get_topk_from_smi(smi)
        print('{} -> {}'.format(smi, lst))

    # model2 = RelevanceTemplatePrioritizer(use_tf=True)
    # model2.load_model()
    # for smi in smis:
    #     lst = model2.get_topk_from_smi(smi)
    #     print('{} -> {}'.format(smi, lst))

    # import time
    # time.sleep(10)
