import makeit.global_config as gc
from makeit.prioritization.prioritizer import Prioritizer
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
from makeit.utilities.io.logging import MyLogger
import math
import sys
import random
import time
import os
import cPickle as pickle
import tensorflow as tf 
import math

relevance_template_prioritizer_loc = 'relevance_template_prioritizer'

def linearND(input_, output_size, scope, reuse=False, init_bias=0.0):
    shape = input_.get_shape().as_list()
    ndim = len(shape)
    stddev = min(1.0 / math.sqrt(shape[-1]), 0.1)
    with tf.variable_scope(scope, reuse=reuse):
        W = tf.get_variable("Matrix", [shape[-1], output_size], tf.float32, tf.random_normal_initializer(stddev=stddev))
    X_shape = tf.gather(tf.shape(input_), range(ndim-1))
    target_shape = tf.concat([X_shape, [output_size]], 0)
    exp_input = tf.reshape(input_, [-1, shape[-1]])
    if init_bias is None:
        res = tf.matmul(exp_input, W)
    else:
        with tf.variable_scope(scope, reuse=reuse):
            b = tf.get_variable("bias", [output_size], initializer=tf.constant_initializer(init_bias))
        res = tf.matmul(exp_input, W) + b
    res = tf.reshape(res, target_shape)
    res.set_shape(shape[:-1] + [output_size])
    return res

class RelevanceTemplatePrioritizer(Prioritizer):
    '''
    Allows to prioritize the templates based on their relevance
    '''

    def __init__(self, retro=True, use_tf=True):
        self.retro = retro
        self.FP_len = 2048
        self.FP_rad = 2
        self.vars = []
        self.template_count = 100
        self.max_cum_prob = 1
        self.batch_size = 1
        self.NK = 100

        if use_tf:
            def load_model(depth=5, hidden_size=300, output_size=61142):
                config = tf.ConfigProto()
                config.gpu_options.allow_growth = True
                self.session = tf.Session(config=config)
                self.input_mol = tf.placeholder(tf.float32, [self.batch_size, self.FP_len])
                self.mol_hiddens = tf.nn.relu(linearND(self.input_mol, hidden_size, scope="encoder0", reuse=None))
                for d in xrange(1, depth):
                    self.mol_hiddens = tf.nn.relu(linearND(self.mol_hiddens, hidden_size, scope="encoder%i"%d, reuse=None))

                self.score = linearND(self.mol_hiddens, output_size, scope="output", reuse=None)
                _, self.topk = tf.nn.top_k(self.score, k=self.NK)

                tf.global_variables_initializer().run(session=self.session)
                size_func = lambda v: reduce(lambda x, y: x*y, v.get_shape().as_list())
                n = sum(size_func(v) for v in tf.trainable_variables())
                print "Model size: %dK" % (n/1000,)

                self.coord = tf.train.Coordinator()
                with open(gc.Relevance_Prioritization['trained_model_path_{}'.format(self.retro)], 'rb') as fid:
                    variables = pickle.load(fid)
                for i, v in enumerate(tf.trainable_variables()):
                    assign_op = tf.assign(v, variables[i])
                    self.session.run(assign_op)
                    del assign_op
                print('Loaded tf model from numpy arrays')

        else:
            def load_model():
                with open(gc.Relevance_Prioritization['trained_model_path_{}'.format(self.retro)], 'rb') as fid:
                    self.vars = pickle.load(fid)
                if gc.DEBUG:
                    MyLogger.print_and_log('Loaded relevance based template prioritization model from {}'.format(
                    gc.Relevance_Prioritization['trained_model_path_{}'.format(self.retro)]), relevance_template_prioritizer_loc)
                return self
        self.load_model = load_model


        if use_tf:
            def get_topk_from_mol(mol, k=100):
                fp = self.mol_to_fp(mol).astype(np.float32).reshape((1, self.FP_len))
                cur_scores, = self.session.run([self.score], feed_dict={
                    self.input_mol: fp,
                })
                indices = list(cur_scores[0,:].argsort()[-k:][::-1])
                cur_scores.sort()
                probs = softmax(cur_scores[0,:])
                return probs[-k:][::-1], indices

        else:
            def get_topk_from_mol(mol, k=100):
                fp = self.mol_to_fp(mol).astype(np.float32)
                cur_scores = self.apply(fp)
                indices = list(cur_scores.argsort()[-k:][::-1])
                cur_scores.sort()
                probs = softmax(cur_scores)
                return probs[-k:][::-1], indices
        self.get_topk_from_mol = get_topk_from_mol

    def mol_to_fp(self, mol):
        if mol is None:
            return np.zeros((self.FP_len,), dtype=np.float32)
        return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, self.FP_rad, nBits=self.FP_len,
                                                              useChirality=True), dtype=np.float32)

    def smi_to_fp(self, smi):
        if not smi:
            return np.zeros((self.FP_len,), dtype=np.float32)
        return self.mol_to_fp(Chem.MolFromSmiles(smi))

    def get_priority(self, input_tuple, **kwargs):
        (templates, target) = input_tuple
        template_count = kwargs.get('template_count', 100)
        max_cum_prob = kwargs.get('max_cum_prob', 0.995)
        # Templates should be sorted by popularity for indices to be correct!
        probs, top_ids = self.get_topk_from_smi(smi=target, k = min(template_count, len(templates)))
        top_templates = []
        cum_score = 0
        for i, id in enumerate(top_ids):
            templates[id]['score'] = probs[i]
            top_templates.append(templates[id])
            cum_score += probs[i]
            #End loop if max cumulative score is exceeded
            if cum_score >= max_cum_prob:
                break
        return top_templates

    # def load_model(self):
    #     with open(gc.Relevance_Prioritization['trained_model_path_{}'.format(self.retro)], 'rb') as fid:
    #         self.vars = pickle.load(fid)
    #     if gc.DEBUG:
    #         MyLogger.print_and_log('Loaded relevance based template prioritization model from {}'.format(
    #         gc.Relevance_Prioritization['trained_model_path_{}'.format(self.retro)]), relevance_template_prioritizer_loc)
    #     return self

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
        return 1 / (1 + math.exp(-x))

def softmax(x):
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()

if __name__ == '__main__':
    model = RelevanceTemplatePrioritizer(use_tf=True)
    model.load_model()
    smis = ['CCCOCCC', 'CCCNc1ccccc1']
    for smi in smis:
        lst = model.get_topk_from_smi(smi)
        print('{} -> {}'.format(smi, lst))

    model2 = RelevanceTemplatePrioritizer(use_tf=True)
    model2.load_model()
    for smi in smis:
        lst = model2.get_topk_from_smi(smi)
        print('{} -> {}'.format(smi, lst))
        
    import time
    time.sleep(10)
