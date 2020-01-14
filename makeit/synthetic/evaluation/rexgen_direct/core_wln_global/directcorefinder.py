import tensorflow as tf
from makeit import global_config as gc
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.nn import linearND, linear
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.models import *
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.ioutils_direct import *
import math, sys, random
from collections import Counter
from optparse import OptionParser
from functools import partial
import threading
from multiprocessing import Queue
import os

NK3 = 80
batch_size = 2 # just fake it, make two 
hidden_size = 300
depth = 3
model_path = gc.TEMPLATE_FREE_FORWARD_PREDICTOR['core_model_path']
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.mol_graph import atom_fdim as adim, bond_fdim as bdim, max_nb, smiles2graph_list as _s2g
smiles2graph_batch = partial(_s2g, idxfunc=lambda x:x.GetIntProp('molAtomMapNumber') - 1)

class DirectCoreFinder():
    def __init__(self, hidden_size=hidden_size, batch_size=batch_size, 
            depth=depth):
        self.hidden_size = hidden_size 
        self.batch_size = batch_size 
        self.depth = depth 

    def load_model(self, model_path=model_path):
        hidden_size = self.hidden_size 
        vbatch_size = self.batch_size 
        depth = self.depth 

        self.graph = tf.Graph()
        with self.graph.as_default():
            input_atom = tf.placeholder(tf.float32, [batch_size, None, adim])
            input_bond = tf.placeholder(tf.float32, [batch_size, None, bdim])
            atom_graph = tf.placeholder(tf.int32, [batch_size, None, max_nb, 2])
            bond_graph = tf.placeholder(tf.int32, [batch_size, None, max_nb, 2])
            num_nbs = tf.placeholder(tf.int32, [batch_size, None])
            node_mask = tf.placeholder(tf.float32, [batch_size, None])
            self.src_holder = [input_atom, input_bond, atom_graph, bond_graph, num_nbs, node_mask]
            self.label = tf.placeholder(tf.int32, [batch_size, None])
            self.binary = tf.placeholder(tf.float32, [batch_size, None, None, binary_fdim])        

            node_mask = tf.expand_dims(node_mask, -1)

            graph_inputs = (input_atom, input_bond, atom_graph, bond_graph, num_nbs, node_mask)
            with tf.variable_scope("encoder"):
                atom_hiddens, _ = rcnn_wl_last(graph_inputs, batch_size=batch_size, hidden_size=hidden_size, depth=depth)

            atom_hiddens1 = tf.reshape(atom_hiddens, [batch_size, 1, -1, hidden_size])
            atom_hiddens2 = tf.reshape(atom_hiddens, [batch_size, -1, 1, hidden_size])
            atom_pair = atom_hiddens1 + atom_hiddens2

            att_hidden = tf.nn.relu(linearND(atom_pair, hidden_size, scope="att_atom_feature", init_bias=None) + linearND(self.binary, hidden_size, scope="att_bin_feature"))
            att_score = linearND(att_hidden, 1, scope="att_scores")
            att_score = tf.nn.sigmoid(att_score)
            att_context = att_score * atom_hiddens1
            att_context = tf.reduce_sum(att_context, 2)

            att_context1 = tf.reshape(att_context, [batch_size, 1, -1, hidden_size])
            att_context2 = tf.reshape(att_context, [batch_size, -1, 1, hidden_size])
            att_pair = att_context1 + att_context2

            pair_hidden = linearND(atom_pair, hidden_size, scope="atom_feature", init_bias=None) + linearND(self.binary, hidden_size, scope="bin_feature", init_bias=None) + linearND(att_pair, hidden_size, scope="ctx_feature")
            pair_hidden = tf.nn.relu(pair_hidden)
            pair_hidden = tf.reshape(pair_hidden, [batch_size, -1, hidden_size])

            score = linearND(pair_hidden, 5, scope="scores")
            score = tf.reshape(score, [batch_size, -1])
            bmask = tf.to_float(tf.equal(self.label, INVALID_BOND)) * 10000
            topk_scores, topk = tf.nn.top_k(score - bmask, k=NK3)
            label_dim = tf.shape(self.label)[1]
            
            # What will be used for inference?
            self.predict_vars = [topk, topk_scores, label_dim, att_score]
            
            # Restore
            self.session = tf.Session()
            saver = tf.train.Saver()
            saver.restore(self.session, model_path)
        
    def predict(self, reactants_smi):

        bo_to_index  = {0.0: 0, 1.0:1, 2.0:2, 3.0:3, 1.5:4}
        bindex_to_o = {val:key for key, val in bo_to_index.items()}
        nbos = len(bo_to_index)

        src_batch, edit_batch = [], []
        m = Chem.MolFromSmiles(reactants_smi)

        if any(not a.HasProp('molAtomMapNumber') for a in m.GetAtoms()):
            mapnum = 1
            for a in m.GetAtoms():
                a.SetIntProp('molAtomMapNumber', mapnum)
                mapnum += 1
        react = Chem.MolToSmiles(m)

        src_batch.append(react)
        src_batch.append(react)
        edit_batch.append('0-1-0.0') # dummy edits
        edit_batch.append('0-1-0.0') # dummy edits

        src_tuple = smiles2graph_batch(src_batch)
        cur_bin, cur_label, sp_label = get_all_batch(zip(src_batch, edit_batch))
        feed_map = {x:y for x,y in zip(self.src_holder, src_tuple)}
        feed_map.update({self.label:cur_label, self.binary:cur_bin})

        cur_topk, cur_sco, cur_dim, cur_att_score = self.session.run(self.predict_vars,
            feed_dict=feed_map)
        cur_dim = int(math.sqrt(cur_dim/5)) # important! changed to divide by 5

        cur_topk = cur_topk[0,:]
        cur_sco = cur_sco[0]
        cur_att_score = cur_att_score[0, :, :]

        bond_preds = []
        bond_scores = []

        # NOTE: we don't filter out molecules known to be reagents, but during training, 
        # molecules known to be reagents/solvents are not allowed to be involved with bond
        # changes.

        for j in range(NK3):
            k = cur_topk[j]
            bindex = k % nbos
            y = ((k - bindex) / nbos) % cur_dim + 1
            x = (k - bindex - (y-1) * nbos) / cur_dim / nbos + 1
            if x < y: # keep canonical
                # x = k / cur_dim + 1 # was for 2D case
                # y = k % cur_dim + 1 # was for 2D case
                bo = bindex_to_o[bindex]
                bond_preds.append("{}-{}-{:.1f}".format(x, y, bo))
                bond_scores.append(cur_sco[j])

        return (react, bond_preds, bond_scores, cur_att_score)
