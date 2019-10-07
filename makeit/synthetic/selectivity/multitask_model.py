import tensorflow as tf
from rdkit import Chem
import pickle as pk
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.nn import linearND
from makeit.synthetic.selectivity.mol_graph import atom_fdim as adim, bond_fdim as bdim, max_nb, smiles2graph_list as _s2g
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.models import *
from makeit.synthetic.selectivity.ioutils_direct import * 
from functools import partial
import os


class tf_predictor():
    def __init__(self, depth=5, hidden_size=300, batch_size=1):
        
        self.depth = depth
        self.hidden_size = hidden_size
        self.batch_size = batch_size
        with open(os.path.join(os.path.dirname(__file__), 'task_dict.pkl'), 'rb') as f:
            self.task_dict = pk.load(f)
        self.task_dict_rev = {v: k for k, v in self.task_dict.items()}
        self.saver = None
        self.num_tasks = len(self.task_dict)
        self.save_path = os.path.join(os.path.dirname(__file__), 'model/model.ckpt-30615') 
        self.smiles2graph_batch = partial(_s2g, idxfunc=lambda x:x.GetIdx())
        self.adim = adim
        self.bdim = bdim
        self.max_nb = max_nb
        
    def build(self):
        # Unpack for convenience
        batch_size = self.batch_size 
        hidden_size = self.hidden_size
        adim = self.adim
        bdim = self.bdim
        max_nb = self.max_nb
        depth = self.depth
        
        self.session = tf.Session()

        input_atom = tf.placeholder(tf.float32, [batch_size, None, adim])
        input_bond = tf.placeholder(tf.float32, [batch_size, None, bdim])
        atom_graph = tf.placeholder(tf.int32, [batch_size, None, max_nb, 2])
        bond_graph = tf.placeholder(tf.int32, [batch_size, None, max_nb, 2])
        num_nbs = tf.placeholder(tf.int32, [batch_size, None])
        node_mask = tf.placeholder(tf.float32, [batch_size, None])
        self._src_holder = [input_atom, input_bond, atom_graph, bond_graph, num_nbs, node_mask]
        self._binary = tf.placeholder(tf.float32, [batch_size, None, None, binary_fdim])

        node_mask_exp = tf.expand_dims(node_mask, -1)
        graph_inputs = (input_atom, input_bond, atom_graph, bond_graph, num_nbs, node_mask_exp)
        
        #WLN-NN model
        with tf.variable_scope("encoder"):
            atom_hiddens, _ = rcnn_wl_last(graph_inputs, batch_size=batch_size, 
                                            hidden_size=hidden_size, depth=depth)

        # For each pair of atoms, compute an attention score
        atom_hiddens1 = tf.reshape(atom_hiddens, [batch_size, 1, -1, hidden_size])
        atom_hiddens2 = tf.reshape(atom_hiddens, [batch_size, -1, 1, hidden_size])
        atom_pair = atom_hiddens1 + atom_hiddens2
        att_hidden = tf.nn.relu(linearND(atom_pair, hidden_size, scope="att_atom_feature", init_bias=None) + 
                                linearND(self._binary, hidden_size, scope="att_bin_feature"))
        att_score = linearND(att_hidden, 1, scope="att_scores")
        att_score = tf.nn.sigmoid(att_score)
        
        # Use the attention scores to compute the "att_context" global features
        att_context = att_score * atom_hiddens1
        att_context = tf.reduce_sum(att_context, 2)
            
        # Compute selectivity toward each atom based on the local and global representations
        atom_logits_alltasks = tf.squeeze(linearND(atom_hiddens, self.num_tasks, scope="local_score") + 
                                            linearND(att_context, self.num_tasks, scope="global_score"))
        
        self.atom_likelihoods_smiles = tf.sigmoid(atom_logits_alltasks)

        self.saver = tf.train.Saver()
        self.saver.restore(self.session, self.save_path)

    def web_predictor(self, smiles):
        '''
        Predictor that will give atom scores for a smiles list
        '''
        cur_batch_size = 1
        src_tuple = self.smiles2graph_batch(smiles)
        cur_bin = binary_features_batch(smiles)
        feed_map = {x:y for x,y in zip(self._src_holder, src_tuple)}
        feed_map.update({self._binary:cur_bin})

        alikelihoods = self.session.run(self.atom_likelihoods_smiles, feed_dict=feed_map)
        
        results = []
        for i,j in enumerate(list(zip(*alikelihoods))):
            results.append({
                'smiles':smiles[0],
                'task': self.task_dict_rev[i],
                'atom_scores': tuple([float(x) for x in j])})

        return results





