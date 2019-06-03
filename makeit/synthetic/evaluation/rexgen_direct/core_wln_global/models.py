import tensorflow as tf
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.mol_graph import max_nb
from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.nn import *

def gated_convnet(graph_inputs, batch_size=64, hidden_size=300, depth=3, res_block=2):
    input_atom, input_bond, atom_graph, bond_graph, num_nbs, node_mask = graph_inputs
    layers = [input_atom]
    atom_features = input_atom
    for i in range(depth):
        fatom_nei = tf.gather_nd(atom_features, atom_graph)
        fbond_nei = tf.gather_nd(input_bond, bond_graph)
        f_nei = tf.concat([fatom_nei, fbond_nei], 3)
        h_nei = linearND(f_nei, hidden_size, "nei_hidden_%d" % i)
        g_nei = tf.nn.sigmoid(linearND(f_nei, hidden_size, "nei_gate_%d" % i))
        f_nei = h_nei * g_nei
        mask_nei = tf.reshape(tf.sequence_mask(tf.reshape(num_nbs, [-1]), max_nb, dtype=tf.float32), [batch_size,-1,max_nb,1])
        f_nei = tf.reduce_sum(f_nei * mask_nei, -2)
        h_self = linearND(atom_features, hidden_size, "self_hidden_%d" % i)
        g_self = tf.nn.sigmoid(linearND(atom_features, hidden_size, "self_gate_%d" % i))
        f_self = h_self * g_self
        atom_features = (f_nei + f_self) * node_mask
        if res_block is not None and i % res_block == 0 and i > 0:
            atom_features = atom_features + layers[-2]
        layers.append(atom_features)
    output_gate = tf.nn.sigmoid(linearND(atom_features, hidden_size, "out_gate")) 
    output = node_mask * (output_gate * atom_features)
    fp = tf.reduce_sum(output, 1)
    return atom_features * node_mask, fp

def rcnn_wl_last(graph_inputs, batch_size, hidden_size, depth, training=True):
    input_atom, input_bond, atom_graph, bond_graph, num_nbs, node_mask = graph_inputs
    atom_features = tf.nn.relu(linearND(input_atom, hidden_size, "atom_embedding", init_bias=None))
    layers = []
    for i in range(depth):
        with tf.variable_scope("WL", reuse=(i>0)) as scope:
            fatom_nei = tf.gather_nd(atom_features, atom_graph)
            fbond_nei = tf.gather_nd(input_bond, bond_graph)
            h_nei_atom = linearND(fatom_nei, hidden_size, "nei_atom", init_bias=None)
            h_nei_bond = linearND(fbond_nei, hidden_size, "nei_bond", init_bias=None)
            h_nei = h_nei_atom * h_nei_bond
            mask_nei = tf.reshape(tf.sequence_mask(tf.reshape(num_nbs, [-1]), max_nb, dtype=tf.float32), [batch_size,-1,max_nb,1])
            f_nei = tf.reduce_sum(h_nei * mask_nei, -2)
            f_self = linearND(atom_features, hidden_size, "self_atom", init_bias=None)
            layers.append(f_nei * f_self * node_mask)
            l_nei = tf.concat([fatom_nei, fbond_nei], 3)
            nei_label = tf.nn.relu(linearND(l_nei, hidden_size, "label_U2"))
            nei_label = tf.reduce_sum(nei_label * mask_nei, -2) 
            new_label = tf.concat([atom_features, nei_label], 2)
            new_label = linearND(new_label, hidden_size, "label_U1")
            atom_features = tf.nn.relu(new_label)
    #kernels = tf.concat(1, layers)
    kernels = layers[-1]
    fp = tf.reduce_sum(kernels, 1)
    return kernels, fp

