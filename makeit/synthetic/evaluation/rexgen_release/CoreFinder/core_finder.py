import tensorflow as tf
from ..utils.nn import linearND, linear
from .mol_graph import atom_fdim as adim, bond_fdim as bdim, max_nb, smiles2graph_batch
from .models import *
from .ioutils import *

class CoreFinder(object):

    def __init__(self, hidden_size, depth, batch_size=20):
        self.hidden_size = hidden_size
        self.batch_size = batch_size
        self.depth = depth

    def load_model(self, model_path):
        hidden_size = self.hidden_size
        batch_size = self.batch_size
        depth = self.depth
        
        self.graph = tf.Graph()
        with self.graph.as_default():
            input_atom = tf.placeholder(tf.float32, [batch_size, None, adim])
            input_bond = tf.placeholder(tf.float32, [batch_size, None, bdim])
            atom_graph = tf.placeholder(tf.int32, [batch_size, None, max_nb, 2])
            bond_graph = tf.placeholder(tf.int32, [batch_size, None, max_nb, 2])
            num_nbs = tf.placeholder(tf.int32, [batch_size, None])
            node_mask = tf.placeholder(tf.float32, [batch_size, None])
            binary = tf.placeholder(tf.float32, [batch_size, None, None, binary_fdim])
            validity = tf.placeholder(tf.float32, [batch_size, None, None])
            core_size = tf.placeholder(tf.int32)

            self.leaf_nodes = [input_atom, input_bond, atom_graph, bond_graph, num_nbs, node_mask, binary, validity, core_size]

            node_mask = tf.expand_dims(node_mask, -1)
            graph_inputs = (input_atom, input_bond, atom_graph, bond_graph, num_nbs, node_mask)
            with tf.variable_scope("encoder"):
                atom_hiddens, _ = rcnn_wl_last(graph_inputs, batch_size=batch_size, hidden_size=hidden_size, depth=depth)

            atom_hiddens1 = tf.reshape(atom_hiddens, [batch_size, 1, -1, hidden_size])
            atom_hiddens2 = tf.reshape(atom_hiddens, [batch_size, -1, 1, hidden_size])
            atom_pair = atom_hiddens1 + atom_hiddens2

            att_hidden = tf.nn.relu(linearND(atom_pair, hidden_size, scope="att_atom_feature", init_bias=None) + linearND(binary, hidden_size, scope="att_bin_feature"))
            att_score = linearND(att_hidden, 1, scope="att_scores")
            att_score = tf.nn.sigmoid(att_score)
            att_context = att_score * atom_hiddens1
            att_context = tf.reduce_sum(att_context, 2)

            att_context1 = tf.reshape(att_context, [batch_size, 1, -1, hidden_size])
            att_context2 = tf.reshape(att_context, [batch_size, -1, 1, hidden_size])
            att_pair = att_context1 + att_context2

            pair_hidden = linearND(atom_pair, hidden_size, scope="atom_feature", init_bias=None) + linearND(binary, hidden_size, scope="bin_feature", init_bias=None) + linearND(att_pair, hidden_size, scope="ctx_feature")
            pair_hidden = tf.nn.relu(pair_hidden)
            score = tf.squeeze(linearND(pair_hidden, 1, scope="scores"), [3]) + validity * 10000

            score = tf.reshape(score, [batch_size, -1])
            _, self.topk = tf.nn.top_k(score, core_size)

            self.session = tf.Session()
            #tf.global_variables_initializer().run(session=self.session)
            saver = tf.train.Saver()
            saver.restore(self.session, tf.train.latest_checkpoint(model_path))

    def predict(self, reactants, num_core):
        reaction_cores = []
        batch_size = self.batch_size
        num_core *= 2

        for it in xrange(0, len(reactants), batch_size):
            src_batch = reactants[it:it + batch_size]
            src_tuple = smiles2graph_batch(src_batch)
            cur_bin, cur_validity = get_all_batch(src_batch)
            leaf_values = src_tuple + (cur_bin, cur_validity, num_core)
            feed_map = {x:y for x,y in zip(self.leaf_nodes, leaf_values)}
            cur_topk = self.session.run(self.topk, feed_dict=feed_map)
            cur_dim = cur_validity.shape[1]
            
            for i in xrange(batch_size):
                res = []
                for j in xrange(num_core):
                    k = cur_topk[i,j]
                    x = k / cur_dim
                    y = k % cur_dim
                    if x < y and cur_validity[i,x,y] == 1:
                        res.append( (x + 1,y + 1) )
                reaction_cores.append(res)

        return reaction_cores

if __name__ == "__main__":
    import sys
    cf = CoreFinder(core_size=10, hidden_size=300, depth=3)
    cf.load_model("uspto-300-3")
    data = []
    for line in sys.stdin:
        data.append(line.split()[0].split('>')[0])
        if len(data) == 40:
            break
    rcores = cf.predict(data, 10)
    for core in rcores:
        print(core)

