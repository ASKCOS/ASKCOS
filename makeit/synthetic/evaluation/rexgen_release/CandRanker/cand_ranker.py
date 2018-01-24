import rdkit
import rdkit.Chem as Chem
import tensorflow as tf
from ..utils.nn import linearND, linear
from .mol_graph import atom_fdim as adim, bond_fdim as bdim, max_nb, smiles2graph_test, bond_types
from .models import *
from edit_mol import edit_mol

class CandRanker(object):

    def __init__(self, hidden_size, depth, MAX_NCAND=2000, TOPK=5):
        self.hidden_size = hidden_size
        self.depth = depth
        self.MAX_NCAND = MAX_NCAND
        self.TOPK = TOPK
        
    def load_model(self, model_path):
        hidden_size = self.hidden_size
        depth = self.depth
        
        self.graph = tf.Graph()
        with self.graph.as_default():
            input_atom = tf.placeholder(tf.float32, [None, None, adim])
            input_bond = tf.placeholder(tf.float32, [None, None, bdim])
            atom_graph = tf.placeholder(tf.int32, [None, None, max_nb, 2])
            bond_graph = tf.placeholder(tf.int32, [None, None, max_nb, 2])
            num_nbs = tf.placeholder(tf.int32, [None, None])
            self.leaf_nodes = [input_atom, input_bond, atom_graph, bond_graph, num_nbs]

            graph_inputs = (input_atom, input_bond, atom_graph, bond_graph, num_nbs) 
            with tf.variable_scope("encoder"):
                _, fp = rcnn_wl_last(graph_inputs, hidden_size=hidden_size, depth=depth)

            reactant = fp[0:1,:]
            candidates = fp[1:,:]
            candidates = candidates - reactant
            candidates = linear(candidates, hidden_size, "candidate")
            match = tf.nn.relu(candidates)
            self.score = tf.squeeze(linear(match, 1, "score"), [1])
            self.scaled_score = tf.nn.softmax(self.score)
            cur_k = tf.minimum(self.TOPK, tf.shape(self.score)[0])
            _, self.topk_cands = tf.nn.top_k(self.score, cur_k)

            self.session = tf.Session()
            #tf.global_variables_initializer().run(session=self.session)
            saver = tf.train.Saver()
            saver.restore(self.session, tf.train.latest_checkpoint(model_path))

    def predict_one(self, r, core, scores=True, top_n=100):
        core = [(x-1,y-1) for x,y in core]
        ncore = len(core)
        while True:
            src_tuple,core_conf = smiles2graph_test(r, core[:ncore])
            if len(core_conf) <= self.MAX_NCAND:
                break
            ncore -= 1
        feed_map = {x:y for x,y in zip(self.leaf_nodes, src_tuple)}
        if scores:
            (cur_scores, cur_probs, candidates) = self.session.run([self.score, self.scaled_score, self.topk_cands], feed_dict=feed_map)
        else:
            candidates = self.session.run(self.topk_cands, feed_dict=feed_map)

        rmol = Chem.MolFromSmiles(r)
        rbonds = {}
        for bond in rmol.GetBonds():
            a1 = bond.GetBeginAtom().GetAtomMapNum()
            a2 = bond.GetEndAtom().GetAtomMapNum()
            t = bond_types.index(bond.GetBondType()) + 1
            a1,a2 = min(a1,a2),max(a1,a2)
            rbonds[(a1,a2)] = t

        cand_smiles = []; cand_scores = []; cand_probs = [];
        for idx in candidates:
            edits = []
            for x,y,t in core_conf[idx]:
                x,y = x+1,y+1
                if ((x,y) not in rbonds and t > 0) or ((x,y) in rbonds and rbonds[(x,y)] != t):
                    edits.append( (x,y,t) )
            cand = edit_mol(rmol, edits)
            cand_smiles.append(cand)
            cand_scores.append(cur_scores[idx])
            cand_probs.append(cur_probs[idx])
        
        outcomes = []
        if scores:
            for i in range(min(len(cand_smiles), top_n)):
                outcomes.append({
                    'rank': i + 1,
                    'smiles': cand_smiles[i],
                    'score': cand_scores[i],
                    'prob': cand_probs[i],
                })
        else:
            for i in range(min(len(cand_smiles), top_n)):
                outcomes.append({
                    'rank': i + 1,
                    'smiles': cand_smiles[i],
                })

        return outcomes

