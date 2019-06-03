from makeit.synthetic.evaluation.rexgen_direct.core_wln_global.directcorefinder import DirectCoreFinder
from makeit.synthetic.evaluation.rexgen_direct.rank_diff_wln.directcandranker import DirectCandRanker
import rdkit.Chem as Chem 
import sys
import os 

class TFFP():
    '''Template-free forward predictor'''
    def __init__(self):
        self.finder = DirectCoreFinder(batch_size=1)
        self.finder.load_model()
        self.ranker = DirectCandRanker()
        self.ranker.load_model()

    def predict(self, smi, top_n=100):
        m = Chem.MolFromSmiles(smi)
        if not m:
            if smi[-1] == '.':
                m = Chem.MolFromSmiles(smi[:-1])
            if not m:
                raise ValueError('Could not parse molecule for TFFP! {}'.format(smi))
        [a.SetIntProp('molAtomMapNumber', i+1) for (i, a) in enumerate(m.GetAtoms())]
        s = Chem.MolToSmiles(m)
        (react, bond_preds, bond_scores, cur_att_score) = self.finder.predict(s)
        outcomes = self.ranker.predict(react, bond_preds, bond_scores, scores=True, top_n=top_n)
        return(outcomes)


if __name__ == "__main__":
    tffp = TFFP()
    if len(sys.argv) < 2:
        react = 'CCCO.CCCBr'
    else:
        react = str(sys.argv[1])

    print(react)
    print(tffp.predict(react))