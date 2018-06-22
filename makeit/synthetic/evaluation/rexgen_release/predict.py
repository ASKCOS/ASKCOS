from makeit.synthetic.evaluation.rexgen_release.CandRanker import CandRanker
from makeit.synthetic.evaluation.rexgen_release.CoreFinder import CoreFinder
import rdkit.Chem as Chem 
import sys
import os 

class TFFP():
    '''Template-free forward predictor'''
    def __init__(self):
        froot = os.path.dirname(__file__)
        self.finder = CoreFinder(hidden_size=300, depth=3, batch_size=1)
        print(os.path.join(froot, 'CoreFinder', 'uspto-300-3'))
        self.finder.load_model(os.path.join(froot, 'CoreFinder', 'uspto-300-3'))
        self.ranker = CandRanker(hidden_size=320, depth=3, TOPK=100)
        self.ranker.load_model(os.path.join(froot, 'CandRanker', 'uspto-320-3'))

    def predict(self, smi, top_n=100, num_core=8):
        m = Chem.MolFromSmiles(smi)
        if not m:
            if smi[-1] == '.':
                m = Chem.MolFromSmiles(smi[:-1])
            if not m:
                raise ValueError('Could not parse molecule for TFFP! {}'.format(smi))
        [a.SetIntProp('molAtomMapNumber', i+1) for (i, a) in enumerate(m.GetAtoms())]
        s = Chem.MolToSmiles(m)
        rcores = self.finder.predict([s], num_core=num_core)[0]
        outcomes = self.ranker.predict_one(s, rcores, scores=True, top_n=top_n)
        return(outcomes)


if __name__ == "__main__":
    tffp = TFFP()
    print(tffp.predict('CCCO.CCCBr'))
