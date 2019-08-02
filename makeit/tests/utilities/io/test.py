import unittest
import makeit.utilities.io.draw as draw
import os
from PIL import Image

class TestDraw(unittest.TestCase):
    def test_01_ReactionStringToImage(self):
        rxn_string = 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1.c1nc[nH]n1.Cl.O=C([O-])O.[Na+]>>OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F'
        result = rxn_image_string = draw.ReactionStringToImage(rxn_string, strip=True)
        expected = Image.open(os.path.join(os.path.dirname(__file__), 'expected/draw_test_rxn_string.png'))
        self.assertEqual(list(expected.getdata()), list(result.getdata()))

    def test_02_ReactionStringToImageRetro(self):
        rxn_string = 'Fc1ccc(C2(Cn3cncn3)CO2)c(F)c1.c1nc[nH]n1.Cl.O=C([O-])O.[Na+]>>OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F'
        result = draw.ReactionStringToImage(rxn_string, strip=True, retro=True)
        expected = Image.open(os.path.join(os.path.dirname(__file__), 'expected/draw_retro_test_rxn_string.png'))
        self.assertEqual(list(expected.getdata()), list(result.getdata()))

    def test_03_TransformStringToImage(self):
        tform = '([O;H0:3]=[C;H0:4](-[C:5])-[NH:2]-[C:1])>>([C:1]-[NH2:2]).([OH:3]-[C;H0:4](=O)-[C:5])'
        result = draw.TransformStringToImage(tform)
        expected = Image.open(os.path.join(os.path.dirname(__file__), 'expected/draw_transform.png'))
        self.assertEqual(list(expected.getdata()), list(result.getdata()))

if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
