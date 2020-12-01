import itertools
import json
from rdkit import Chem
from makeit import global_config as gc


with open(gc.BAN_LIST_PATH) as f:
    ban_list = json.load(f)

BANNED_SMILES = {Chem.MolToSmiles(Chem.MolFromSmiles(smi), isomericSmiles=True)
                 for smi in itertools.chain(*ban_list.values()) if smi is not None}
