from __future__ import print_function
import argparse
from numpy.random import shuffle # for random selection
import rdkit.Chem as Chem          # molecule building
from rdkit.Chem import AllChem
from collections import defaultdict
from rdkit import RDLogger
import datetime # for info files
import json # for dumping
import sys  # for commanad line
import os   # for file paths
import re 
import itertools
from tqdm import tqdm 


if __name__ == '__main__':
    print('Loading transformer...')

    # DATABASE
    from pymongo import MongoClient    # mongodb plugin
    db_client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)

    # Database
    db = db_client['uspto']
    RETRO_DB = db['transforms_retro_v1_allunmapped']

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    # Import retro transformer class and load
    # v3 uses rdchiral, for chiral transformations
    import makeit.webapp.transformer_v3 as transformer
    RetroTransformer = transformer.Transformer()
    RetroTransformer.load(RETRO_DB, mincount=5, get_retro=True, 
        get_synth=False, mincount_chiral=5)
    print('Chiral retro transformer loaded {} retro templates'.format(RetroTransformer.num_templates))

    # Also - get the Pricer and load it into the RetroTransformer object.
    db = db_client['reaxys_v2']
    BUYABLE_DB = db['buyables']
    print('Loading prices...')
    import makeit.retro.pricer as pricer
    RetroTransformer.Pricer = pricer.Pricer()
    RetroTransformer.Pricer.load(None, BUYABLE_DB)
    print('Loaded known prices')

    while True:
        try:
            prompt = raw_input('Enter SMILES: ')
            if prompt.strip() == 'quit':
                break
            smiles = prompt.strip()
            result = RetroTransformer.perform_retro(smiles)
            outcomes = result.return_top(3)
            print(outcomes)

        except Exception as e:
            print(e)