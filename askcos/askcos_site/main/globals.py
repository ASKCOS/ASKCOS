# Setting logging low
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
from database import db_client
from django.conf import settings
import makeit.utilities.io.pickle as pickle
import os
import makeit.global_config as gc

# Chiral Retro Transformer
import makeit.retrosynthetic.transformer as transformer
RetroTransformer = transformer.RetroTransformer(template_prioritizer=None, precursor_prioritizer=None, fast_filter=None)
RetroTransformer.load()
RETRO_CHIRAL_FOOTNOTE = 'Using {} chiral retrosynthesis templates from {}/{}'.format(
    gc.RELEVANCE_TEMPLATE_PRIORITIZATION['reaxys']['output_size'],
    gc.RETRO_TEMPLATES['database'],
    gc.RETRO_TEMPLATES['collection']
)

### Databases
db = db_client[gc.REACTIONS['database']]
REACTION_DB = db[gc.REACTIONS['collection']]

db = db_client[gc.CHEMICALS['database']]
CHEMICAL_DB = db[gc.CHEMICALS['collection']]

db = db_client[gc.BUYABLES['database']]
BUYABLE_DB = db[gc.BUYABLES['collection']]

db = db_client[gc.SOLVENTS['database']]
SOLVENT_DB = db[gc.SOLVENTS['collection']]


### Prices
print('Loading prices...')
import makeit.utilities.buyable.pricer as pricer
Pricer = pricer.Pricer()
Pricer.load()
print('Loaded known prices')

TransformerOnlyKnown = None

PREDICTOR_FOOTNOTE = ''

# Keeping track of what reactions have already been done
DONE_SYNTH_PREDICTIONS = {}

# TEMPLATE_BACKUPS = []
# for (dbname, collname) in settings.TEMPLATE_BACKUPS:
#     TEMPLATE_BACKUPS.append(db_client[dbname][collname])


# Historian
# from makeit.utilities.historian.reactions import ReactionHistorian
# reactionhistorian = ReactionHistorian()
# reactionhistorian.load_from_file()
reactionhistorian = None

# from makeit.utilities.historian.chemicals import ChemHistorian
# chemhistorian = ChemHistorian()
# chemhistorian.load_from_file()
chemhistorian = None

from makeit.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
scscorer = SCScorePrecursorPrioritizer()
scscorer.load_model(model_tag='1024bool')
print('Loaded SCScorer on website')
print(scscorer.get_score_from_smiles('CCCC', noprice=True))


# Solvent choices - the save file is created by the template-based forward predictor
solvent_choices = []
file_path = gc.SOLVENTS['file_name']
if os.path.isfile(file_path):
    with open(file_path, 'rb') as fid:
        solvent_name_to_smiles = pickle.load(fid)
    solvent_choices = [{'smiles': v, 'name': k} for (k, v) in solvent_name_to_smiles.items()]
else:
    db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                    'id'], connect=gc.MONGO['connect'])
    db = db_client[gc.SOLVENTS['database']]
    SOLVENT_DB = db[gc.SOLVENTS['collection']]
    for doc in SOLVENT_DB.find({'_id': {'$ne': 'default'}}):
        solvent_choices.append({
            'smiles': doc['smiles'],
            'name': doc['name'],
        })
