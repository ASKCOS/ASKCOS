# Setting logging low
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
from database import db_client
from django.conf import settings
import makeit.utilities.io.pickle as pickle
import os 

# Retro transformer
import makeit.retrosynthetic.transformer as transformer 
RetroTransformer = transformer.RetroTransformer()
RetroTransformer.load(chiral=False, refs=True, rxns=False) 
RetroTransformer.reorder()
RETRO_FOOTNOTE = 'Using {} retrosynthesis templates (mincount {}) from {}/{}'.format(len(RetroTransformer.templates),
    settings.RETRO_TRANSFORMS['mincount'], settings.RETRO_TRANSFORMS['database'], settings.RETRO_TRANSFORMS['collection'])

# Chiral Retro Transformer
import makeit.retrosynthetic.transformer as transformer 
RetroTransformerChiral = transformer.RetroTransformer()
RetroTransformerChiral.load(chiral=True, refs=True, rxns=False) 
RetroTransformerChiral.reorder()
RETRO_CHIRAL_FOOTNOTE = 'Using {} chiral retrosynthesis templates (mincount {} if achiral, mincount {} if chiral) from {}/{}'.format(len(RetroTransformerChiral.templates),
    settings.RETRO_TRANSFORMS_CHIRAL['mincount'], 
    settings.RETRO_TRANSFORMS_CHIRAL['mincount_chiral'], 
    settings.RETRO_TRANSFORMS_CHIRAL['database'], 
    settings.RETRO_TRANSFORMS_CHIRAL['collection'])
RetroTransformer.templates += RetroTransformerChiral.templates[:]
del RetroTransformerChiral
print('Merged two retrotransformers into one, since this is just for template look-up')
print('{} total templates available'.format(len(RetroTransformer.templates)))
RetroTransformer.reorder() # rebuilds id->template dictionary

### Forward transformer 
import makeit.synthetic.enumeration.transformer as transformer
SynthTransformer = transformer.ForwardTransformer()
SynthTransformer.load(chiral=False, rxns=False)
SynthTransformer.reorder()
SYNTH_FOOTNOTE = 'Using {} forward templates (mincount {}) from {}/{}'.format(SynthTransformer.num_templates,
    settings.SYNTH_TRANSFORMS['mincount'], settings.SYNTH_TRANSFORMS['database'], settings.SYNTH_TRANSFORMS['collection'])

### Databases
db = db_client[settings.REACTIONS['database']]
REACTION_DB = db[settings.REACTIONS['collection']]
# RETRO_LIT_FOOTNOTE = 'Searched {} known reactions from literature'.format(REACTION_DB.count())

db = db_client[settings.INSTANCES['database']]
INSTANCE_DB = db[settings.INSTANCES['collection']]

db = db_client[settings.CHEMICALS['database']]
CHEMICAL_DB = db[settings.CHEMICALS['collection']]

db = db_client[settings.BUYABLES['database']]
BUYABLE_DB = db[settings.BUYABLES['collection']]

db = db_client[settings.SOLVENTS['database']]
SOLVENT_DB = db[settings.SOLVENTS['collection']]

db = db_client[settings.REACTIONS_OLD['database']]
REACTION_DB_OLD = db[settings.REACTIONS_OLD['collection']]

db = db_client[settings.INSTANCES_OLD['database']]
INSTANCE_DB_OLD = db[settings.INSTANCES_OLD['collection']]

db = db_client[settings.CHEMICALS_OLD['database']]
CHEMICAL_DB_OLD = db[settings.CHEMICALS_OLD['collection']]


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

TEMPLATE_BACKUPS = []
for (dbname, collname) in settings.TEMPLATE_BACKUPS:
    TEMPLATE_BACKUPS.append(db_client[dbname][collname])


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
from makeit.utilities.io.files import get_abraham_solvents_path
file_path = get_abraham_solvents_path()
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