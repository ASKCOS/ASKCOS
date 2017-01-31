### Retro transformer
from db import db_client
from django.conf import settings
database = db_client[settings.RETRO_TRANSFORMS['database']]
RETRO_DB = database[settings.RETRO_TRANSFORMS['collection']]
import makeit.retro.transformer as transformer 
RetroTransformer = transformer.Transformer(
	parallel = settings.RETRO_TRANSFORMS['parallel'], 
	nb_workers = settings.RETRO_TRANSFORMS['nb_workers'],
)
mincount_retro = settings.RETRO_TRANSFORMS['mincount']
RetroTransformer.load(RETRO_DB, mincount = mincount_retro, get_retro = True, get_synth = False)
print('Loaded {} retro templates'.format(RetroTransformer.num_templates))
RETRO_FOOTNOTE = 'Using {} retrosynthesis templates (mincount {}) from {}/{}'.format(RetroTransformer.num_templates,
	mincount_retro, settings.RETRO_TRANSFORMS['database'], settings.RETRO_TRANSFORMS['collection'])

### Forward transformer 
database = db_client[settings.SYNTH_TRANSFORMS['database']]
SYNTH_DB = database[settings.SYNTH_TRANSFORMS['collection']]
SynthTransformer = transformer.Transformer()
mincount_synth = settings.SYNTH_TRANSFORMS['mincount']
SynthTransformer.load(SYNTH_DB, mincount = 100000000000000, get_retro = False, get_synth = True)
print('Loaded {} forward templates'.format(SynthTransformer.num_templates))
SYNTH_FOOTNOTE = 'Using {} forward templates (mincount {}) from {}/{}'.format(SynthTransformer.num_templates,
	mincount_synth, settings.SYNTH_TRANSFORMS['database'], settings.SYNTH_TRANSFORMS['collection'])

### Databases
db = db_client[settings.REACTIONS['database']]
REACTION_DB = db[settings.REACTIONS['collection']]
RETRO_LIT_FOOTNOTE = 'Searched {} known reactions from literature'.format(REACTION_DB.count())

db = db_client[settings.INSTANCES['database']]
INSTANCE_DB = db[settings.INSTANCES['collection']]

db = db_client[settings.CHEMICALS['database']]
CHEMICAL_DB = db[settings.CHEMICALS['collection']]

db = db_client[settings.BUYABLES['database']]
BUYABLE_DB = db[settings.BUYABLES['collection']]

db = db_client[settings.SOLVENTS['database']]
SOLVENT_DB = db[settings.SOLVENTS['collection']]

### Prices
print('Loading prices...')
import makeit.retro.pricer as pricer
Pricer = pricer.Pricer()
Pricer.load(CHEMICAL_DB, BUYABLE_DB)
print('Loaded known prices')

### Literaturue transformer
import makeit.retro.transformer_onlyKnown as transformer_onlyKnown
TransformerOnlyKnown = transformer_onlyKnown.TransformerOnlyKnown()
TransformerOnlyKnown.load(CHEMICAL_DB, REACTION_DB)

# Builder
from makeit.webapp.treeBuilder import TreeBuilder 
builder = TreeBuilder(Pricer = Pricer, RetroTransformer = RetroTransformer)

# Intelligent predictor
from makeit.webapp.forwardPredictor import ForwardPredictor 
predictor = ForwardPredictor(nb_workers = 2, TRANSFORM_DB = SYNTH_DB, SOLVENT_DB = SOLVENT_DB)
predictor.load_templates(mincount = mincount_synth)
predictor.load_model(settings.PREDICTOR['trained_model_path'])
PREDICTOR_FOOTNOTE = 'Results generated using {} forward synthetic templates (mincount {}) from {}/{}, scored by a trained machine learning model: '.format(predictor.num_templates,	mincount_synth, settings.SYNTH_TRANSFORMS['database'], settings.SYNTH_TRANSFORMS['collection']) + settings.PREDICTOR['info']

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)