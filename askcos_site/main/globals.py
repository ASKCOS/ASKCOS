# Setting logging low
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

### Retro transformer
from database import db_client
from django.conf import settings
database = db_client[settings.RETRO_TRANSFORMS['database']]
RETRO_DB = database[settings.RETRO_TRANSFORMS['collection']]
import makeit.retro.transformer as transformer 
RetroTransformer = transformer.Transformer(
	parallel = False, 
	nb_workers = 1,
)
mincount_retro = settings.RETRO_TRANSFORMS['mincount']
RetroTransformer.load(RETRO_DB, mincount=mincount_retro, get_retro=False, get_synth=False, refs=True)
print('Loaded {} retro templates'.format(RetroTransformer.num_templates))
RETRO_FOOTNOTE = 'Using {} retrosynthesis templates (mincount {}) from {}/{}'.format(RetroTransformer.num_templates,
	mincount_retro, settings.RETRO_TRANSFORMS['database'], settings.RETRO_TRANSFORMS['collection'])

### Chiral Retro Transformer
database = db_client[settings.RETRO_TRANSFORMS_CHIRAL['database']]
RETRO_DB = database[settings.RETRO_TRANSFORMS_CHIRAL['collection']]
import makeit.webapp.transformer_v2 as transformer_v2
RetroTransformerChiral = transformer_v2.Transformer()
mincount_retro = settings.RETRO_TRANSFORMS_CHIRAL['mincount']
mincount_retro_chiral = settings.RETRO_TRANSFORMS_CHIRAL['mincount_chiral']
RetroTransformerChiral.load(RETRO_DB, mincount=mincount_retro, get_retro=False, 
    get_synth=False, refs=True, mincount_chiral=mincount_retro_chiral)
print('Loaded {} retro templates'.format(RetroTransformerChiral.num_templates))
RETRO_CHIRAL_FOOTNOTE = 'Using {} chiral retrosynthesis templates (mincount {} if achiral, mincount {} if chiral) from {}/{}'.format(RetroTransformerChiral.num_templates,
    mincount_retro, mincount_retro_chiral, settings.RETRO_TRANSFORMS_CHIRAL['database'], 
    settings.RETRO_TRANSFORMS_CHIRAL['collection'])
RetroTransformer.templates += RetroTransformerChiral.templates[:]
del RetroTransformerChiral
print('Merged two retrotransformers into one, since this is just for template look-up')
print('{} total templates available'.format(len(RetroTransformer.templates)))
RetroTransformer.reorder() # rebuilds id->template dictionary

### Forward transformer 
database = db_client[settings.SYNTH_TRANSFORMS['database']]
SYNTH_DB = database[settings.SYNTH_TRANSFORMS['collection']]
SynthTransformer = transformer.Transformer()
mincount_synth = settings.SYNTH_TRANSFORMS['mincount']
SynthTransformer.load(SYNTH_DB, mincount = mincount_synth, get_retro = False, get_synth = True)
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

# ### Literaturue transformer
# import makeit.retro.transformer_onlyKnown as transformer_onlyKnown
# TransformerOnlyKnown = transformer_onlyKnown.TransformerOnlyKnown()
# TransformerOnlyKnown.load(CHEMICAL_DB, REACTION_DB)
TransformerOnlyKnown = None 

# Builder
from makeit.webapp.treeBuilder import TreeBuilder 
builder = TreeBuilder(Pricer = Pricer, RetroTransformer = RetroTransformer)

# Intelligent predictor
from makeit.webapp.forwardPredictor import ForwardPredictor 
predictor = ForwardPredictor(nb_workers = 0, TRANSFORM_DB = SYNTH_DB, SOLVENT_DB = SOLVENT_DB)
predictor.load_templates(mincount = mincount_synth)
# predictor.load_model(settings.PREDICTOR['trained_model_path'])
PREDICTOR_FOOTNOTE = 'Results generated using <= {} forward synthetic templates (mincount >= {}) from {}/{}, scored by a trained machine learning model: '.format(predictor.num_templates,	mincount_synth, settings.SYNTH_TRANSFORMS['database'], settings.SYNTH_TRANSFORMS['collection']) + settings.PREDICTOR['info']

# NN predictor
from sklearn.externals import joblib
from askcos_site.functions.nnPredictor import NNConditionPredictor
print('Imported joblib and NNConditionPredictor')
# # Load all the instance IDs from the test model
# rxd_ids = []
# rxn_ids = []
# with open(settings.CONTEXT_REC['info_path'], 'r') as infile:
#     rxn_ids.append(infile.readlines()[1:])  # a list of str(rxn_ids) with '\n'
# for id in rxn_ids[0]:
#     rxd_ids.append(id.replace('\n', ''))
# print('Read context recommendation info file from {}'.format(settings.CONTEXT_REC['info_path']))
# lshf_nn = joblib.load(settings.CONTEXT_REC['model_path'])
# print('Loaded context recommendation nearest-neighbor model from {}'.format(settings.CONTEXT_REC['model_path']))
NN_PREDICTOR = NNConditionPredictor(nn_model=None, rxn_ids=None)
print('Initialized NN_PREDICTOR - unnecessary')
# Keeping track of what reactions have already been done
DONE_SYNTH_PREDICTIONS = {}
