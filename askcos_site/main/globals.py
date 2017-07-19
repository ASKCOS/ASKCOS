# Setting logging low
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
from database import db_client
from django.conf import settings
import cPickle as pickle
import os 

def get_retrotransformer_achiral_path(dbname, collname, mincount_retro):
    return os.path.join(settings.LOCAL_STORAGE['dir'], 
        'retrotransformer_achiral_using_%s-%s_mincount%i.pkl' % (dbname, collname, mincount_retro))

def get_retrotransformer_chiral_path(dbname, collname, mincount_retro, mincount_retro_chiral):
    return os.path.join(settings.LOCAL_STORAGE['dir'], 
        'retrotransformer_chiral_using_%s-%s_mincount%i_mincountchiral%i.pkl' % (dbname, collname, mincount_retro, mincount_retro_chiral))

def get_synthtransformer_path(dbname, collname, mincount):
    return os.path.join(settings.LOCAL_STORAGE['dir'], 
        'synthtransformer_using_%s-%s_mincount%i.pkl' % (dbname, collname, mincount))

def get_pricer_path(chem_dbname, chem_collname, buyable_dbname, buyable_collname):
    return os.path.join(settings.LOCAL_STORAGE['dir'], 
        'pricer_using_%s-%s_and_%s-%s.pkl' % (chem_dbname, chem_collname, buyable_dbname, buyable_collname))

### Retro transformer
import makeit.retro.transformer as transformer 
RetroTransformer = transformer.Transformer()
save_path = get_retrotransformer_achiral_path(
    settings.RETRO_TRANSFORMS['database'],
    settings.RETRO_TRANSFORMS['collection'],
    settings.RETRO_TRANSFORMS['mincount'],
)
if os.path.isfile(save_path):
    with open(save_path, 'rb') as fid:
        RetroTransformer.templates = pickle.load(fid)
    RetroTransformer.reorder()
else:
    database = db_client[settings.RETRO_TRANSFORMS['database']]
    RETRO_DB = database[settings.RETRO_TRANSFORMS['collection']]
    mincount_retro = settings.RETRO_TRANSFORMS['mincount']
    RetroTransformer.load(RETRO_DB, mincount=mincount_retro, 
        get_retro=False, get_synth=False, refs=True)
    print('Saving achiral retro transformer for the (only?) first time')
    with open(save_path, 'wb') as fid:
        pickle.dump(RetroTransformer.templates, fid, -1)
print('Loaded {} retro templates'.format(RetroTransformer.num_templates))
RETRO_FOOTNOTE = 'Using {} retrosynthesis templates (mincount {}) from {}/{}'.format(RetroTransformer.num_templates,
    settings.RETRO_TRANSFORMS['mincount'], settings.RETRO_TRANSFORMS['database'], settings.RETRO_TRANSFORMS['collection'])

### Chiral Retro Transformer
database = db_client[settings.RETRO_TRANSFORMS_CHIRAL['database']]
RETRO_DB = database[settings.RETRO_TRANSFORMS_CHIRAL['collection']]
import makeit.webapp.transformer_v3 as transformer_v3
RetroTransformerChiral = transformer_v3.Transformer()
save_path = get_retrotransformer_chiral_path(
    settings.RETRO_TRANSFORMS_CHIRAL['database'],
    settings.RETRO_TRANSFORMS_CHIRAL['collection'],
    settings.RETRO_TRANSFORMS_CHIRAL['mincount'],
    settings.RETRO_TRANSFORMS_CHIRAL['mincount_chiral'],
)
if os.path.isfile(save_path):
    with open(save_path, 'rb') as fid:
        RetroTransformerChiral.templates = pickle.load(fid)
    RetroTransformerChiral.reorder()
else:
    mincount_retro = settings.RETRO_TRANSFORMS_CHIRAL['mincount']
    mincount_retro_chiral = settings.RETRO_TRANSFORMS_CHIRAL['mincount_chiral']
    RetroTransformerChiral.load(RETRO_DB, mincount=mincount_retro, get_retro=False, 
        get_synth=False, refs=True, mincount_chiral=mincount_retro_chiral)
    print('Saving chiral retro transformer for the (only?) first time')
    with open(save_path, 'wb') as fid:
        pickle.dump(RetroTransformerChiral.templates, fid, -1)
print('Loaded {} retro templates'.format(RetroTransformerChiral.num_templates))
RETRO_CHIRAL_FOOTNOTE = 'Using {} chiral retrosynthesis templates (mincount {} if achiral, mincount {} if chiral) from {}/{}'.format(RetroTransformerChiral.num_templates,
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
SynthTransformer = transformer.Transformer()
save_path = get_synthtransformer_path(
    settings.SYNTH_TRANSFORMS['database'],
    settings.SYNTH_TRANSFORMS['collection'],
    settings.SYNTH_TRANSFORMS['mincount'],
)
if os.path.isfile(save_path):
    with open(save_path, 'rb') as fid:
        SynthTransformer.templates = pickle.load(fid)
    SynthTransformer.reorder()
else:
    database = db_client[settings.SYNTH_TRANSFORMS['database']]
    SYNTH_DB = database[settings.SYNTH_TRANSFORMS['collection']]
    mincount_synth = settings.SYNTH_TRANSFORMS['mincount']
    SynthTransformer.load(SYNTH_DB, mincount = mincount_synth, get_retro = False, get_synth = True)
    print('Loaded {} forward templates'.format(SynthTransformer.num_templates))
    print('Saving synth transformer for the (only?) first time')
    with open(save_path, 'wb') as fid:
        pickle.dump(SynthTransformer.templates, fid, -1)
SYNTH_FOOTNOTE = 'Using {} forward templates (mincount {}) from {}/{}'.format(SynthTransformer.num_templates,
    settings.SYNTH_TRANSFORMS['mincount'], settings.SYNTH_TRANSFORMS['database'], settings.SYNTH_TRANSFORMS['collection'])

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


db = db_client[settings.REACTIONS_OLD['database']]
REACTION_DB_OLD = db[settings.REACTIONS_OLD['collection']]

db = db_client[settings.INSTANCES_OLD['database']]
INSTANCE_DB_OLD = db[settings.INSTANCES_OLD['collection']]

db = db_client[settings.CHEMICALS_OLD['database']]
CHEMICAL_DB_OLD = db[settings.CHEMICALS_OLD['collection']]


### Prices
print('Loading prices...')
import makeit.retro.pricer as pricer
Pricer = pricer.Pricer()
save_path = get_pricer_path(
    settings.CHEMICALS['database'], 
    settings.CHEMICALS['collection'], 
    settings.BUYABLES['database'], 
    settings.BUYABLES['collection'],
)
if os.path.isfile(save_path):
    with open(save_path, 'rb') as fid:
        Pricer.prices = pickle.load(fid) 
        Pricer.prices_flat = pickle.load(fid) 
        Pricer.prices_by_xrn = pickle.load(fid) 
else:
    Pricer.load(CHEMICAL_DB, BUYABLE_DB)
    with open(save_path, 'wb') as fid:
        pickle.dump(Pricer.prices, fid, -1)
        pickle.dump(Pricer.prices_flat, fid, -1)
        pickle.dump(Pricer.prices_by_xrn, fid, -1)
print('Loaded known prices')

# ### Literaturue transformer
# import makeit.retro.transformer_onlyKnown as transformer_onlyKnown
# TransformerOnlyKnown = transformer_onlyKnown.TransformerOnlyKnown()
# TransformerOnlyKnown.load(CHEMICAL_DB, REACTION_DB)
TransformerOnlyKnown = None 


PREDICTOR_FOOTNOTE = 'Results generated using <= {} forward synthetic templates \
(mincount >= {}) from {}/{}, scored by a trained machine learning model: '.format(
    SynthTransformer.num_templates,  settings.SYNTH_TRANSFORMS['mincount'], 
    settings.SYNTH_TRANSFORMS['database'], 
    settings.SYNTH_TRANSFORMS['collection']) + settings.PREDICTOR['info']

# Keeping track of what reactions have already been done
DONE_SYNTH_PREDICTIONS = {}

TEMPLATE_BACKUPS = []
for (dbname, collname) in settings.TEMPLATE_BACKUPS:
    TEMPLATE_BACKUPS.append(db_client[dbname][collname])