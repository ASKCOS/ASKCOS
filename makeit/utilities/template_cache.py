from cachetools import LFUCache
from pymongo import MongoClient
from bson.objectid import ObjectId
from rdchiral.initialization import rdchiralReaction
from rdkit.Chem import AllChem
from makeit.utilities.io.logger import MyLogger
import makeit.global_config as gc


def doc_to_template(document, chiral):
    """Returns a template given a document from the database or file.

    Args:
        document (dict): Document of template from database or file.
        chiral (bool): Whether to pay attention to chirality.

    Returns:
        dict: Retrosynthetic template.
    """
    if 'reaction_smarts' not in document:
        return
    reaction_smarts = str(document['reaction_smarts'])
    if not reaction_smarts:
        return

    # different thresholds for chiral and non chiral reactions
    chiral_rxn = False
    for c in reaction_smarts:
        if c in ('@', '/', '\\'):
            chiral_rxn = True
            break

    # Define dictionary
    template = {
        'name':                 document['name'] if 'name' in document else '',
        'reaction_smarts':      reaction_smarts,
        'incompatible_groups':  document['incompatible_groups'] if 'incompatible_groups' in document else [],
        'reference':            document['reference'] if 'reference' in document else '',
        'references':           document['references'] if 'references' in document else [],
        'rxn_example':          document['rxn_example'] if 'rxn_example' in document else '',
        'explicit_H':           document['explicit_H'] if 'explicit_H' in document else False,
        '_id':                  document['_id'] if '_id' in document else -1,
        'product_smiles':       document['product_smiles'] if 'product_smiles' in document else [],
        'necessary_reagent':    document['necessary_reagent'] if 'necessary_reagent' in document else '',
        'efgs':                 document['efgs'] if 'efgs' in document else None,
        'intra_only':           document['intra_only'] if 'intra_only' in document else False,
        'dimer_only':           document['dimer_only'] if 'dimer_only' in document else False,
    }
    template['chiral'] = chiral_rxn

    # Frequency/popularity score
    template['count'] = document.get('count', 1)

    # Define reaction in RDKit and validate
    try:
        # Force reactants and products to be one pseudo-molecule (bookkeeping)
        reaction_smarts_one = '(' + reaction_smarts.replace('>>', ')>>(') + ')'

        if chiral:
            rxn = rdchiralReaction(str(reaction_smarts_one))
            template['rxn'] = rxn
        else:
            rxn = AllChem.ReactionFromSmarts(
                str(reaction_smarts_one))
            if rxn.Validate()[1] == 0:
                template['rxn'] = rxn
            else:
                template['rxn'] = None

    except Exception as e:
        if gc.DEBUG:
            MyLogger.print_and_log('Couldnt load : {}: {}'.format(
                reaction_smarts_one, e), 'template_cache', level=1)
        template['rxn'] = None
        template['rxn_f'] = None
    return template

class TemplateCache(LFUCache):
    def __init__(self, maxsize, chiral=True, getsizeof=None):
        self.TEMPLATE_DB = None
        self.chiral = chiral
        super(TemplateCache, self).__init__(maxsize, getsizeof=getsizeof)

    def __missing__(self, key):
        if self.TEMPLATE_DB is None:
            db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                                    'id'], connect=gc.MONGO['connect'])

            db_name = gc.RETRO_TRANSFORMS_CHIRAL['database']
            collection = gc.RETRO_TRANSFORMS_CHIRAL['collection']
            self.TEMPLATE_DB = db_client[db_name][collection]

        doc = self.TEMPLATE_DB.find_one({'_id': ObjectId(key)})
        template = doc_to_template(doc, self.chiral)
        self[key] = template
        return template
