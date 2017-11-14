import time
start_time = time.time()
from pymongo import MongoClient
from bson import ObjectId
import re
import global_config as gc
balance = 3
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from utilities.fingerprinting import get_input_condition
from utilities.strings import string_or_range_to_float

template_id = ObjectId('57bc470dfbff5063fbb7b6fd')

def get_usable_data():
    client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    db = client[gc.INSTANCES['database']]
    instances = db[gc.INSTANCES['collection']]
    chemicals = client[gc.CHEMICALS['database']][gc.CHEMICALS['collection']]
    yield_database = db[gc.YIELDS['collection']]
    yield_database.remove() 
    forward_temps = client[gc.SYNTH_TRANSFORMS['database']][gc.SYNTH_TRANSFORMS['collection']]
    reaction_count = 0
    failed_criteria = 0
    missing = 0
    template = forward_temps.find_one({'_id': ObjectId('57bc470dfbff5063fbb7b6fd')})
    for id in template['references']:
        instance = instances.find_one({'_id':id})
        if not instance:
            print 'Instance {} not found'.format(id)
            missing += 1
            continue
        solvent = instance['RXD_SOLXRN']
        reagent = instance['RXD_RGTXRN']
        catalyst = instance['RXD_CATXRN']
        solrgtcat = get_input_condition(instance, chemicals, asone = True)
        conditions_mol = Chem.MolFromSmiles(solrgtcat)
        not_criteria = not instance['RXD_YD'] \
                        or not string_or_range_to_float(instance['RXD_T']) \
                        or (not solvent and not reagent and not catalyst) \
                        or not string_or_range_to_float(instance['RXD_TIM']) \
                        or not conditions_mol
        
        if not_criteria:
            failed_criteria += 1
        else:
            yield_database.insert(instance)
            print 'Added instance {}'.format(id)
            reaction_count +=1
                   
    print 'Total number of reactions added: {}'.format(reaction_count)
    print 'Total number of reactions left out for criteria: {}'.format(failed_criteria)
    print 'Total number of reactions missing: {}'.format(missing)
    '''
    print 'Total number of flow reactions added: {}'.format(flow_count)
    print 'Total number of reactions with excessive temperature: {}'.format(ex_temp)
    print 'Total number of reactions without flow conditions specified: {}'.format(flow_cond)
    print 'Total number of associates that are missing: {}'.format(missing_count)
    print "Total number of associates that don't have conditions specified: {}".format(missing_solvent_count)
    '''
if __name__ == '__main__':
    get_usable_data()
    '''
    client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    db = client[gc.INSTANCES['database']]
    instances = db[gc.INSTANCES['collection']]
    i = instances.find({'RX_ID':280360})
    print i.count()
    '''