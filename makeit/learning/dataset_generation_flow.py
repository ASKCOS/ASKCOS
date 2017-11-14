import time
start_time = time.time()
from pymongo import MongoClient
import re
import global_config as gc
balance = 3
import rdkit.Chem as Chem
from utilities.fingerprinting import get_input_condition_as_smiles, get_reaction_as_smiles
from rdkit.Chem import AllChem


def get_all_data():
    print '@{}: Initializing databases...'.format(round((time.time()-start_time)*1000)/1000)
    client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    db = client[gc.REACTIONS['database']]
    reactions = db[gc.REACTIONS['collection']]
    instances = db[gc.INSTANCES['collection']]
    flow_database = db[gc.FLOW_CONDITIONS['collection']]
    chemicals = db[gc.CHEMICALS['collection']]
    flow_database.remove()    

    conditions_flow_reactor = [{'RXD_COND':{'$regex':"[Ff]low [Rr]eactor"}}, {'RXD_COND':{'$regex':"[Ff]low[Rr]eactor"}}]
    
    conditions_continuous_flow = [{'RXD_COND':{'$regex':"[cC]ontinuous [Ff]low"}},{'RXD_COND':{'$regex':"[Cc]ontinuous[Ff]low"}}]
    
    conditions_micro_reactor = [{'RXD_COND':{'$regex':"[Mm]icro[Rr]eactor"}},{'RXD_COND':{'$regex':"[Mm]icro [Rr]eactor"}}]
    
    all_conditions = conditions_flow_reactor + conditions_continuous_flow + conditions_micro_reactor
    flowsearch = re.compile('[Ff]low [Rr]eactor|[Ff]low[Rr]eactor|[Cc]ontinuous [Ff]low|[Cc]ontinuous[Ff]low|[Mm]icro[Rr]eactor|[Mm]icro [Rr]eactor|[Mm]icro-[Rr]eactor')
    instance = instances.find({'$or':all_conditions})
    
    print '@{}: Databases initialized.'.format(round((time.time()-start_time)*1000)/1000)

    reaction_count = 0
    missing_count = 0
    missing_smiles = 0
    encountered = []
    missing_solvent_count = 0
    flow_count = 0
    ex_temp = 0
    flow_cond = 0
    print '@{}: Starting retreival of associated instances...'.format(round((time.time()-start_time)*1000)/1000)

    for doc in instance:
        vars = doc['_id'].split("-")
        var = int(vars[1])
        rx = int(vars[0])
        reaction = reactions.find_one({'_id':rx})
        number = reaction['RX_NVAR']
        #filter out reactions that don't have smiles.
        smiles = get_reaction_as_smiles(doc, reactions, chemicals)
        if smiles == 'NONE':
            missing_smiles += 1
            print '@{}: Reaction {} does not have SMILES information.'.format(round((time.time()-start_time)*1000)/1000, rx)
            continue
        split_smiles = smiles.split('>>')
        for smile in split_smiles:
            if not smile:
                missing_smiles += 1
                print '@{}: Reaction {} does not have SMILES information.'.format(round((time.time()-start_time)*1000)/1000, rx)
                continue
        temp = doc['RXD_T']
        # use lower bound of temperature
        try:
            temp = float(temp)
        except Exception as e:
            if '-' in temp:
                temp = float(temp[0:1+temp.strip()[1:].index('-')].strip())
        solvent = doc['RXD_SOLXRN']
        reagent = doc['RXD_RGTXRN']
        catalyst = doc['RXD_CATXRN']
        solrgtcat = get_input_condition_as_smiles(doc, chemicals, asone = True, check = True)
        conditions_mol = Chem.MolFromSmiles(solrgtcat)
        if (not solvent and not reagent and not catalyst) or not solrgtcat or not conditions_mol:
            #skip this one: we need conditions to assess compatibility.
            print '@{}: Reaction {} does not have conditional information.'.format(round((time.time()-start_time)*1000)/1000, rx)
            flow_cond += 1
            continue
        if temp and temp > 200:
            #skip this one: likely to be large scale process, not micro.
            print '@{}: Reaction {} takes place at too high temperature: {}.'.format(round((time.time()-start_time)*1000)/1000, rx, temp)
            ex_temp +=1
            #skip this one: likely to be large scale process, not micro.
            continue
        #only go through each reaction once.
        if rx not in encountered:
            encountered.append(rx)
            doc['flow'] = True
            flow_database.insert(doc)
            flow_count +=1
            reaction_count += 1
            for i in range(number):
                index = i + 1
                if index == var:
                    continue
                else:
                    tag = "-".join((str(rx),str(index)))
                    associate = instances.find_one({'_id':tag})
                    tag = "-".join((str(rx),str(index)))
                    associate = instances.find_one({'_id':tag})
                    if not associate:
                        #print 'ERROR: Expected associate {} not found'.format(tag)
                        missing_count += 1
                    else:
                        solvent_a = associate['RXD_SOLXRN']
                        reagent_a = associate['RXD_RGTXRN']
                        catalyst_a = associate['RXD_CATXRN']
                        solrgtcat = get_input_condition_as_smiles(associate, chemicals, asone = True, check = True)
                        conditions_mol = Chem.MolFromSmiles(solrgtcat)
                        temp = doc['RXD_T']
                        # use lower bound of temperature
                        try:
                            temp = float(temp)
                        except Exception as e:
                            if '-' in temp:
                                temp = float(temp[0:1+temp.strip()[1:].index('-')].strip()) 
                       
                        if (not solvent_a and not reagent_a and not catalyst_a) or not solrgtcat or not conditions_mol:
                            missing_solvent_count += 1
                        elif temp and temp > 200:
                            ex_temp +=1
                        else:
                            #print '{}, {}'.format(tag, associate)
                            flow = False
                            
                            #extract conditions as a single string.
                            condition = str(associate['RXD_COND'])
    
                            if flowsearch.search(condition):
                                flow = True
                            reaction_count +=1
                            associate['flow_condition'] = solrgtcat
                            associate['flow'] = flow
                            flow_database.insert(associate)
        else:
            print '@{}: Reaction {} has already been encountered.'.format(round((time.time()-start_time)*1000)/1000, rx)


    print 'Total number of reactions added: {}'.format(reaction_count)
    print 'Total number of flow reactions added: {}'.format(flow_count)
    print 'Total number of reactions without smiles specified: {}'.format(missing_smiles)
    print 'Total number of reactions with excessive temperature: {}'.format(ex_temp)
    print 'Total number of reactions without flow conditions specified: {}'.format(flow_cond)
    print 'Total number of associates that are missing: {}'.format(missing_count)
    print "Total number of associates that don't have conditions specified: {}".format(missing_solvent_count)
def get_50_50_data():
    print '@{}: Initializing databases...'.format(round((time.time()-start_time)*1000)/1000)
    client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    db = client[gc.REACTIONS['database']]
    reactions = db[gc.REACTIONS['collection']]
    instances = db[gc.INSTANCES['collection']]
    flow_database = db[gc.FLOW_CONDITIONS_50['collection']]
    chemicals = db[gc.CHEMICALS['collection']]
    flow_database.remove()    

    conditions_flow_reactor = [{'RXD_COND':{'$regex':"[Ff]low [Rr]eactor"}}, {'RXD_COND':{'$regex':"[Ff]low[Rr]eactor"}}]
    
    conditions_continuous_flow = [{'RXD_COND':{'$regex':"[cC]ontinuous [Ff]low"}},{'RXD_COND':{'$regex':"[Cc]ontinuous[Ff]low"}}]
    
    conditions_micro_reactor = [{'RXD_COND':{'$regex':"[Mm]icro[Rr]eactor"}},{'RXD_COND':{'$regex':"[Mm]icro [Rr]eactor"}}]
    
    all_conditions = conditions_flow_reactor + conditions_continuous_flow + conditions_micro_reactor
    flowsearch = re.compile('[Ff]low [Rr]eactor|[Ff]low[Rr]eactor|[Cc]ontinuous [Ff]low|[Cc]ontinuous[Ff]low|[Mm]icro[Rr]eactor|[Mm]icro [Rr]eactor|[Mm]icro-[Rr]eactor')
    instance = instances.find({'$or':all_conditions})
    
    print '@{}: Databases initialized.'.format(round((time.time()-start_time)*1000)/1000)

    reaction_count = 0
    missing_count = 0
    missing_smiles = 0
    encountered = []
    missing_solvent_count = 0
    flow_count = 0
    ex_temp = 0
    flow_cond = 0
    print '@{}: Starting retrieval of associated instances...'.format(round((time.time()-start_time)*1000)/1000)
    
    for doc in instance:
        vars = doc['_id'].split("-")
        var = int(vars[1])
        rx = int(vars[0])
        reaction = reactions.find_one({'_id':rx})
        number = reaction['RX_NVAR']
        #filter out reactions that don't have smiles.
        smiles = get_reaction_as_smiles(doc, reactions, chemicals)
        if smiles == 'NONE':
            missing_smiles += 1
            print '@{}: Reaction {} does not have SMILES information.'.format(round((time.time()-start_time)*1000)/1000, rx)
            continue
        split_smiles = smiles.split('>>')
        for smile in split_smiles:
            if not smile:
                missing_smiles += 1
                print '@{}: Reaction {} does not have SMILES information.'.format(round((time.time()-start_time)*1000)/1000, rx)
                continue
        temp = doc['RXD_T']
        # use lower bound of temperature
        try:
            temp = float(temp)
        except Exception as e:
            if '-' in temp:
                temp = float(temp[0:1+temp.strip()[1:].index('-')].strip())
        solvent = doc['RXD_SOLXRN']
        reagent = doc['RXD_RGTXRN']
        catalyst = doc['RXD_CATXRN']
        solrgtcat = get_input_condition_as_smiles(doc, chemicals, asone = True, check = True)
        conditions_mol = Chem.MolFromSmiles(solrgtcat)
        if (not solvent and not reagent and not catalyst) or not solrgtcat or not conditions_mol:
            #skip this one: we need conditions to assess compatibility.
            print '@{}: Reaction {} does not have conditional information.'.format(round((time.time()-start_time)*1000)/1000, rx)
            flow_cond += 1
            continue
        if temp and temp > 200:
            #skip this one: likely to be large scale process, not micro.
            print '@{}: Reaction {} takes place at too high temperature: {}.'.format(round((time.time()-start_time)*1000)/1000, rx, temp)
            ex_temp +=1
            #skip this one: likely to be large scale process, not micro.
            continue
        #only go through each reaction once.
        if rx not in encountered:
            encountered.append(rx)
            doc['flow'] = True
            flow_database.insert(doc)
            flow_count += 1
            reaction_count += 1
            #try to get to 50-50 ratio.
            flow_entry_balance = balance
            for i in range(number):
                index = i + 1
                if index == var:
                    continue
                else:
                    tag = "-".join((str(rx),str(index)))
                    associate = instances.find_one({'_id':tag})
                    if not associate:
                        #print 'ERROR: Expected associate {} not found'.format(tag)
                        missing_count += 1
                    
                    else:
                        solvent_a = associate['RXD_SOLXRN']
                        reagent_a = associate['RXD_RGTXRN']
                        catalyst_a = associate['RXD_CATXRN']
                        solrgtcat = get_input_condition_as_smiles(associate, chemicals, asone = True, check = True)
                        conditions_mol = Chem.MolFromSmiles(solrgtcat)
                        temp = doc['RXD_T']
                        # use lower bound of temperature
                        try:
                            temp = float(temp)
                        except Exception as e:
                            if '-' in temp:
                                temp = float(temp[0:1+temp.strip()[1:].index('-')].strip()) 
                       
                        if (not solvent_a and not reagent_a and not catalyst_a) or not solrgtcat or not conditions_mol:
                            missing_solvent_count += 1
                        elif temp and temp > 200:
                            ex_temp +=1
                        else:
                            #print '{}, {}'.format(tag, associate)
                            flow = False
                            
                            #extract conditions as a single string.
                            condition = str(associate['RXD_COND'])
                            if flow_entry_balance >= 0:
                                if flowsearch.search(condition):
                                    flow_entry_balance += 1
                                    flow = True
                                else:
                                    flow_entry_balance -= 1
                                
                                #Info exceeds doc size
                                #associate['flow_condition'] = solrgtcat
                                associate['flow'] = flow
                                #only add entry if equal or more flow entries present.
                                flow_database.insert(associate)
                                reaction_count += 1
        else:
            print '@{}: Reaction {} has already been encountered.'.format(round((time.time()-start_time)*1000)/1000, rx)


    print 'Total number of reactions added: {}'.format(reaction_count)
    print 'Total number of flow reactions added: {}'.format(flow_count)
    print 'Total number of reactions without smiles specified: {}'.format(missing_smiles)
    print 'Total number of reactions with excessive temperature: {}'.format(ex_temp)
    print 'Total number of reactions without flow conditions specified: {}'.format(flow_cond)
    print 'Total number of associates that are missing: {}'.format(missing_count)
    print "Total number of associates that don't have conditions specified: {}".format(missing_solvent_count)
if __name__ == '__main__':
    #get_50_50_data()
    get_all_data()