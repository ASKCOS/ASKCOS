import global_config as gc
from pymongo import MongoClient
from utilities.fingerprinting import get_input_condition_as_smiles
client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
db1 = client[gc.SOLVENTS['database']]
solvents = db1[gc.SOLVENTS['collection']]
reactions = client[gc.REACTIONS['database']][gc.REACTIONS['collection']]
db2 = client[gc.FLOW_CONDITIONS2['database']]
flow_database = db2[gc.FLOW_CONDITIONS2['collection']]
chemicals = db2[gc.CHEMICALS['collection']]
all = flow_database.find()
encounters = 0
total = 0
i=0
corr = 0
wrong = 0
for one in all:
    a = one['RXD_SOLXRN']
    b = one['RXD_RGTXRN']
    c = one['RXD_CATXRN']
    for id in a:
        four = '++++'
        three = '+++'
        subs = chemicals.find_one({'_id':id})
        try:
            print subs['SMILES_new']
            if four in subs['SMILES_new']:
                updated = subs['SMILES_new'].replace(four, '+4')
                print updated
                raw_input("Update OK?")
                chemicals.update({'_id':id},{'$set':{'SMILES_new':updated}})
            if three in subs['SMILES_new']:
                updated = subs['SMILES_new'].replace(three, '+3')
                print updated
                raw_input("Update OK?")
                chemicals.update({'_id':id},{'$set':{'SMILES_new':updated}})
        except:
            pass
    for id in b:
        four = '++++'
        three = '+++'
        subs = chemicals.find_one({'_id':id})
        try:
            print subs['SMILES_new']
            if four in subs['SMILES_new']:
                updated = subs['SMILES_new'].replace(four, '+4')
                print updated
                raw_input("Update OK?")
                chemicals.update({'_id':id},{'$set':{'SMILES_new':updated}})
            if three in subs['SMILES_new']:
                updated = subs['SMILES_new'].replace(three, '+3')
                print updated
                raw_input("Update OK?")
                chemicals.update({'_id':id},{'$set':{'SMILES_new':updated}})
        except:
            pass
    
    for id in c:
        four = '++++'
        three = '+++'
        subs = chemicals.find_one({'_id':id})
        try:
            print subs['SMILES_new']
            if four in subs['SMILES_new']:
                updated = subs['SMILES_new'].replace(four, '+4')
                print updated
                raw_input("Update OK?")
                chemicals.update({'_id':id},{'$set':{'SMILES_new':updated}})
            if three in subs['SMILES_new']:
                updated = subs['SMILES_new'].replace(three, '+3')
                print updated
                raw_input("Update OK?")
                chemicals.update({'_id':id},{'$set':{'SMILES_new':updated}})
        except:
            pass
'''
print encounters
print chemicals.find_one({'_id':11342940})
print '{}'.format(total+encounters)
for chem in chemicals.find():
    if chem['SMILES'] == 'O.[HH].[NaH]':
        wrong+=1
    if chem['SMILES'] == '[Na+].[OH-]' or chem['SMILES'] == '[OH-].[Na+]':
        corr +=1
print wrong
print corr
'''