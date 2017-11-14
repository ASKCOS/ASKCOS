import global_config as gc
import numpy as np
import h5py
from keras.models import Sequential, Model
from keras.layers import Dense, Activation, Flatten, Dropout, Input, merge, concatenate
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.utils import np_utils
from keras.datasets import mnist
from keras import optimizers
from matplotlib import pyplot as plt
from utilities.i_o.logging import MyLogger
from utilities.fingerprinting import get_condition_input_from_smiles, get_condition_input_from_instance, get_input_condition_as_smiles,get_reaction_as_smiles,get_reaction_input_from_instance 
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import cPickle as pickle
import os
flowNN_loc = 'flowNN4'
all_data=True
FILE = False
SAVE = True
FullTest = False
training_bias = 2
if all_data:
    training_bias = 15

np.random.seed(123)

def set_up_model_structure(layer_nodes = None, activation = 'softmax', input = None, loss = 'binary_crossentropy'):
    
    reaction = Input(shape = (1, 256), name = "reaction fingerprint")
    solvent = Input(shape = (1, 256), name = "solvent fingerprint")
    reagent = Input(shape = (1, 256), name = "reagent fingerprint")
    
    reaction_2 = Dense(128, activation = activation)(reaction)
    solvent_2 = Dense(128, activation = activation)(solvent)
    reagent_2 = Dense(128, activation = activation)(reagent)
    
    reaction_3 = Dense(64, activation = activation)(reaction_2)
    conditions = concatenate([solvent_2, reagent_2])
    conditions_2 = Dense(64, activation = activation)(conditions)
    
    decision = concatenate([reaction_3, conditions_2])
    decision_2 = Dense(64, activation = activation)(decision)
    result = Dense(1, activation = activation)(decision)
    
    model = Model(inputs = [reaction, solvent, reagent], outputs = [result])
    model.compile(loss = loss, optimizer = 'adam', metrics = ['accuracy'])
    
    return model
def get_data_workaround(reactions = None, chemicals = None, train_test_split = 0.85, write_to_file = False):
    input_docs = []
    if os.path.isfile(gc.FLOW_CONDITIONS['data_loc']):
        with open(gc.FLOW_CONDITIONS['data_loc'], 'rb') as file:
            input_test = pickle.load(file)
            input_train = pickle.load(file)
            output_test = pickle.load(file)
            output_train = pickle.load(file)
            input_test_doc = pickle.load(file)
            input_train_doc = pickle.load(file)
    input_docs.extend(input_test_doc)
    input_docs.extend(input_train_doc)   
    input_reaction_train = []
    input_solvent_train = []
    input_reagent_train = []
    output_train = []
    input_reaction_test = []
    input_solvent_test = []
    input_reagent_test = []
    output_test = []
    input_train_doc = []
    input_test_doc = []
    weights = [] 
    for doc in input_docs:
        split = np.random.random_sample()
        isFlow = doc['flow']
        
        fps = get_condition_input_from_instance(doc, chemicals, astwo = True, split = True, use_new = True)
        reac = get_reaction_input_from_instance(doc, reactions, chemicals)
        
        if(split < train_test_split):
            input_reaction_train.append(reac)
            input_solvent_train.append(fps[0])
            input_reagent_train.append(fps[1])
            output_train.append(np.array([np.array([1.0] if isFlow else [0.0])]))
        else:
            input_reaction_test.append(reac)
            input_solvent_test.append(fps[0])
            input_reagent_test.append(fps[1])
            output_test.append(np.array([np.array([1.0] if isFlow else [0.0])]))
            
    if write_to_file:
        if all_data:
            with open(gc.FLOW_CONDITIONS4['data_loc'], 'wb') as file:
                pickle.dump(input_reaction_test, file, gc.protocol)
                pickle.dump(input_reaction_train, file, gc.protocol)
                pickle.dump(input_solvent_test, file, gc.protocol)
                pickle.dump(input_solvent_train, file, gc.protocol)
                pickle.dump(input_reagent_test, file, gc.protocol)
                pickle.dump(input_reagent_train, file, gc.protocol)
                pickle.dump(output_test, file, gc.protocol)
                pickle.dump(output_train, file, gc.protocol)
                pickle.dump(input_test_doc, file, gc.protocol)
                pickle.dump(input_train_doc, file, gc.protocol)
                MyLogger.print_and_log('Flow condition data written to {}'.format(gc.FLOW_CONDITIONS4['data_loc']), flowNN_loc)
        else:
            with open(gc.FLOW_CONDITIONS4_50['data_loc'], 'wb') as file:
                pickle.dump(input_reaction_test, file, gc.protocol)
                pickle.dump(input_reaction_train, file, gc.protocol)
                pickle.dump(input_solvent_test, file, gc.protocol)
                pickle.dump(input_solvent_train, file, gc.protocol)
                pickle.dump(input_reagent_test, file, gc.protocol)
                pickle.dump(input_reagent_train, file, gc.protocol)
                pickle.dump(output_test, file, gc.protocol)
                pickle.dump(output_train, file, gc.protocol)
                pickle.dump(input_test_doc, file, gc.protocol)
                pickle.dump(input_train_doc, file, gc.protocol)
                MyLogger.print_and_log('Flow condition data written to {}'.format(gc.FLOW_CONDITIONS4_50['data_loc']), flowNN_loc)
    
    input_reaction_test =  np.array(input_reaction_test).astype('float32')
    input_reaction_train = np.array(input_reaction_train).astype('float32')
    input_solvent_test = np.array(input_solvent_test).astype('float32')
    input_solvent_train = np.array(input_solvent_train).astype('float32')
    input_reagent_test = np.array(input_reagent_test).astype('float32')
    input_reagent_train = np.array(input_reagent_train).astype('float32')
    output_test = np.array(output_test)
    output_train = np.array(output_train)
    output_test = output_test.astype('float32')
    output_train = output_train.astype('float32')
    for data in output_train:
        for set in data:
            if set[0] == 1:
                weights.append(training_bias)
            else:
                weights.append(1)
    return (input_test, input_train, output_test, output_train,input_test_doc,input_train_doc, np.array(weights))

def get_data(flow_database = None,reactions = None, chemicals = None, train_test_split = 0.85, write_to_file = False, read_from_file = False):
    
    input_reaction_train = []
    input_solvent_train = []
    input_reagent_train = []
    output_train = []
    input_reaction_test = []
    input_solvent_test = []
    input_reagent_test = []
    output_test = []
    input_train_doc = []
    input_test_doc = []
    weights = []
    if read_from_file:
        if all_data:
            if os.path.isfile(gc.FLOW_CONDITIONS4['data_loc']):
                with open(gc.FLOW_CONDITIONS4['data_loc'], 'rb') as file:
                    input_reaction_test = pickle.load(file)
                    input_reaction_train = pickle.load(file)
                    input_solvent_test = pickle.load(file)
                    input_solvent_train = pickle.load(file)
                    input_reagent_test = pickle.load(file)
                    input_reagent_train = pickle.load(file)
                    output_test = pickle.load(file)
                    output_train = pickle.load(file)
                    input_test_doc = pickle.load(file)
                    input_train_doc = pickle.load(file)
                    MyLogger.print_and_log('Flow condition data read from {}'.format(gc.FLOW_CONDITIONS4['data_loc']), flowNN_loc)
            else:
                MyLogger.print_and_log('Cannot load data from non-existent file. Exiting...', flowNN_loc, level = 3)
        else:
            if os.path.isfile(gc.FLOW_CONDITIONS4_50['data_loc']):
                with open(gc.FLOW_CONDITIONS4_50['data_loc'], 'rb') as file:
                    input_reaction_test = pickle.load(file)
                    input_reaction_train = pickle.load(file)
                    input_solvent_test = pickle.load(file)
                    input_solvent_train = pickle.load(file)
                    input_reagent_test = pickle.load(file)
                    input_reagent_train = pickle.load(file)
                    output_test = pickle.load(file)
                    output_train = pickle.load(file)
                    input_test_doc = pickle.load(file)
                    input_train_doc = pickle.load(file)
                    MyLogger.print_and_log('Flow condition data read from {}'.format(gc.FLOW_CONDITIONS4_50['data_loc']), flowNN_loc)
            else:
                MyLogger.print_and_log('Cannot load data from non-existent file. Exiting...', flowNN_loc, level = 3)
    else:
        if not flow_database:
            MyLogger.print_and_log('Cannot retrieve data without a flow database.', flowNN_loc, level =3)
        if not chemicals:
            MyLogger.print_and_log('Cannot retrieve data without chemicals database', flowNN_loc, level = 3)
        i = 0
        
        for doc in flow_database.find():
            '''
            if i>1000:
                break
            i+=1
            '''
            split = np.random.random_sample()
            isFlow = doc['flow']
            
            fps = get_condition_input_from_instance(doc, chemicals, astwo = True, split = True, use_new = True)
            reac = get_reaction_input_from_instance(doc, reactions, chemicals)
            
            if(split < train_test_split):
                input_reaction_train.append(reac)
                input_solvent_train.append(fps[0])
                input_reagent_train.append(fps[1])
                output_train.append(np.array([np.array([1.0] if isFlow else [0.0])]))
            else:
                input_reaction_test.append(reac)
                input_solvent_test.append(fps[0])
                input_reagent_test.append(fps[1])
                output_test.append(np.array([np.array([1.0] if isFlow else [0.0])]))
                
    if write_to_file:
        if all_data:
            with open(gc.FLOW_CONDITIONS4['data_loc'], 'wb') as file:
                pickle.dump(input_reaction_test, file, gc.protocol)
                pickle.dump(input_reaction_train, file, gc.protocol)
                pickle.dump(input_solvent_test, file, gc.protocol)
                pickle.dump(input_solvent_train, file, gc.protocol)
                pickle.dump(input_reagent_test, file, gc.protocol)
                pickle.dump(input_reagent_train, file, gc.protocol)
                pickle.dump(output_test, file, gc.protocol)
                pickle.dump(output_train, file, gc.protocol)
                pickle.dump(input_test_doc, file, gc.protocol)
                pickle.dump(input_train_doc, file, gc.protocol)
                MyLogger.print_and_log('Flow condition data written to {}'.format(gc.FLOW_CONDITIONS4['data_loc']), flowNN_loc)
        else:
            with open(gc.FLOW_CONDITIONS4_50['data_loc'], 'wb') as file:
                pickle.dump(input_reaction_test, file, gc.protocol)
                pickle.dump(input_reaction_train, file, gc.protocol)
                pickle.dump(input_solvent_test, file, gc.protocol)
                pickle.dump(input_solvent_train, file, gc.protocol)
                pickle.dump(input_reagent_test, file, gc.protocol)
                pickle.dump(input_reagent_train, file, gc.protocol)
                pickle.dump(output_test, file, gc.protocol)
                pickle.dump(output_train, file, gc.protocol)
                pickle.dump(input_test_doc, file, gc.protocol)
                pickle.dump(input_train_doc, file, gc.protocol)
                MyLogger.print_and_log('Flow condition data written to {}'.format(gc.FLOW_CONDITIONS4_50['data_loc']), flowNN_loc)
    
    input_reaction_test =  np.array(input_reaction_test).astype('float32')
    input_reaction_train = np.array(input_reaction_train).astype('float32')
    input_solvent_test = np.array(input_solvent_test).astype('float32')
    input_solvent_train = np.array(input_solvent_train).astype('float32')
    input_reagent_test = np.array(input_reagent_test).astype('float32')
    input_reagent_train = np.array(input_reagent_train).astype('float32')
    output_test = np.array(output_test)
    output_train = np.array(output_train)
    output_test = output_test.astype('float32')
    output_train = output_train.astype('float32')
    for data in output_train:
        for set in data:
            if set[0] == 1:
                weights.append(training_bias)
            else:
                weights.append(1)
    return (input_reaction_test, input_reaction_train, input_solvent_test,input_solvent_train, input_reagent_test, input_reagent_train, output_test, output_train,input_test_doc,input_train_doc, np.array(weights))

def train_model(model,X_train, Y_train, weights, batch_size, epochs):
    model.fit(X_train, Y_train, batch_size=batch_size, epochs=epochs, sample_weight = weights, verbose=1)

def test_model(model, input_test, output_test):
    return model.evaluate(input_test, output_test, verbose=0)

def model_to_file(model,full = True):
    
    if full:
        model.save(gc.FLOW_CONDITIONS4['model_loc'])
    else:
        model.save(gc.FLOW_CONDITIONS4_50['model_loc'])
def model_from_file(full = True):
    if full:
        return load_model(gc.FLOW_CONDITIONS4['model_loc'])
    else:
        return load_model(gc.FLOW_CONDITIONS4_50['model_loc'])
        
def tests():
    
    from pymongo import MongoClient
    
    MyLogger.initialize_logFile()
    client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    db2 = client[gc.FLOW_CONDITIONS4['database']]
    flow_database = None
    if all_data:
        flow_database = db2[gc.FLOW_CONDITIONS4['collection']]
    else:
        flow_database = db2[gc.FLOW_CONDITIONS4_50['collection']]
    instance_database = db2[gc.INSTANCES['collection']]
    chemicals = db2[gc.CHEMICALS['collection']]
    reactions = db2[gc.REACTIONS['collection']]
    (input_reaction_test, input_reaction_train, input_solvent_test,input_solvent_train, input_reagent_test, input_reagent_train, output_test, output_train,input_test_doc,input_train_doc, weights) = get_data(flow_database = flow_database,reactions = reactions,chemicals=chemicals, train_test_split=0.9, write_to_file=True)
    #(input_test, input_train, output_test, output_train,input_test_doc,input_train_doc, weights) = get_data_workaround(reactions = reactions,chemicals=chemicals, train_test_split=0.9, write_to_file=True)
    #(input_test, input_train, output_test, output_train,input_test_doc,input_train_doc, weights) = get_data(train_test_split=0.9, read_from_file=True)
    input_train = [input_reaction_train, input_solvent_train, input_reagent_train]
    input_test = [input_reaction_test, input_solvent_test, input_reagent_test]
    layer_nodes = [128,64,32,4, 2]
    #layer_nodes = [256,128,32,2]
    if FILE:
        model = model_from_file(full = all_data)
    else:
        model = set_up_model_structure(layer_nodes = layer_nodes, input = gc.fingerprint_bits*2)
        train_model(model, input_train, output_train, weights, 32, 10)
        print test_model(model, input_test, output_test)
    if SAVE:
        try:
            model_to_file(model, full = all_data)
        except Exception as e:
            print e
        
    
    #print model.predict(np.array([input_test[165]]))
    print output_test[165]
    print input_test_doc[165]
    #instance = instance_database.find_one({'RX_ID':3124853})
    #trial = get_condition_input_from_instance(instance, chemicals, astwo = True, use_new = True)
    #print model.predict(trial)
    '''
    smiles = ['O','[Cs+].[F-]']
    print 'Should be plausible: {}'.format(model.predict(np.array([get_condition_input_from_smiles(smiles)])))
    smiles = ['CC#N','[Cs+].[F-]']
    print 'Should not be plausible: {}'.format(model.predict(np.array([get_condition_input_from_smiles(smiles)])))
    smiles = ['O','']
    print 'Should be close to certain: {}'.format(model.predict(np.array([get_condition_input_from_smiles(smiles)])))
    instance = flow_database.find_one({'flow':True})
    id = instance['RX_ID']
    cont = True
    inp = None
    try:
        while cont:
            instance = flow_database.find_one({'$and':[{'flow':True},{'RX_ID':{'$gt':id}}]})
            id = instance['RX_ID']
            inp = get_condition_input_from_instance(instance, chemicals, astwo = True, use_new = True)
            if inp == None:
                cont= True
            else:
                cont =False
    except ValueError:
        pass
    print '{} @ T:{} should be true: {}'.format(id, instance['RXD_T'], model.predict(np.array([inp])))
    
    inp = None
    try:
        while cont:
            instance = flow_database.find_one({'$and':[{'flow':False},{'RX_ID':{'$gt':id}}]})
            id = instance['RX_ID']
            inp = get_condition_input_from_instance(instance, chemicals, astwo = True, use_new = True)
            if inp == None:
                cont= True
            else:
                cont =False
    except ValueError:
        pass
    print '{} @ T:{} should be false: {}'.format(id,instance['RXD_T'], model.predict(np.array([inp])))
    
    
    if FullTest:
        FP = 0
        FN = 0
        UD = 0
        CP = 0
        CN = 0
        #Use full database for this test
        flow_database = db2[gc.FLOW_CONDITIONS4['collection']]
        reaction_database = db2[gc.REACTIONS['collection']]
        for i,instance in enumerate(flow_database.find()):
            score = model.predict(np.array([get_condition_input_from_instance(instance,chemicals, astwo = True, use_new = True)]))
            #print 'Should be {}, is {}'.format(instance['flow'], score)
            if abs(score[0][0][0] - (1 if instance['flow'] else 0)) > 0.6666:
                if abs(score[0][0][0] - (1 if instance['flow'] else 0)) > 0.85:
                    print 'Reaction #{} with Reaxys ID {} should be {}, is {}'.format(i, instance['_id'], instance['flow'], score)
                    conditions = get_input_condition_as_smiles(instance, chemicals)
                    reaction_smiles = get_reaction_as_smiles(instance, reaction_database, chemicals)
                    print '\tSmiles of failed reaction: {}'.format(reaction_smiles)
                    print '\tSolvent: {}\tReagent: {}\tCatalyst: {}'.format(conditions[0],conditions[1],conditions[2])
                if instance['flow']:
                    FN += 1
                else:
                    FP += 1
            if abs(score[0][0][0] - (1 if instance['flow'] else 0)) < 0.6666 and abs(score[0][0][0] - (1 if instance['flow'] else 0)) > 0.3333:
                UD += 1
            
            if abs(score[0][0][0] - (1 if instance['flow'] else 0)) < 0.6666 and abs(score[0][0][0] - (1 if instance['flow'] else 0)) < 0.3333:
                if instance['flow']:
                    CP += 1
                else:
                    CN += 1
                
        print 'False positives: {}'.format(FP)
        print 'False negatives: {}'.format(FN)      
        print 'Correct positives: {}'.format(CP)
        print 'Correct negatives: {}'.format(CN)
        print 'Number of undecided outcomes: {}'.format(UD)
        '''
if __name__ == '__main__':
    tests()
    