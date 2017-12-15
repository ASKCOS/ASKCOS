import global_config as gc
import numpy as np
import h5py
from keras.models import Sequential, load_model
from keras.layers import Dense, Activation, Flatten, Dropout, Input
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.utils import np_utils
from keras.datasets import mnist
from keras import optimizers
from matplotlib import pyplot as plt
from utilities.i_o.logging import MyLogger
from utilities.fingerprinting import get_condition_input_from_smiles, get_condition_input_from_instance, get_input_condition_as_smiles,get_reaction_as_smiles 
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import cPickle as pickle
import os
flowNN_loc = 'flowNN3'
all_data=True
FILE = True
SAVE = True
FullTest = False
training_bias = 2
if all_data:
    training_bias = 15

np.random.seed(123)
'''
Neural network for the recognition of flow compatible reaction conditions.
Input is the concatenation of a 256 bit fingerprint for the solvent and a 256 bit fingerprint for the reagents and catalyst
'''
def set_up_model_structure(layer_nodes = None, activation = 'softmax', input = None, loss = 'categorical_crossentropy'):
    if not layer_nodes:
        MyLogger.print_and_log('Cannot build model without information on layer node counts.', flowNN_loc, level =3)
    if not input:
        MyLogger.print_and_log('Cannot build model without information on input structure.', flowNN_loc, level =3)
    
    adam = optimizers.Adam(lr=0.01)
    model = Sequential()
    layer = 1
    for number in layer_nodes:
        if layer == 1:
            model.add(Dense(number, activation=activation, input_shape = (1,input)))
        else:
            model.add(Dense(number, activation=activation))
        #model.add(Dropout(0.1))
        layer +=1
    
    model.compile(loss=loss, optimizer=adam, metrics=['accuracy'])
    return model
    
def get_data(flow_database = None,chemicals = None, train_test_split = 0.85, write_to_file = False, read_from_file = False):
    
    input_train = []
    output_train = []
    input_test = []
    output_test = []
    input_train_doc = []
    input_test_doc = []
    weights = []
    if read_from_file:
        if all_data:
            if os.path.isfile(gc.FLOW_CONDITIONS3['data_loc']):
                with open(gc.FLOW_CONDITIONS3['data_loc'], 'rb') as file:
                    input_test = pickle.load(file)
                    input_train = pickle.load(file)
                    output_test = pickle.load(file)
                    output_train = pickle.load(file)
                    input_test_doc = pickle.load(file)
                    input_train_doc = pickle.load(file)
                    MyLogger.print_and_log('Flow condition data read from {}'.format(gc.FLOW_CONDITIONS3['data_loc']), flowNN_loc)
            else:
                MyLogger.print_and_log('Cannot load data from non-existent file. Exiting...', flowNN_loc, level = 3)
        else:
            if os.path.isfile(gc.FLOW_CONDITIONS3_50['data_loc']):
                with open(gc.FLOW_CONDITIONS3_50['data_loc'], 'rb') as file:
                    input_test = pickle.load(file)
                    input_train = pickle.load(file)
                    output_test = pickle.load(file)
                    output_train = pickle.load(file)
                    input_test_doc = pickle.load(file)
                    input_train_doc = pickle.load(file)
                    MyLogger.print_and_log('Flow condition data read from {}'.format(gc.FLOW_CONDITIONS3_50['data_loc']), flowNN_loc)
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
            
            fp = get_condition_input_from_instance(doc, chemicals, astwo = True, use_new = True)
            
            if(split < train_test_split):
                output_train.append(np.array([np.array([1.0,0.0]) if isFlow else np.array([0.0,1.0])]))
                input_train.append(fp)
                input_train_doc.append(doc)
            else:
                output_test.append(np.array([np.array([1.0,0.0] if isFlow else [0.0,1.0])]))
                input_test.append(fp)
                input_test_doc.append(doc)
                
    if write_to_file:
        if all_data:
            with open(gc.FLOW_CONDITIONS3['data_loc'], 'wb') as file:
                pickle.dump(input_test, file, gc.protocol)
                pickle.dump(input_train, file, gc.protocol)
                pickle.dump(output_test, file, gc.protocol)
                pickle.dump(output_train, file, gc.protocol)
                pickle.dump(input_test_doc, file, gc.protocol)
                pickle.dump(input_train_doc, file, gc.protocol)
                MyLogger.print_and_log('Flow condition data written to {}'.format(gc.FLOW_CONDITIONS3['data_loc']), flowNN_loc)
        else:
            with open(gc.FLOW_CONDITIONS3_50['data_loc'], 'wb') as file:
                pickle.dump(input_test, file, gc.protocol)
                pickle.dump(input_train, file, gc.protocol)
                pickle.dump(output_test, file, gc.protocol)
                pickle.dump(output_train, file, gc.protocol)
                pickle.dump(input_test_doc, file, gc.protocol)
                pickle.dump(input_train_doc, file, gc.protocol)
                MyLogger.print_and_log('Flow condition data written to {}'.format(gc.FLOW_CONDITIONS3_50['data_loc']), flowNN_loc)
    input_test =  np.array(input_test)
    input_train = np.array(input_train)
    input_test = input_test.astype('float32')
    input_train = input_train.astype('float32')
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

def train_model(model,X_train, Y_train, weights, batch_size, epochs):
    model.fit(X_train, Y_train, batch_size=batch_size, epochs=epochs, sample_weight = weights, verbose=1)

def test_model(model, input_test, output_test):
    return model.evaluate(input_test, output_test, verbose=0)

def model_to_file(model,full = True):
    
    if full:
        model.save(gc.FLOW_CONDITIONS3['model_loc'])
    else:
        model.save(gc.FLOW_CONDITIONS3_50['model_loc'])
def model_from_file(full = True):
    if full:
        return load_model(gc.FLOW_CONDITIONS3['model_loc'])
    else:
        return load_model(gc.FLOW_CONDITIONS3_50['model_loc'])
        
def tests():
    
    from pymongo import MongoClient
    
    MyLogger.initialize_logFile()
    client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', gc.MONGO['id'], connect = gc.MONGO['connect'])
    db2 = client[gc.FLOW_CONDITIONS2['database']]
    flow_database = None
    if all_data:
        flow_database = db2[gc.FLOW_CONDITIONS3['collection']]
    else:
        flow_database = db2[gc.FLOW_CONDITIONS3_50['collection']]
    instance_database = db2[gc.INSTANCES['collection']]
    chemicals = db2[gc.CHEMICALS['collection']]
    #(input_test, input_train, output_test, output_train,input_test_doc,input_train_doc, weights) = get_data(flow_database = flow_database,chemicals=chemicals, train_test_split=0.9, write_to_file=True)
    (input_test, input_train, output_test, output_train,input_test_doc,input_train_doc, weights) = get_data(train_test_split=0.9, read_from_file=True)
    #layer_nodes = [128,64,32,4, 2]
    layer_nodes = [256,128,32,2]
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
        
    
    print model.predict(np.array([input_test[165]]))
    print output_test[165]
    print input_test_doc[165]
    #instance = instance_database.find_one({'RX_ID':3124853})
    #trial = get_condition_input_from_instance(instance, chemicals, astwo = True, use_new = True)
    #print model.predict(trial)
    smiles = [('solv','O'),('cata','[Cs+].[F-]')]
    print 'Should be plausible: {}'.format(model.predict(np.array([get_condition_input_from_smiles(smiles)])))
    smiles = [('solv','CC#N'),('cata','[Cs+].[F-]')]
    print 'Should not be plausible: {}'.format(model.predict(np.array([get_condition_input_from_smiles(smiles)])))
    smiles = [('solv','O'),('cata',"")]
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
    
    from AUROC import get_AUROC
    test_true = []
    test_pred = []
    for i,input in enumerate(input_test):
        score = model.predict(np.array([input]))[0][0][0]
        test_true.append(1.0 if input_test_doc[i]['flow'] else 0.0)
        test_pred.append(score)
    print 'AUC test = {}'.format(get_AUROC(test_true,test_pred))
    train_true = []
    train_pred = []
    for i,input in enumerate(input_train):
        score = model.predict(np.array([input]))[0][0][0]
        train_true.append(1.0 if input_train_doc[i]['flow'] else 0.0)
        train_pred.append(score)
    print 'AUC train = {}'.format(get_AUROC(train_true, train_pred))
    
    if FullTest:
        FP = 0
        FN = 0
        UD = 0
        CP = 0
        CN = 0
        #Use full database for this test
        flow_database = db2[gc.FLOW_CONDITIONS3['collection']]
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
if __name__ == '__main__':
    tests()
    