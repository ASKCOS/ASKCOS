##############################################################################
#Settings:
all_data=True #Use the full data set or use a more balanced set

read_inputs_from_file = False #Read the preprocessed data from a file. Only use if no changes to the settings 
                              #have been made. Otherwise read from raw data
                              
read_raw_data_from_file = False #Read the raw data from a saved file or use MongoDB

write_raw_data_to_file = True #Write the extracted data to a local file (data is saved as smiles strings)

write_input_data_to_file = True #Write the extracted data to a local file (actual input data)

FILE = False #Read a trained and saved model from the file

SAVE = True #Save the trained model to the file

FullTest = True #Test again on all reactions in the database (including training)

print_text = False

#Data input settings
rxn_fingerprint_size = 1024 #Size of the reaction fingerprint
s_fingerprint_size = 256 #Size of the solvent fingerprint
r_fingerprint_size = 1024 #Size of the reagent fingerprint
rxn_compression = 1 #Compression of the reaction fingerprint. Make sure that rxn_fingerprint_size/rxn_compression is int.
training_bias = 2 #Bias when using the balanced data set
if all_data:
    training_bias = 7 #Bia when using the full database
    
#Model settings
activation_function = 'sigmoid' #specify activation function for the network
optimizer = 'adam' #specify the solver for the network
loss = 'binary_crossentropy' #specify the loss function for the network

#Training settings
batch_size = 32
epochs = 25

#Testing settings
lower_threshold = 1.0/3.0 #Threshold for considering a prediction accurate.
upper_threshold = 2.0/3.0 #Threshold for considering a prediction incorrect.
#Region between the thresholds is the undecided region

################################################################################

#Imports
import global_config as gc
import numpy as np
import time
import h5py
from keras.models import Sequential, Model, load_model
from keras.layers import Dense, Activation, Flatten, Dropout, Input, merge, concatenate
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.utils import np_utils
from keras.datasets import mnist
from keras import optimizers
from matplotlib import pyplot as plt
from pymongo import MongoClient
from utilities.io.logging import MyLogger
from utilities.fingerprinting import get_condition_input_from_smiles, get_condition_input_from_instance, get_input_condition_as_smiles,get_reaction_as_smiles,get_reaction_input_from_instance, get_reaction_input_from_smiles 
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import cPickle as pickle
import os
flowNN_loc = 'flowNN4'
np.random.seed(123)
if rxn_fingerprint_size%rxn_compression==0:
    rxn_compressed_size = int(rxn_fingerprint_size/rxn_compression)
else:
    MyLogger.print_and_log('Compression factor must be divisor of reaction fingerprint length. Exiting... ', flowNN_loc, level = 3)

def set_up_model_structure(activation = 'sigmoid', optimizer = 'adam', loss = 'binary_crossentropy'):
    
    reaction = Input(shape = (1, rxn_compressed_size), name = "reaction fingerprint")
    solvent = Input(shape = (1, s_fingerprint_size), name = "solvent fingerprint")
    reagent = Input(shape = (1, r_fingerprint_size), name = "reagent fingerprint")
    
    reaction_2 = Dense(128, activation = activation)(reaction)
    solvent_2 = Dense(128, activation = activation)(solvent)
    reagent_2 = Dense(128, activation = activation)(reagent)
    
    reaction_3 = Dense(64, activation = activation)(reaction_2)
    
    conditions = concatenate([solvent_2, reagent_2])
    conditions_2 = Dense(64, activation = activation)(conditions)
    
    decision = concatenate([reaction_3, conditions_2])
    decision_2 = Dense(64, activation = activation)(decision)
    decision_3 = Dense(32, activation = activation)(decision_2)
    result = Dense(1, activation = activation)(decision_3)
    
    model = Model(inputs = [reaction, solvent, reagent], outputs = [result])
    model.compile(loss = loss, optimizer = optimizer, metrics = ['accuracy'])
    
    return model

'''
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
'''

def get_data(flow_database = None,reactions = None, chemicals = None, train_test_split = 0.85, write_raw_to_file = False, 
             write_input_to_file = False, read_raw_from_file = False, read_inputs_from_file = False):
    
    smiles_train = []
    smiles_test = []
    input_reaction_train = []
    input_solvent_train = []
    input_reagent_train = []
    output_train = []
    input_reaction_test = []
    input_solvent_test = []
    input_reagent_test = []
    output_test = []
    weights = []
    ids_test = []
    ids_train = []
    smiles_data = []
    #############################################################################
    #If told to read the inputs directly from a file: load from file
    #############################################################################
    if read_inputs_from_file:
        if all_data:
            if os.path.isfile(gc.FLOW_CONDITIONS4['data_loc']):
                with open(gc.FLOW_CONDITIONS4['data_loc'], 'rb') as file:
                    data_set = pickle.load(file)
                    MyLogger.print_and_log('Flow condition data read from {}'.format(gc.FLOW_CONDITIONS4['data_loc']), flowNN_loc)
            else:
                MyLogger.print_and_log('Cannot load data from non-existent file. Exiting...', flowNN_loc, level = 3)
        else:
            if os.path.isfile(gc.FLOW_CONDITIONS4_50['data_loc']):
                with open(gc.FLOW_CONDITIONS4_50['data_loc'], 'rb') as file:
                    data_set = pickle.load(file)
                    MyLogger.print_and_log('Flow condition data read from {}'.format(gc.FLOW_CONDITIONS4_50['data_loc']), flowNN_loc)
            else:
                MyLogger.print_and_log('Cannot load data from non-existent file. Exiting...', flowNN_loc, level = 3)
    
    #############################################################################
    #Otherwise:
    #############################################################################
    else:
        
        #############################################################################
        #If told to read the raw, unprocessed smiles data from a file: load from file
        #############################################################################
        if read_raw_from_file:
            if all_data:
                if os.path.isfile(gc.FLOW_CONDITIONS4['raw_data_loc']):
                    with open(gc.FLOW_CONDITIONS4['raw_data_loc'], 'rb') as file:
                        smiles_data = pickle.load(file)
                        MyLogger.print_and_log('Raw flow condition data read from {}'.format(gc.FLOW_CONDITIONS4['raw_data_loc']), flowNN_loc)
                else:
                    MyLogger.print_and_log('Cannot load data from non-existent file. Exiting...', flowNN_loc, level = 3)
            else:
                if os.path.isfile(gc.FLOW_CONDITIONS4_50['raw_data_loc']):
                    with open(gc.FLOW_CONDITIONS4_50['raw_data_loc'], 'rb') as file:
                        smiles_data = pickle.load(file)
                        MyLogger.print_and_log('Raw flow condition data read from {}'.format(gc.FLOW_CONDITIONS4_50['raw_data_loc']), flowNN_loc)
                else:
                    MyLogger.print_and_log('Cannot load data from non-existent file. Exiting...', flowNN_loc, level = 3)
        
        #############################################################################
        #Otherwise: Read data from pymongo database
        #############################################################################
        else:
            if not flow_database:
                MyLogger.print_and_log('Cannot retrieve data without a flow database.', flowNN_loc, level =3)
            if not chemicals:
                MyLogger.print_and_log('Cannot retrieve data without chemicals database', flowNN_loc, level = 3)
            i = 0
            for doc in flow_database.find():
                if (i%1000 == 0):
                    MyLogger.print_and_log('Processed {} reactions.'.format(i), flowNN_loc)
                i+=1
                split = np.random.random_sample()
                isFlow = doc['flow']
                condition_smiles = get_input_condition_as_smiles(doc, chemicals, astwo = True, use_new = True)
                reaction_smiles = get_reaction_as_smiles(doc, reactions, chemicals)
                smiles_data.append({
                                    'isFlow': 1.0 if isFlow else 0.0, 
                                    'condition_smiles': condition_smiles, 
                                    'reaction_smiles': reaction_smiles, 
                                    '_id': doc['_id']
                                    })
        
        #############################################################################
        #If desired: write the raw data to a local file 
        #############################################################################
        if write_raw_to_file:
            if all_data:
                with open(gc.FLOW_CONDITIONS4['raw_data_loc'], 'wb') as file:
                    pickle.dump(smiles_data, file, gc.protocol)
                    MyLogger.print_and_log('Raw flow condition data written to {}'.format(gc.FLOW_CONDITIONS4['data_loc']), flowNN_loc)
            else:
                with open(gc.FLOW_CONDITIONS4_50['raw_data_loc'], 'wb') as file:
                    pickle.dump(smiles_data, file, gc.protocol)
                    MyLogger.print_and_log('Raw flow condition data written to {}'.format(gc.FLOW_CONDITIONS4_50['data_loc']), flowNN_loc)
        
        #############################################################################
        #Start preprocessing the raw data
        #############################################################################
        i = 0
        #Get the fingerprints
        for data in smiles_data:
            fps = get_condition_input_from_smiles(data['condition_smiles'], split = True, r_fp = r_fingerprint_size, s_fp = s_fingerprint_size)
            reac = get_reaction_input_from_smiles(data['reaction_smiles'], r_fp = rxn_fingerprint_size, c_f = rxn_compression)
    
            #Normalize (for compressed reaction fingerprints)
            max_val = max([np.max(reac), abs(np.min(reac))])
            reac = reac/float(1.0 if max_val == 0 else max_val)
            
            split = np.random.random_sample()
            #Make the training/testing split
            if(split < train_test_split):
                smiles_train.append(data['reaction_smiles'])
                input_reaction_train.append(reac)
                input_solvent_train.append(fps[0])
                input_reagent_train.append(fps[1])
                ids_train.append(data['_id'])
                output_train.append(np.array([np.array([data['isFlow']])]))
                if data['isFlow']:
                    weights.append(training_bias)
                else:
                    weights.append(1)
            else:
                smiles_test.append(data['reaction_smiles'])
                input_reaction_test.append(reac)
                input_solvent_test.append(fps[0])
                input_reagent_test.append(fps[1])
                ids_test.append(data['_id'])
                output_test.append(np.array([np.array([data['isFlow']])]))
                
            if i%1000 ==0:
                MyLogger.print_and_log('Input generated for {} reactions'.format(i),flowNN_loc)
            i+=1
        #Make sure the vectors are in the correct format
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
        data_set = {
                    'smiles_test':smiles_test,
                    'smiles_train':smiles_train,
                    'input_test_reaction':input_reaction_test,
                    'input_train_reaction': input_reaction_train,
                    'input_test_solvent': input_solvent_test,
                    'input_train_solvent': input_solvent_train,
                    'input_test_reagent': input_reagent_test,
                    'input_train_reagent': input_reagent_train,
                    'output_test': output_test,
                    'output_train': output_train,
                    'weights': np.array(weights),
                    'ids_test': ids_test,
                    'ids_train':ids_train
                    }
        
        #############################################################################
        #If desired, write the processed data to a file
        #############################################################################
        if write_input_to_file:
            if all_data:
                with open(gc.FLOW_CONDITIONS4['data_loc'], 'wb') as file:
                    pickle.dump(data_set, file, gc.protocol)
                    MyLogger.print_and_log('Flow condition data written to {}'.format(gc.FLOW_CONDITIONS4['data_loc']), flowNN_loc)
            else:
                with open(gc.FLOW_CONDITIONS4_50['data_loc'], 'wb') as file:
                    pickle.dump(data_set, file, gc.protocol)
                    MyLogger.print_and_log('Flow condition data written to {}'.format(gc.FLOW_CONDITIONS4_50['data_loc']), flowNN_loc)
                    
    return data_set

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

def full_test(data, model, no_write=False):
    FP = 0
    FN = 0
    UD = 0
    CP = 0
    CN = 0
    test_no = len(data['ids_test'])
    train_no = len(data['ids_train'])
    #Use full database for this test
    client = MongoClient(gc.MONGO['path'], gc.MONGO[ 'id'], connect=gc.MONGO['connect'])
    db2 = client[gc.FLOW_CONDITIONS4['database']]
    flow_database = db2[gc.FLOW_CONDITIONS4['collection']]
    reaction_database = db2[gc.REACTIONS['collection']]
    data_test = []
    data_train = []
    true_test = []
    true_train = []
    pred_test = []
    pred_train = []
    for i in range(test_no + train_no):
        reac = None
        solv = None
        reag = None
        id = None
        data_line_test = []
        data_line_train = []
        if i<test_no:
            reac = data['input_test_reaction'][i].reshape(1,1,rxn_compressed_size)
            solv = data['input_test_solvent'][i].reshape(1,1,s_fingerprint_size)
            reag = data['input_test_reagent'][i].reshape(1,1,r_fingerprint_size)
            flow = data['output_test'][i][0][0]
            true_test.append(flow)
            id = data['ids_test'][i]
            smiles = data['smiles_test'][i]
            data_line_test.append('{}'.format(id))
            data_line_test.append('{}'.format(flow))
        else:
            reac = data['input_train_reaction'][i-test_no].reshape(1,1,rxn_compressed_size)
            solv = data['input_train_solvent'][i-test_no].reshape(1,1,s_fingerprint_size)
            reag = data['input_train_reagent'][i-test_no].reshape(1,1,r_fingerprint_size)
            flow = data['output_train'][i-test_no][0][0]
            true_train.append(flow)
            id = data['ids_train'][i-test_no]
            smiles = data['smiles_train'][i-test_no]
            data_line_train.append('{}'.format(id))
            data_line_train.append('{}'.format(flow))
            
            
        score = model.predict([reac, solv, reag])[0][0][0] 
        if i<test_no:
            pred_test.append(score)
            data_line_test.append('{}'.format(score))
        else:
            pred_train.append(score)
            data_line_train.append('{}'.format(score))
        #Output format to export to excel:
        #Format:
        #Reaction ID reaxys /t Reaction smiles /t Reported compatibility /t Predicted Compatibility /t Undecided /t False Negative /t False Positive
        text_line = ['{}'.format(id), '{}'.format(smiles), '{}'.format(flow), '{}'.format(score), 'No', 'No', 'No']
        
        #Wrong prediction
        if abs(score - flow) > upper_threshold:
            if flow:
                FN += 1
                text_line[5] = 'Yes'
            else:
                FP += 1
                text_line[6] = 'Yes'
        
        #Undecided prediction 
        elif abs(score - flow) > lower_threshold:
            UD += 1
            text_line[4] = 'Yes'
        
        #Correct prediction
        else:
            if flow:
                CP += 1
            else:
                CN += 1
        if print_text and not no_write:
            MyLogger.print_and_log('\t'.join(text_line),flowNN_loc)
        if data_line_test != []:
            data_test.append(data_line_test)
        elif data_line_train != []:
            data_train.append(data_line_train)
        else:
            pass
    
    if not no_write:
        MyLogger.print_and_log('Reporting individual training results:',flowNN_loc)
        for data in data_train:
            MyLogger.print_and_log('\t'.join(data), flowNN_loc)
        MyLogger.print_and_log('Reporting individual test results:',flowNN_loc)
        for data in data_test:
            MyLogger.print_and_log('\t'.join(data), flowNN_loc)
    from AUROC import get_AUROC
    MyLogger.print_and_log('False positives: {}'.format(FP),flowNN_loc)
    MyLogger.print_and_log('False negatives: {}'.format(FN),flowNN_loc)
    MyLogger.print_and_log('Correct positives: {}'.format(CP),flowNN_loc)
    MyLogger.print_and_log('Correct negatives: {}'.format(CN),flowNN_loc)
    MyLogger.print_and_log('Number of undecided outcomes: {}'.format(UD),flowNN_loc)
    (auc_train, TPR_train, FPR_train) = get_AUROC(true_train, pred_train)
    (auc_test, TPR_test, FPR_test) = get_AUROC(true_test, pred_test)
    MyLogger.print_and_log('AUROC for training data: {}'.format(auc_train),flowNN_loc)
    MyLogger.print_and_log(FPR_train,flowNN_loc)
    MyLogger.print_and_log(TPR_train,flowNN_loc)
    MyLogger.print_and_log('AUROC for test data: {}'.format(auc_test),flowNN_loc)
    MyLogger.print_and_log(FPR_test,flowNN_loc)
    MyLogger.print_and_log(TPR_test,flowNN_loc)
    
    

def basic_tests():
    
    #Obtain data
    MyLogger.print_and_log('Starting to read data', flowNN_loc)
    if read_inputs_from_file:
        data = get_data(train_test_split=0.9, read_inputs_from_file=True)
    elif read_raw_data_from_file:
        data = get_data(train_test_split=0.9, read_raw_from_file=True, write_input_to_file = write_input_data_to_file, write_raw_to_file=write_raw_data_to_file)
    else:
        import makeit.global_config as gc
        client = MongoClient(gc.MONGO['path'], gc.MONGO[ 'id'], connect=gc.MONGO['connect'])
        db2 = client[gc.FLOW_CONDITIONS4['database']]
        flow_database = None
        MyLogger.print_and_log('Initialized databases', flowNN_loc)
        
        if all_data:
            flow_database = db2[gc.FLOW_CONDITIONS4['collection']]
        else:
            flow_database = db2[gc.FLOW_CONDITIONS4_50['collection']]
        instance_database = db2[gc.INSTANCES['collection']]
        chemicals = db2[gc.CHEMICALS['collection']]
        reactions = db2[gc.REACTIONS['collection']]
        data = get_data(flow_database = flow_database, reactions = reactions,chemicals=chemicals, 
                        train_test_split=0.9, read_raw_from_file=False, write_raw_to_file=write_raw_data_to_file, 
                        write_input_to_file=write_input_data_to_file)
    
    input_train = [data['input_train_reaction'], data['input_train_solvent'], data['input_train_reagent']]
    input_test = [data['input_test_reaction'], data['input_test_solvent'], data['input_test_reagent']]
    
    MyLogger.print_and_log('Data read, setting up model.', flowNN_loc)
    
    #Get a trained model: either by loading from saved file or by retraining.
    if FILE:
        model = model_from_file(full = all_data)
    else:
        model = set_up_model_structure(activation=activation_function, loss=loss, optimizer=optimizer)
        train_model(model, input_train, data['output_train'], data['weights'], batch_size, epochs)
        print test_model(model, input_test, data['output_test'])
    if SAVE:
        try:
            model_to_file(model, full = all_data)
        except Exception as e:
            print e
    
    MyLogger.print_and_log('Trained model loaded.', flowNN_loc)
    
    #Test on a specific reaction
    MyLogger.print_and_log('Prediction is: {}'.
                           format(model.predict(
                                    [data['input_test_reaction'][165].reshape(1,1,rxn_compressed_size), 
                                    data['input_test_solvent'][165].reshape(1,1,s_fingerprint_size), 
                                    data['input_test_reagent'][165].reshape(1,1,r_fingerprint_size)]
                                )),
                           flowNN_loc
                           )
    
    MyLogger.print_and_log('Prediction should be: {}'.format(data['output_test'][165]),flowNN_loc)
    MyLogger.print_and_log('For reaction: {}'.format(data['ids_test'][165]),flowNN_loc)
    
    return (data,model)
        
if __name__ == '__main__':
    MyLogger.initialize_logFile()
    (data,model) = basic_tests()
    # Run test on the full reactions database
    if FullTest:
        full_test(data,model, no_write=True)
    