import global_config as gc
import numpy as np
import math
from keras.models import Sequential, Model
from keras.layers import Dense, Activation, Flatten, Dropout, Input, merge, concatenate
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.utils import np_utils
from keras.datasets import mnist
from keras import optimizers
from matplotlib import pyplot as plt
from utilities.io.logging import MyLogger
from utilities.fingerprinting import get_reaction_input_from_smiles, get_condition_input_from_smiles, get_reaction_input_from_instance, get_condition_input_from_instance
from utilities.strings import string_or_range_to_float
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import cPickle as pickle
from pymongo import MongoClient
import os
yieldNN_loc = 'yieldNN'

def get_data(yield_database = None,chemicals = None, reactions = None, train_test_split = 0.85, write_to_file = False, read_from_file = False):
    
    input_RP_train = []
    input_CO_train = []
    input_TE_train = []
    input_TI_train = []
    output_train = []
    input_RP_test = []
    input_CO_test = []
    input_TE_test = []
    input_TI_test = []
    output_test = []
    input_train_doc = []
    input_test_doc = []
    if read_from_file:
        if os.path.isfile(gc.YIELDS['data_loc']):
            with open(gc.YIELDS['data_loc'], 'rb') as file:
                input_RP_train = pickle.load(file)
                input_CO_train = pickle.load(file)
                input_TE_train = pickle.load(file)
                input_TI_train = pickle.load(file)
                output_train = pickle.load(file)
                input_RP_test = pickle.load(file)
                input_CO_test = pickle.load(file)
                input_TE_test = pickle.load(file)
                input_TI_test = pickle.load(file)
                output_test = pickle.load(file)
                input_train_doc = pickle.load(file)
                input_test_doc = pickle.load(file)
                MyLogger.print_and_log('Flow condition data read from {}'.format(gc.YIELDS['data_loc']), yieldNN_loc)
        else:
            MyLogger.print_and_log('Cannot load data from non-existent file. Exiting...', yieldNN_loc, level = 3)
    else:
        if not yield_database:
            MyLogger.print_and_log('Cannot retrieve data without a yield database.', yieldNN_loc, level =3)
        if not chemicals:
            MyLogger.print_and_log('Cannot retrieve data without chemicals database', yieldNN_loc, level = 3)
        if not reactions:
            MyLogger.print_and_log('Cannot retrieve data without reactions database', yieldNN_loc, level = 3)
        
        for doc in yield_database.find():
            arr = np.array([0])
            split = np.random.random_sample()
            try:
                fp_C = get_condition_input_from_instance(doc, chemicals, use_new = True)
                fp_R = get_reaction_input_from_instance(doc, reactions, chemicals)
            except Exception as e:
                print e
                continue
            if not (type(fp_C) is type(arr)):
                if fp_C == None:
                    MyLogger.print_and_log('Unparseable conditions encountered. Skipping', yieldNN_loc)
                    continue
            if not(type(fp_R) is type(arr)):
                if fp_R == None:
                    MyLogger.print_and_log('Unparseable reaction encountered. Skipping.', yieldNN_loc)
                    continue
            T = doc['RXD_T']
            T0 = string_or_range_to_float(T)
            T = string_or_range_to_float(T)+273.15 #Use Kelvin!
            t = doc['RXD_TIM']
            t0 = string_or_range_to_float(t)
            t = string_or_range_to_float(t)*3600.0 #Time in seconds!
            y = doc['RXD_NYD']
            y0 = y
            y = y/100.0
            if y0 == -1:
                MyLogger.print_and_log('No yield available', yieldNN_loc)
                continue
            if t0 == -1:
                MyLogger.print_and_log('No time available', yieldNN_loc)
                continue
            if T0 == -1:
                MyLogger.print_and_log('No temperature available', yieldNN_loc)
                continue
            if(split < train_test_split):
                output_train.append(np.array([np.array([y])]))
                input_RP_train.append(fp_R)
                input_CO_train.append(fp_C)
                input_TE_train.append(np.array([np.array([T])]))
                input_TI_train.append(np.array([np.array([t])]))
                input_train_doc.append(doc)
            else:
                output_test.append(np.array([np.array([y])]))
                input_RP_test.append(fp_R)
                input_CO_test.append(fp_C)
                input_TE_test.append(np.array([np.array([T])]))
                input_TI_test.append(np.array([np.array([t])]))
                input_test_doc.append(doc)
            
    if write_to_file:
        with open(gc.YIELDS['data_loc'], 'wb') as file:
            pickle.dump(input_RP_train, file, gc.protocol)
            pickle.dump(input_CO_train, file, gc.protocol)
            pickle.dump(input_TE_train, file, gc.protocol)
            pickle.dump(input_TI_train, file, gc.protocol)
            pickle.dump(output_train, file, gc.protocol)
            pickle.dump(input_RP_test, file, gc.protocol)
            pickle.dump(input_CO_test, file, gc.protocol)
            pickle.dump(input_TE_test, file, gc.protocol)
            pickle.dump(input_TI_test, file, gc.protocol)
            pickle.dump(output_test, file, gc.protocol)
            pickle.dump(input_train_doc, file, gc.protocol)
            pickle.dump(input_test_doc, file, gc.protocol)
            MyLogger.print_and_log('Yield estimation data written to {}'.format(gc.YIELDS['data_loc']), yieldNN_loc)
    input_RP_train = np.array(input_RP_train).astype('float32')
    input_CO_train = np.array(input_CO_train).astype('float32')
    input_TE_train = np.array(input_TE_train).astype('float32')
    input_TI_train = np.array(input_TI_train).astype('float32')
    output_train = np.array(output_train).astype('float32')
    input_RP_test = np.array(input_RP_test).astype('float32')
    input_CO_test = np.array(input_CO_test).astype('float32')
    input_TE_test = np.array(input_TE_test).astype('float32')
    input_TI_test = np.array(input_TI_test).astype('float32')
    output_test = np.array(output_test).astype('float32')

    return (input_RP_test, input_RP_train,input_CO_test, input_CO_train,input_TE_test, input_TE_train,input_TI_test, input_TI_train, output_test, output_train,input_test_doc,input_train_doc)

def set_up_model_structure():
    
    reactprod = Input(shape = (1, 1024), name = "Reaction fingerprint")
    conditions = Input(shape = (1, 3*256), name = "Condition fingerprints")
    temperature = Input(shape = (1,1), name = "Temperature")
    time = Input(shape = (1,1), name = "Time")
    template_1 = Dense(256, activation = 'tanh')(reactprod)
    template_2 = Dense(64, activation = 'tanh')(template_1)
    condition_1 = Dense(256, activation = 'tanh')(conditions)
    condition_2 = Dense(256, activation = 'tanh')(condition_1)
    template_condition = concatenate([template_2, condition_2])
    tc_1 = Dense(256, activation = 'tanh')(template_condition)
    tc_2 = Dense(128, activation = 'tanh')(tc_1)
    
    #final_condition = concatenate([tc_2, temperature, time])
    
    #fc_1 = Dense(8, activation = 'tanh')(final_condition)
    fc_2 = Dense(5, activation = 'tanh')(tc_2)
    final_condition = concatenate([fc_2, temperature, time])
    
    o_1 = Dense(32, activation = 'tanh')(final_condition)
    o_2 = Dense(32, activation = 'tanh')(o_1)
    yi = Dense(1, activation = 'linear')(o_1)
    #yi = 1.0-math.exp(fc_2[0]*temperature[0][0]**fc_2[1]*math.exp(-fc_2[2]/8.314/temperature[0][0]*time[0][0]))
    #yi = Dense(1, activation = 'linear')(fc_2)
    model = Model(inputs = [reactprod, conditions, temperature, time], outputs = [yi])
    model.compile(loss = 'mean_squared_error', optimizer = 'adam', metrics = ['accuracy'])
    return model

def train_model(model, input_train, output_train, batch_size, epochs):
    model.fit(input_train, output_train, batch_size=batch_size, epochs=epochs, verbose=1)
    
def test_model(model, input_test, output_test):
    return model.evaluate(input_test, output_test, verbose=0)
    
if __name__ == "__main__":
    MyLogger.initialize_logFile()
    model = set_up_model_structure()
    client = MongoClient(gc.MONGO['path_yield_flow'], gc.MONGO['id'], connect = gc.MONGO['connect'])
    yield_database = client[gc.YIELDS['database']][gc.YIELDS['collection']]
    chemicals = client[gc.CHEMICALS['database']][gc.CHEMICALS['collection']]
    reactions = client[gc.REACTIONS['database']][gc.REACTIONS['collection']]
    (input_RP_test, input_RP_train,input_CO_test, input_CO_train,input_TE_test, input_TE_train,input_TI_test, input_TI_train, output_test, output_train,input_test_doc,input_train_doc) = get_data(yield_database, chemicals, reactions, read_from_file=True, write_to_file=False)
    reaction = "CCCC(=O)O.CCCO>>CCCC(=O)OCCC.O"
    conditions = [('solv',"OCCCC"),('reag',"[Na+].[Cl-]"),('cata',"")]
    inpa = np.array([get_reaction_input_from_smiles(reaction)])
    inpb = np.array([get_condition_input_from_smiles(conditions)])
    inpc = np.array([np.array([np.array([345])])])
    inpd = np.array([np.array([np.array([10])])])
    print model.predict([inpa,inpb,inpc,inpd])
    print model.evaluate([input_RP_test,input_CO_test,input_TE_test,input_TI_test],output_test)
    model.fit([input_RP_train, input_CO_train, input_TE_train, input_TI_train], output_train, batch_size = 64, epochs = 25, verbose = 1)
    
    print model.predict([inpa,inpb,inpc,inpd])
    for i,a in enumerate(input_RP_test):
        score = model.predict([np.array([input_RP_test[i]]), np.array([input_CO_test[i]]), np.array([input_TE_test[i]]),np.array([input_TI_test[i]])])
        print '{}\t{}\t{}\t{}'.format(score,output_test[i],input_TI_test[i][0][0],input_TE_test[i][0][0])