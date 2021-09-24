import ROOT as r 
import numpy as np
import pickle,math,os
from keras.models import Sequential
from keras.layers import Dense, Dropout

import tensorflow as tf
from tensorflow.compat.v1.keras.backend import set_session
from tensorflow.keras.models import Model, load_model

from keras import optimizers


def getCompiledModelA():
    # optimal so far ( 0.980, 0.966)
    model = Sequential()
    model.add(Dense(10,input_dim=9, activation='relu'))
    model.add(Dense(7, activation='relu'))
    model.add(Dense(3, activation='softmax'))
    adam = optimizers.Adam(lr=1e-4) 
    model.compile(loss='categorical_crossentropy', optimizer=adam, metrics=['accuracy','categorical_crossentropy'])
    return model

def getCompiledModelB():
    model = Sequential()
    model.add(Dense(10,input_dim=9, activation='relu'))
    model.add(Dropout(0.4))
    model.add(Dense(5, activation='relu'))
    model.add(Dropout(0.4))
    model.add(Dense(5, activation='relu'))
    model.add(Dense(3, activation='softmax'))
    adam = optimizers.Adam(lr=1e-4) 
    model.compile(loss='categorical_crossentropy', optimizer=adam, metrics=['accuracy','categorical_crossentropy'])
    return model

def saveTFmodel(keras_model_fullpath):
    model = load_model(keras_model_fullpath)
    dirname = os.path.dirname(keras_model_fullpath)
    tfmodelname = os.path.basename(keras_model_fullpath).replace('.h5','_tf')
    model.save(dirname+'/'+tfmodelname)


if __name__ == "__main__":

    data = pickle.load( open('vars.pkl','rb'))
    sums = np.sum(data['train_y'],axis=0)
    print(sums)

    sig = sums[0]
    bkg = sums[1] + sums[2]

    class_weight = { 0 : float((sig+bkg)/sig),
                     1 : float((sig+bkg)/bkg),
                     2 : float((sig+bkg)/bkg)}
    print ('weights will be', class_weight)

    with tf.compat.v1.Session(config=tf.compat.v1.ConfigProto(
            intra_op_parallelism_threads=4,
            inter_op_parallelism_threads=4)) as sess:
        set_session(sess)
        
        model = getCompiledModelA()
        #model = getCompiledModelB()

        history = model.fit( data['train_x'], data['train_y'], epochs=30, batch_size=70, validation_data=(data['test_x'], data['test_y']), class_weight=class_weight)

        # keras model (H5)
        kerasmodelname = os.getcwd()+'/trained_model_A.h5'
        model.save(kerasmodelname)
        # tf model (PB)
        saveTFmodel(kerasmodelname)
        
        pickle_out = open('history_A.pkl','wb')
        pickle.dump( history.history, pickle_out)
        pickle_out.close()
        
