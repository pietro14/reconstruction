import numpy as np
import pickle
import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt

from keras.models import load_model
from sklearn.metrics import roc_curve
from sklearn.metrics import auc

if __name__ == "__main__":

    data = pickle.load( open('vars.pkl','rb'))

    # fig, ax = plt.subplots()
    # plt.hist(data['train_x'])
    # fig.savefig('hist.png')

    #model = load_model('trained_model_A.h5')
    model = load_model('trained_model_A_BEST_DATA_31-03-2021.h5')

    prediction = model.predict(data['test_x'])
    x = data['test_x']
    y = np.argmax(data['test_y'], axis=1)

    classifier = np.sum( prediction[:,[0,1]], axis=1)

    fpr_keras, tpr_keras, thresholds_keras = roc_curve(y<=1, classifier)
    auc_keras = auc(fpr_keras, tpr_keras)

    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label='Keras (area = {:.3f})'.format(auc_keras))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')

    plt.show()
    plt.savefig('roc.png')


    for var in range(9):
        together = np.dstack( (x[:,var], y ) )[0]
        
        class1 =  (together[together[:,1] == 0]) [:,0]
        class2 =  (together[together[:,1] == 1]) [:,0]
        class3 =  (together[together[:,1] == 2]) [:,0]

        bins = 20
        plt.clf()
        plt.hist(class1, bins, alpha=0.5,density=True, label='nR')
        plt.hist(class2, bins, alpha=0.5,density=True, label='eR')
        plt.hist(class3, bins, alpha=0.5,density=True, label='other')
        plt.legend(loc='upper right')
        plt.show()
    
        plt.savefig('input_%d.png'%var)


    for node in range(3):
        together = np.dstack( (prediction[:,node], y ) )[0]
        
        class1 =  (together[together[:,1] == 0]) [:,0]
        class2 =  (together[together[:,1] == 1]) [:,0]
        class3 =  (together[together[:,1] == 2]) [:,0]
        
        bins = 20
        plt.clf()
        plt.hist(class1, bins, alpha=0.5,density=True, label='nR')
        plt.hist(class2, bins, alpha=0.5,density=True, label='eR')
        plt.hist(class3, bins, alpha=0.5,density=True, label='other')
        plt.legend(loc='upper right')
        plt.show()

        plt.savefig('output_%d.png'%node)

    together = np.dstack( (classifier,y))[0]
    class1 =  (together[together[:,1] == 0]) [:,0]
    class2 =  (together[together[:,1] == 1]) [:,0]
    class3 =  (together[together[:,1] == 2]) [:,0]
    bins = 20
    plt.clf()
    plt.hist(class1, bins, alpha=0.5,density=True, label='nR')
    plt.hist(class2, bins, alpha=0.5,density=True, label='eR')
    plt.hist(class3, bins, alpha=0.5,density=True, label='other')
    plt.legend(loc='upper right')
    plt.show()
    
    plt.savefig('output_combined.png')

