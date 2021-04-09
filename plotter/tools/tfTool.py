import os, ROOT
import numpy as np
from keras.models import load_model


class TFTool:
    def __init__(self, name, h5, vars, classes, varorder):
        self.name = name
        self.h5   = h5
        self.vars = vars
        self.classes = classes
        self.varorder= varorder
        self.debug   = False

        # set tensorflow interface
        variables_ = ROOT.vector('string')()
        classes_   = ROOT.vector('string')()
        for var in self.varorder: variables_.push_back( var ) 
        for cla in self.classes : classes_.push_back(cla)
        self.model = load_model(h5)
        self.outbranches = [ '%s_%s'%(x,self.name) for x in self.classes]
        

    def __call__(self, ev, isc):
        inp = []
        ret = {}
        for key in self.varorder:
            var = self.vars[key]
            if self.debug: print (key, var(ev,isc))
            inp.append(var(ev,isc))
        inp = np.asarray(inp)
        inp = np.array([ tuple(inp) ])
        # https://www.tensorflow.org/api_docs/python/tf/keras/Model#predict
        #res = self.model.predict(inp)[0] # this would be faster when passing an array of arrays (e.g. for all the events together)
        res = self.model(inp)[0] # this is faster than predict(inp) in our case hwre we pass just 1 array of variables per event
        for i,cla in enumerate(self.classes):
            if self.debug: print (cla, res[i])
            ret['%s_%s'%(self.name,cla)] = res[i]
        
        return ret

            
