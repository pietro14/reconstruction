import ROOT, itertools
import math, joblib
import numpy as np

from framework.datamodel import Collection 
from framework.eventloop import Module

class ClusterVarsLime(Module):
    def __init__(self,weightsFile):
        self.vars = ["sc_regr_integral"]
        self.model = joblib.load(weightsFile)
        
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for v in self.vars:
            self.out.branch(v, "F", lenVar="nSc")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def analyze(self,event):
        # prepare output
        ret = {}
        for V in self.vars: ret[V] = []

        clusters = Collection(event,"sc","nSc")
        for c in clusters:
            if c.integral>0 and c.nhits>150 and c.rms>10 and c.xmean>250:
                inputs = np.array([c.xmean,c.ymean,c.tgaussamp/c.tgausssigma,c.lgaussamp/c.lgausssigma])
                y_pred = self.model.predict(inputs.reshape(1,-1))
                #print ("X = ", inputs)
                #print ("y = ",y_pred)
                #print ("====")
                ret["sc_regr_integral"].append(c.integral/y_pred)
            else:
                ret["sc_regr_integral"].append(c.integral)

        for V in self.vars:
            self.out.fillBranch(V,ret[V])

        return True

regressedEnergy = lambda : ClusterVarsLime("../data/gbrLikelihood_440V_mse.sav")

