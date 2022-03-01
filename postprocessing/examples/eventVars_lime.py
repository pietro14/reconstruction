import ROOT, itertools
import math, joblib
import numpy as np

from framework.datamodel import Collection 
from framework.eventloop import Module

class ClusterVarsLime(Module):
    def __init__(self,weightsFiles):
        self.vars = []
        self.models = {}
        for regr,wfile in weightsFiles.items():
            print("loading regression model from file: ",wfile)
            self.vars.append("sc_{r}_integral".format(r=regr))
            self.models[regr] = joblib.load(wfile)
        
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
            for regr,model in self.models.items():
                if c.integral>0 and c.nhits>150 and c.rms>10 and c.xmean>250:
                    inputs = np.array([c.xmean,c.ymean,c.tgaussamp/c.tgausssigma,c.lgaussamp/c.lgausssigma])
                    y_pred = model.predict(inputs.reshape(1,-1))
                    #print ("X = ", inputs)
                    #print ("y = ",y_pred)
                    #print ("====")
                    ret["sc_{m}_integral".format(m=regr)].append(c.integral/y_pred)
                else:
                    ret["sc_{m}_integral".format(m=regr)].append(c.integral)
                    
        for V in self.vars:
            self.out.fillBranch(V,ret[V])

        return True

regressedEnergy = lambda : ClusterVarsLime({"regr"     : "../data/gbrLikelihood_440V_mse.sav",
                                            "qregr"    : "../data/gbrLikelihood_440V_q0.50.sav",
                                            "qregr_up" : "../data/gbrLikelihood_440V_q0.05.sav",
                                            "qregr_dn" : "../data/gbrLikelihood_440V_q0.95.sav"})

