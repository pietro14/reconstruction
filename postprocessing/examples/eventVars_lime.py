import ROOT, itertools
import math, joblib
import numpy as np

from framework.datamodel import Collection 
from framework.eventloop import Module

class ClusterVarsLime(Module):
    def __init__(self,weightsFile,scale):
        self.vars = ["sc_regr_integral"]
        self.model = joblib.load(weightsFile)
        self.scale = scale
        
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
            if c.integral/self.scale < 1.6 and c.nhits>150 and c.rms>10 and c.xmean>250:
                inputs = np.array([c.longrms,
                                   c.latrms,
                                   c.lfullrms,
                                   c.tfullrms,
                                   c.xmean,
                                   c.ymean,
                                   c.tgausssigma,
                                   c.lgausssigma,
                                   c.width,
                                   c.length,
                                   c.tgaussmean/c.width,
                                   c.lgaussmean/c.length,
                                   c.tgaussamp/c.integral,
                                   c.lgaussamp/c.integral])
                y_pred = self.model.predict(inputs.reshape(1,-1))
                # print ("X = ", inputs)
                # print ("y = ",y_pred)
                # print ("====")
                ret["sc_regr_integral"].append(y_pred * self.scale)
            else:
                ret["sc_regr_integral"].append(c.integral)

        for V in self.vars:
            self.out.fillBranch(V,ret[V])

        return True

regressedEnergy = lambda : ClusterVarsLime("../data/gbrLikelihood_440V_mse.sav",scale=8000)

