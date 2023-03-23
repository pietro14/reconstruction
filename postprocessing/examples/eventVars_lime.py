import ROOT, itertools
import math, joblib
import numpy as np

from framework.datamodel import Collection 
from framework.eventloop import Module

class ClusterVarsLime(Module):
    def __init__(self,weightsFiles,verbose=False):
        self.vars = []
        self.models = {}
        for regr,wfile in weightsFiles.items():
            print("loading regression model from file: ",wfile)
            self.vars.append("sc_{r}".format(r=regr))
            self.models[regr] = joblib.load(wfile)
        self.verbose=verbose
        
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
                y_pred = c.integral if 'integral' in regr else -1
                if 4000<c.integral<14000 and 8<c.rms<16 and c.tgausssigma*0.152>0.6 and np.hypot(c.xmean-2304./2.,c.ymean-2304./2.)<900 and c.length*0.152<50: # training phase space
                    if 'integral' in regr:
                        inputs = np.array([c.integral,c.xmean,c.ymean,c.rms,c.width,c.integral/c.nhits,c.tgausssigma])
                    else:
                        inputs = np.array([c.xmean,c.ymean,c.width,c.integral/c.nhits,c.tgausssigma,c.lgausssigma,c.tfullrms,c.lfullrms,c.rms])
                    y_pred = model.predict(inputs.reshape(1,-1))
                    if self.verbose:
                        print ("REGR = ",regr)
                        print ("X = ", inputs)
                        if 'integral' in regr:
                            print ("(yraw,ypred: ratio pred/raw) = (%f, %f: %f) " % (c.integral,y_pred,y_pred/c.integral))
                        else:
                            print ("ypred = ",y_pred)
                        print ("====")
                ret["sc_{m}".format(m=regr)].append(y_pred)
                    
        for V in self.vars:
            self.out.fillBranch(V,ret[V])

        return True

regressedEnergy = lambda : ClusterVarsLime({"regr_integral"     : "../data/gbrLikelihood_440V_mse.sav",
                                            "qregr_integral"    : "../data/gbrLikelihood_440V_q0.50.sav",
                                            "qregr_integral_up" : "../data/gbrLikelihood_440V_q0.05.sav",
                                            "qregr_integral_dn" : "../data/gbrLikelihood_440V_q0.95.sav"}, verbose=False)

regressedZ = lambda : ClusterVarsLime({"z" : "../data/gbrLikelihood_Z_440V_mse.sav"}, verbose=False)
