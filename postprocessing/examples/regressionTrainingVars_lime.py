import ROOT, itertools
import math, joblib
import numpy as np

from framework.datamodel import Collection 
from framework.eventloop import Module

class RegressionTrainingVarsLime(Module):
    def __init__(self):
        self.vars = ["sc_trueint"]
        #                 (runmin,runmax) : [(peak value, integral_min,  integral_max) ... repeated for 1 or 2 peaks]
        self.mapenergy = {(5801,5805) : [(9.1e3,  6e3,  12e3)],                        # Cu  8 keV
                          (5806,5810) : [(9.1e3,  6e3,   9e3), (15e3,   9e3,  18e3)],  # Rb 14 keV
                          (5811,5820) : [(9.1e3,  6e3,  12e3), (19e3,  12e3,  25e3)],  # Mo 18 keV
                          (5821,5830) : [(9.1e3,  6e3,  12e3), (24e3,  12e3,  34e3)],  # Ag 23 keV
                          (5831,5845) : [(9.1e3,  6e3,  12e3), (49e3,  20e3,  60e3)],  # Ba 49 keV
                          (5867,5911) : [(  8e3,  2e3,  12e3)]                         # Fe  6 keV  (multi-Z)
                          }
        self.sigma0 = 0.1  
        
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

        peaksInRun = []
        for runs,energies in self.mapenergy.items():
            if runs[0] <= event.run <= runs[1]:
                peaksInRun = energies
                break
            
        clusters = Collection(event,"sc","nSc")
        for c in clusters:
            etrue = -1
            for p in peaksInRun:
                if p[1] < c.integral < p[2]:
                    etrue = p[0] * np.random.normal(1,self.sigma0)
                    break
            ret["sc_trueint"].append(etrue)
        for V in self.vars:
            self.out.fillBranch(V,ret[V])

        return True

trueEnergy = lambda : RegressionTrainingVarsLime()

