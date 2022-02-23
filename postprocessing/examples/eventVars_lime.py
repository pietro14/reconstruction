import ROOT, itertools

from framework.datamodel import Collection 
from framework.eventloop import Module

class ClusterVarsLime(Module):
    def __init__(self,weightsFile):
        self.vars = ["sc_regr_integral"]

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
            ret["sc_regr_integral"].append(c.integral)

        for V in self.vars:
            self.out.fillBranch(V,ret[V])

        return True

regressedEnergy = lambda : ClusterVarsLime("pippo.txt")

