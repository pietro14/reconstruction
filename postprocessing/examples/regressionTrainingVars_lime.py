import ROOT, itertools
import math, joblib
import numpy as np
import pandas as pd

from framework.datamodel import Collection 
from framework.eventloop import Module

def getly(df,hv,z):
    return np.mean(df[(df['vgem']==hv)&(df['z']==z)]['ly'].values)

class RegressionTrainingVarsLime(Module):
    def __init__(self):
        self.vars = ["sc_trueint","sc_truez","sc_hv"]
        self.runmap = {
            9364 :[ 440 , 5],
            9378 :[ 440 , 15],
            9388 :[ 440 , 15],
            9390 :[ 440 , 25],
            9402 :[ 440 , 25],
            9412 :[ 440 , 25],
            9422 :[ 440 , 25],
            9432 :[ 440 , 25],
            9442 :[ 440 , 25],
            9735 :[ 440 , 48],
            9745 :[ 440 , 36],

            9365 :[ 430 , 5],
            9379 :[ 430 , 15],
            9391 :[ 430 , 25],
            9403 :[ 430 , 25],
            9413 :[ 430 , 25],
            9423 :[ 430 , 25],
            9433 :[ 430 , 25],
            9443 :[ 430 , 25],
            9736 :[ 430 , 48],
            9746 :[ 430 , 36],

            9366 :[ 420 , 5],
            9380 :[ 420 , 15],
            9392 :[ 420 , 25],
            9404 :[ 420 , 25],
            9414 :[ 420 , 25],
            9424 :[ 420 , 25],
            9434 :[ 420 , 25],
            9444 :[ 420 , 25],
            9737 :[ 420 , 48],
            9747 :[ 420 , 36],

            9367 :[ 410 , 5],
            9381 :[ 410 , 15],
            9393 :[ 410 , 25],
            9405 :[ 410 , 25],
            9415 :[ 410 , 25],
            9425 :[ 410 , 25],
            9435 :[ 410 , 25],
            9445 :[ 410 , 25],
            9738 :[ 410 , 48],
            9748 :[ 410 , 36],

            9368 :[ 400 , 5],
            9382 :[ 400 , 15],
            9394 :[ 400 , 25],
            9406 :[ 400 , 25],
            9416 :[ 400 , 25],
            9426 :[ 400 , 25],
            9436 :[ 400 , 25],
            9446 :[ 400 , 25],
            9739 :[ 400 , 48],
            9749 :[ 400 , 36],

            9369 :[ 390 , 5],
            9383 :[ 390 , 15],
            9395 :[ 390 , 25],
            9407 :[ 390 , 25],
            9417 :[ 390 , 25],
            9427 :[ 390 , 25],
            9437 :[ 390 , 25],
            9740 :[ 390 , 48],
            9750 :[ 390 , 36],
            
            9370 :[ 380 , 5],
            9384 :[ 380 , 15],
            9396 :[ 380 , 25],
            9408 :[ 380 , 25],
            9418 :[ 380 , 25],
            9428 :[ 380 , 25],
            9438 :[ 380 , 25],
            9741 :[ 380 , 48],
            9751 :[ 380 , 36],

            9371 :[ 370 , 5],
            9385 :[ 370 , 15],
            9397 :[ 370 , 25],
            9409 :[ 370 , 25],
            9419 :[ 370 , 25],
            9429 :[ 370 , 25],
            9439 :[ 370 , 25],
            9742 :[ 370 , 48],
            9752 :[ 370 , 36],

            9372 :[ 360 , 5],
            9386 :[ 360 , 15],
            9398 :[ 360 , 25],
            9410 :[ 360 , 25],
            9420 :[ 360 , 25],
            9430 :[ 360 , 25],
            9440 :[ 360 , 25],
            9743 :[ 360 , 48],
            9753 :[ 360 , 36],
        }

        df = pd.read_pickle("data/fitm-zhvscans-ext.pkl")
        # the "target" variable is the value without saturation, i.e. the value at max z
        self.energy_range_map = {360: (1500,5000,getly(df,360,48)),
                                 370: (2000,6000,getly(df,370,48)),
                                 380: (2000,7000,getly(df,380,48)),
                                 390: (2500,8000,getly(df,390,48)),
                                 400: (2500,10000,getly(df,400,48)),
                                 410: (2500,11000,getly(df,410,48)),
                                 420: (3000,12000,getly(df,420,48)),
                                 430: (4000,14000,getly(df,430,48)),
                                 440: (4000,17000,getly(df,440,48)),
                                 }
            
        self.sigma0 = 0.001 # relative to energy
        self.sigmaz0 = 0.5 # absolute in Z (cm)
        
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

        if event.run not in self.runmap:
            clusters = Collection(event,"sc","nSc")
            for c in clusters:
                ret["sc_trueint"].append(-1)
                ret["sc_truez"].append(-1)
                ret["sc_hv"].append(-1)
        else:
            hv = self.runmap[event.run][0]
            z = self.runmap[event.run][1]
            clusters = Collection(event,"sc","nSc")
            for c in clusters:
                etrue = -1
                ztrue = -1
                emin,emax,etrue0 = self.energy_range_map[hv]
                if emin < c.integral < emax:
                    etrue = etrue0 #* np.random.normal(1,self.sigma0)
                    if 5<z<48:
                        ztrue = z #np.random.uniform(z-2,z+2)  #z+np.random.normal(0,self.sigmaz0)
                    elif z<=5:
                        ztrue = np.random.uniform(z-5,z+0.5)
                    elif z>=48:
                        ztrue = np.random.uniform(z-0.5,z+5)
                    else:
                        print ("Waring : z = %f not foreseen in measurement map ",z)
                ret["sc_trueint"].append(etrue)
                ret["sc_truez"].append(ztrue)
                ret["sc_hv"].append(hv)
        for V in self.vars:
            self.out.fillBranch(V,ret[V])

        return True

trueEnergy = lambda : RegressionTrainingVarsLime()

