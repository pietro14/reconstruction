# USAGE: python gbr_trainer.py ~/Work/data/cygnus/RECO/lime2021/v2/fe55/reco_run04455_3D.root params_gbrtrain.txt

import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
import uproot
import pandas as pd

firstEvent = 0
lastEvent = 100

rfile = "~/cernbox/CYGNUS/reco/lime/merged_winter22/reco_runs5X.root"
events = uproot.open("~/cernbox/CYGNUS/reco/lime/merged_winter22/reco_runs5X.root")
variables = ["sc_xmean","sc_ymean","sc_xmean/sc_ymean"]
data_main = events["Events"].arrays(variables,library="pd",entry_start=firstEvent,entry_stop=lastEvent)

print (data_main)

#data_sel = data[(data['sc_xmean']>1000)]

frrfile = rfile.replace("/reco_","/friends/reco_").replace(".root","_Friend.root")
events_friends = uproot.open(frrfile)
variables_friend = ["sc_trueint"]
data_friend = events_friends["Friends"].arrays(variables_friend,library="pd",entry_start=firstEvent,entry_stop=lastEvent)

print (data_friend)

data = pd.concat([data_main,data_friend],axis=1)

print ("===data before")
print(data)
print (" ===== ")

#print ("===data after")
#data = data[(data['sc_xmean']>1)]
#print(data)
#print (" ===== ")

print(data.iloc[[0]])

print (" AAAA ===== ")
data_arr = data.to_numpy()
Nevents = len(data_arr)
conc = np.stack(data.to_numpy()[0],axis=-1)
print ("conc 0")
print (conc)
for i in range(1,Nevents):
    event = np.stack(data.to_numpy()[i],axis=-1)
    conc = np.vstack([conc,event])
print ("conc alla fine")
print (conc)
print ("--->")

data_pd = pd.DataFrame(conc, columns=variables+variables_friend)

print (" DDDDDDD ===== ")
print (data_pd)
data_sel = data_pd[(data_pd['sc_xmean/sc_ymean']>1)]
print(data_sel)

#print ( np.concatenate(np.concatenate(data.values)).reshape(N,-1).T )

#conc = np.concatenate(np.concatenate(data.values)).reshape(N,-1).T
#conc = np.concatenate(np.concatenate(data.values)).reshape(N,-1).T
#print (conc)

#conc = np.stack(data.iloc[[0]],axis=0)
#print (conc)


#print (data_sel)
#conc = pd.concat([data,df], ignore_index=True, axis=0, sort=False)

        
# tfile = ROOT.TFile.Open(rfile)
# tree = tfile.Get(self.tree_name)
# if friendrfile:
#     tree.AddFriend("Friends",friendrfile)

# # so target is always the first variable
# variables = [self.target] + self.var.split("|")
# print("List of variables = ",variables)
# dataset = tree2array(tree,variables,object_selection={self.cuts_base : variables})
# tfile.Close()

# conc = np.stack(dataset[0],axis=-1)
# for i in range(1,len(dataset)):
#     conc = np.append(conc,np.stack(dataset[i],axis=-1),axis=0)
# X = conc[:,1:]
# y = conc[:,0]
# return X,y
    
