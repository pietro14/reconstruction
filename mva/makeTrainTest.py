from __future__ import print_function

import ROOT as r 
import numpy as np
import pickle
import math
from keras.utils import np_utils
import tqdm

from multiprocessing import Pool


testCut   = lambda ev: ev.event%5==0
trainCut  = lambda ev: ev.event%5!=0


featureList = [
    "density"          ,
    "length"           ,
    "slimness"         ,
    "lgsigma"      ,
    "tgsigma"      ,
    "latrms"           ,
    "longrms"          ,
    "size"             ,
    "nhits"            ,
]


features = {

    "density"          : lambda ev,isc : ev.sc_integral[isc]/ev.sc_nhits[isc] if  ev.sc_integral[isc]>0 else -9, # 0
    "length"           : lambda ev,isc : ev.sc_length[isc],                                                      # 1
    "slimness"         : lambda ev,isc : ev.sc_width[isc]/ev.sc_length[isc],                                     # 2
    "lgsigma"          : lambda ev,isc : ev.sc_lgausssigma[isc],                                                 # 3
    "tgsigma"          : lambda ev,isc : ev.sc_tgausssigma[isc],                                                 # 4
    "latrms"           : lambda ev,isc : ev.sc_latrms[isc],                                                      # 5
    "longrms"          : lambda ev,isc : ev.sc_longrms[isc],                                                     # 6
    "size"             : lambda ev,isc : ev.sc_size[isc],                                                        # 7 
    "nhits"            : lambda ev,isc : ev.sc_nhits[isc],                                                       # 8

    }

def preselection(ev,isc):
    pixw = 0.152 # pixel width in mm
    NX = 2304
    framesize = 500
    ## first select a strict non-noisy region
    good = ev.sc_xmean[isc]>framesize and ev.sc_xmean[isc]<NX-framesize and ev.sc_ymean[isc]>framesize and ev.sc_ymean[isc]<NX-framesize
    good = good and ev.sc_lgausssigma[isc]<10
    if ev.sc_length[isc]*pixw > 50. or ev.sc_width[isc]/ev.sc_length[isc]<0.4 or ev.sc_width[isc]*pixw>6.54 or pixw*ev.sc_tgausssigma[isc]>0.3:
        good = False
    return good
    
classes = {
    ## NR add a preselection to cut all the very low energy stuff that I don't understand. The rest looks like data
    'nr'      : { 'cut': lambda ev,isc: preselection(ev,isc) and ev.sc_integral[isc]/ev.sc_nhits[isc]>14, 'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
    'er'      : { 'cut': lambda ev,isc: preselection(ev,isc) and ev.sc_integral[isc]/ev.sc_nhits[isc]>12, 'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
    'other'   : { 'cut': lambda ev,isc: preselection(ev,isc), 'lst_train' : [], 'lst_test' : [] , 'lst_y_train' : [], 'lst_y_test' : [] },
}

sampleDir='/Users/emanuele/Work/data/cygnus/RECO/lime2020/v6/'

nrStuff = []
erStuff = []
othStuff = []

# this is data, but it contains in both cases other bkg
nrStuff.extend( ['ambe_lime_v6_runs3737_3791.root'] )
erStuff.extend( ['fe_lime_v6_runs3686_3691.root'] )

# this is simulation with pure signal
#nrStuff.extend( ['CYGNO_60_40_He_NR_10_keV_30cm_SRIM_IDAO_lime2021.root'] )
#erStuff.extend( ['CYGNO_60_40_ER_6_keV_30cm_IDAO_lime2021.root'] )
othStuff.extend( ['cosmics_afterambe_lime_v6_runs3792_3794.root'] )

def toNumpy(maxEntries,task):
    print('starting', task)
    fil, typs = task
    print('List of features for', featureList + eval('featureList'))
    tfile = r.TFile(fil); ttree = tfile.Events
    results = {}
    for ty in typs: 
        results[ty + '_test']  = []
        results[ty + '_train'] = []

    print("start looping on file ",fil)
    for iev,ev in enumerate(ttree):
        if iev%1000==0: print('Processing event ',iev,'...')
        if iev>maxEntries: break
        tstr = 'test' if testCut(ev) else 'train'
        for isc in range(ev.nSc):
           for ty in typs:
               if classes[ty]['cut'](ev,isc):
                   results[ty+'_'+tstr].append([ features[s](ev,isc) for s in (featureList) ])
    tfile.Close()
    print('finishing', task)
    return results




if __name__ == "__main__":

    from functools import partial

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    # common options, independent of the flavour chosen
    parser.add_option("--max-entries",     dest="maxEntries", default=1000000000, type="int", help="Max entries to process in each tree")
    (options, args) = parser.parse_args()

    tasks = []
    print('Setting up the tasks')
    for samp in nrStuff:
        tasks.append( (sampleDir+'/'+samp, ['nr']) )

    for samp in erStuff:
        tasks.append( (sampleDir+'/'+samp, ['er']) )
    
    for samp in othStuff:
        tasks.append( (sampleDir+'/'+samp, ['other']) )
    
    print('Going to run the big thing')
    print("max entries = ",options.maxEntries)
    p =  Pool(4)
    func = partial(toNumpy,options.maxEntries)
    results = list(tqdm.tqdm(p.imap(func, tasks), total=len(tasks)))

    print('Now putting everything together')

    types = ['nr', 'er', 'other']
    for result in results: 
        for ty in types:
            if ty+'_train' in result:
                classes[ty]['lst_train'].extend( result[ty+'_train'])
                classes[ty]['lst_test' ].extend( result[ty+'_test'])

            
    print('Setting the indices')
    toDump = {} 
    for i, ty in enumerate(types):
        classes[ty]['lst_train'  ] = np.asarray(classes[ty]['lst_train'])
        classes[ty]['lst_y_train'] = i*np.ones((classes[ty]['lst_train'].shape[0],1))
        classes[ty]['lst_test'   ] = np.asarray(classes[ty]['lst_test'])
        classes[ty]['lst_y_test' ] = i*np.ones((classes[ty]['lst_test'].shape[0],1))

    train_x = np.concatenate( tuple( [classes[ty]['lst_train'] for ty in types] ), axis=0)
    train_y = np_utils.to_categorical( np.concatenate( tuple( [classes[ty]['lst_y_train'] for ty in types] ), axis=0), len(classes))
    test_x = np.concatenate( tuple( [classes[ty]['lst_test'] for ty in types] ), axis=0)
    test_y = np_utils.to_categorical( np.concatenate( tuple( [classes[ty]['lst_y_test'] for ty in types] ), axis=0), len(classes))
    toDump['train_x'] = train_x
    toDump['train_y'] = train_y
    toDump['test_x' ] = test_x
    toDump['test_y' ] = test_y

    ### dump to file
    print ('dump to vars.pkl now...')
    pickle_out = open('vars.pkl','wb')
    pickle.dump( toDump, pickle_out)
    pickle_out.close()
