import numpy as np

def getContours(xbox,ybox):
    
    ymi = []
    yma = []
    xmi = []
    xma = []
    xu = np.unique(xbox)
    for i in xu:
        if ybox[xbox == i].size > 0:
            mi = np.min(ybox[xbox == i])
            ymi.append(mi)
            xmi.append(i)
            ma = np.max(ybox[xbox == i])
            yma.append(ma)
            xma.append(i)
            
    xri = np.concatenate((xmi, xma[::-1],[xmi[0]]))
    yri = np.concatenate((ymi, yma[::-1], [ymi[0]]))
    return xri,yri

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def inputFile(numrun,filedir,formattype):

    ## Loading image for analysis
    # otherwise use 'mid'
    if formattype == 'h5':
        imagename     = 'histogram_Run'
    else:
        imagename     = 'histograms_Run'

    filename = '%s%s%s.root' % (filedir,imagename,numrun)
    return filename