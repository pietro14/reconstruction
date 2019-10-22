import numpy as np
import glob, os
import re
import sys


def cluInfo(clusters,points,Ci,image,m_image,scale):
    
    labels = clusters.labels_
    pos    = points[labels == Ci]
    tag    = clusters.tag_[clusters.labels_== Ci][0]
    Xi     = list(pos[:,0].astype(int))
    Yi     = list(pos[:,1].astype(int))
    
    Is = []
    Ib = []
    ## to get X and Y in full dimension
    Xn = np.zeros(len(Xi)*scale*scale,dtype=int)
    Yn = np.zeros(len(Yi)*scale*scale,dtype=int)
    for i in range(0, len(Xi)):    # loop on single cluster value
        y0      = int(Xi[i]*scale)
        x0      = int(Yi[i]*scale)
        factor = scale*scale
        
        Yn[factor*i:factor*(i+1)] = np.reshape(np.reshape(np.arange((y0),(y0+scale)).repeat(scale,axis = 0),[scale,scale]).T,factor)
        Xn[factor*i:factor*(i+1)] = np.arange((x0),(x0+scale)).repeat(scale,axis = 0)
        
        Is.extend(list(image[(y0):(y0+scale),(x0):(x0+scale)].reshape(scale*scale).astype(int)))
        Ib.extend(list(m_image[(y0):(y0+scale),(x0):(x0+scale)].reshape(scale*scale).astype(int)))
    
    Xf = list(Xn.astype(int))
    Yf = list(Yn.astype(int))
    del Xi,Yi,Xn,Yn
      
    return Xf, Yf, Is, Ib, tag

        
def openTable(file):
    with open(file) as f:
        datain = eval(f.read())
    return datain

def saveTable(file,data):
    with open(file, 'w') as f:
        f.write(repr(data))
        

def colorbar(mappable):
    # plot colorbars        
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

def plot2hist(vari, bins = [20,40,100], liml = 0,limr = 50, label = '', scale = '', unity = '', density = True, logx = False, logy = False):
    ## Function to show in the same histogram the three categories
    
    # vari     = is a 3xN List which row has the information of one category       - List
    # nins     = is the number of bins to construct the histogram                  - int
    # liml     = Xlim min                                                          - float
    # limr     = Xlim max                                                          - float
    # label    = is the Xlabel to show in the plot                                 - string
    # scale    = is 'm' - mili, 'u' - micro, 'n' - nano,..                         - char
    # unity    = is the metric unit                                                - char
    # density  = is the flag to set the histogram to show the density or not       - Boolean
    # logx     = is the flag to set the X axis to log on the histogram or not      - Boolean
    # logy     = is the flag to set the Y axis to log on the histogram or not      - Boolean
    
    import matplotlib.pyplot as plt
    v = vari.copy()
    
    if scale == 'k':
        for i in range(0,3):
            v[i] = v[i]/(10**3)
    elif scale == 'M':
        for i in range(0,3):
            v[i] = v[i]/(10**6)
    elif scale == 'G':
        for i in range(0,3):
            v[i] = v[i]/(10**9)
    elif scale == 'T':
        for i in range(0,3):
            v[i] = v[i]/(10**12)
    elif scale == 'm':
        for i in range(0,3):
            v[i] = v[i]*(10**3)
    elif scale == 'u':
        for i in range(0,3):
            v[i] = v[i]*(10**6)
    elif scale == 'n':
        for i in range(0,3):
            v[i] = v[i]*(10**9)
    elif scale == 'p':
        for i in range(0,3):
            v[i] = v[i]*(10**12)
    else:
        for i in range(0,3):
            v[i] = v[i]
    
    
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111)
    
    plt.hist(v[0], bins=bins[0], fc='r', alpha = 0.7, density=density)
    plt.hist(v[1], bins=bins[1], fc='b', alpha = 0.7, density=density)
    plt.hist(v[2], bins=bins[2], fc='darkorange', alpha = 0.7, density=density)
    plt.xlim([liml, limr])
    
    if logx:
        plt.xscale("log")
    if logy:
        plt.yscale("log")
    if density:
        plt.ylabel('Probability')
    else:
        plt.ylabel('Counts')
    plt.xlabel(label + '(' + scale + unity + ')',fontsize=18)
    plt.legend(['Recoils', 'Soft Electrons', 'MeV Electrons'],prop={'size': 18})
    plt.show()
    plt.close
    
def getTaggedVariable(vari,col):
    ## Function get the three categories on the same variable
    #    INPUT    
    # vari = is the DataFrame Pandas with all the information               - DataFrame
    # col  = is the name of the coloumn wanted                              - String
    #    OUTPUT
    # v = is a 3xN List which row has the information of one category       - List
    
    v = [np.array(vari[col][vari.Tag == 'l']), np.array(vari[col][vari.Tag == 'm']), np.array(vari[col][vari.Tag == 's'])]
    return v

def plotMesh(X,Y,Z,el,graus):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.view_init(elev=el, azim=graus)
    plt.show()
    return ax

#### Tools to rotate the cluster

def rotate(oX, oY, pX, pY, angle):
    from math import sin
    from math import cos
    from numpy import array
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox = oX
    oy = oY
    px = array(pX)
    py = array(pY)

    qx = ox + cos(angle) * (px - ox) - sin(angle) * (py - oy)
    qy = oy + sin(angle) * (px - ox) + cos(angle) * (py - oy)
    
    return qx.tolist(), qy.tolist()

def getAngle(X,Y):
    from numpy import polyfit
    from numpy import poly1d
    from numpy import arctan
    
    # - - - - Reta 0,0
    xo = [-2048,2048]
    yo = [0.001,0.001]
    
    zo = polyfit(yo,xo, 1)
    fo = poly1d(zo)    
    m1 = fo.c[0] 
    
    z = polyfit(X,Y, 1)
    func = poly1d(z) 
    m2 = func.c[0]
    
    
    
    angle = arctan(m1-m2/(1-m1*m2))
    
    return angle

def getC(X,Y):
    from numpy import polyfit
    from numpy import poly1d
    from numpy import arctan
    
    z = polyfit(X,Y, 1)
    func = poly1d(z)    
    
    return func.c[0]

def get_sliceleng(X,Y,pieces):
    # Function to get the mean length of the cluster
    # in X or Y direction.
    pieces = pieces
    
    newX = np.array(X) # Direction of the slices
    newY = np.array(Y) # Direction of the Mean Length
    
    slices = np.linspace(np.min(newX),np.max(newX),pieces)
    meanLY = np.zeros([(pieces-1),],dtype=float)

    for i in range(0,(pieces-1)):
    
        y = newY[(newX > slices[i]) & (newX < slices[i+1])]
        meanLY[i] = np.max(y) - np.min(y)
    return meanLY

def plot1hist(variable, bins, liml = 0, limr = 50, label = '', density = True, logx = False, logy = False):
    ## Function to show one variable on a histogram
    
    # variable = is a 1xN List with the variable information                       - List
    # nins     = is the number of bins to construct the histogram                  - int
    # nsd      = is the multiplication factor of the sigma to define the Xlim min  - float
    # nse      = is the multiplication factor of the sigma to define the Xlim max  - float
    # label    = is the Xlabel to show in the plot                                 - string
    # density  = is the flag to set the histogram to show the density or not       - Boolean
    # logx     = is the flag to set the X axis to log on the histogram or not      - Boolean
    # logy     = is the flag to set the Y axis to log on the histogram or not      - Boolean
    
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10,7))
    ax  = fig.add_subplot(111)
    
    v2  = variable
    e   = np.size(v2)
    m   = np.mean(v2[(v2 != 0) & (np.isnan(v2) == False)])
    s   = np.std(v2[(v2 != 0) & (np.isnan(v2) == False)])
    
    plt.hist(variable, bins=bins, fc='r', alpha = 0.7, density=density)
    plt.xlim([liml, limr])
    
    if logx:
        plt.xscale("log")
    if logy:
        plt.yscale("log")
    if density:
        plt.ylabel('Probability')
    else:
        plt.ylabel('Counts')
    plt.xlabel(label,fontsize=18)
    #plt.legend(['Recoils', 'Soft Electrons', 'MeV Electrons'],prop={'size': 18})

    textstr = '\n'.join((
        r'Entries $=%d$' % (e, ),
        r'Mean $=%.2f$' % (m, ),
        r'Std Dev$ =%.2f$' % (s, )))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)

    # place a text box in upper left in axes coords
    plt.text(0.7, 0.9, textstr, fontsize=14,
            verticalalignment='top',transform=ax.transAxes, bbox=props)
    
    plt.show()

def pl3d(X,Y,Z,azim=0, bottom = 80):

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt



    fig = plt.figure(figsize = (10,7))
    ax  = fig.add_subplot(111, projection='3d')

    ax.scatter(X, Y, Z, c = 'r', marker = 'o')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_zlim(bottom = bottom)
    ax.view_init(elev = 0., azim = azim)

    plt.show()
    
    
def plotLineCluster(ax,sx,ex,sxy,exy):
    import matplotlib.pyplot as plt
    
    #fig = plt.figure(figsize=(7,7))
    ax.plot([sx, ex], [sxy, exy], 'b-')
    ax.plot(ex, exy, 'rs')
    ax.plot(sx, sxy, 'k<')
    
    circle1=plt.Circle((895,920),900,color='r',fill=False)
    plt.gcf().gca().add_artist(circle1)
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_xlim([-50,2078])
    ax.set_ylim([-50,2078])
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.axis('square')
    
    ax.legend(['Line','Start Point', 'End Point'])
    
def getSignalCircle(featuresL, zX = 895, zY = 920, r = 900, mr = 20, gr = 10, flag = False):
    
    ind   = np.array(featuresL.index[featuresL.Image < 210])
    ex    = np.array(featuresL.EndX[featuresL.Image < 210])
    exy   = np.array(featuresL.EndXy[featuresL.Image < 210])
    
    teste = np.sqrt((ex-zX)**2 + (exy-zY)**2)
    indn  = np.where((teste > (r-mr)) & (teste < (r+gr)))
    if flag:
        return ind[indn],indn[0]
    else:
        return ind[indn]

def getSignalEndInCircle(featuresL, zX = 895, zY = 920, r = 900, gr = 10):
    
    ind   = np.array(featuresL.index[featuresL.Image < 210])
    ex    = np.array(featuresL.StartX[featuresL.Image < 210])
    exy   = np.array(featuresL.StartXy[featuresL.Image < 210])
    
    teste = np.sqrt((ex-zX)**2 + (exy-zY)**2)
    indin  = np.where(teste < (r-gr))
    indout  = np.where(teste >= (r-gr))
    
    return ind[indin],ind[indout]

def plot_shapeprofile(X,Y,L,P = 0, px = 10, p = 60, debug = False, bp = False):
    from scipy import stats
    from scipy.optimize import curve_fit
    from astropy.modeling.models import Gaussian1D
    import matplotlib.pyplot as plt
    from toolslib import get_shapeprofile

    
    L    = np.array(L) - np.array(P)
    newX = np.array(X) # Direction of the slices
    newY = np.array(Y) # Direction of the Mean Length

    pieces = np.int(np.round((np.max(newX)-np.min(newX))/px))
    slices = np.linspace(np.min(newX),np.max(newX),pieces)
    piecesY = np.int(np.round(np.max(newY))-np.round(np.min(newY)))
    
    xm = np.zeros([pieces-1,],dtype=float)
    zm = np.zeros([pieces-1,],dtype=float)
    zs = np.zeros([pieces-1,],dtype=float)
    zM = np.zeros([pieces-1,],dtype=float)
    
    matrixY = np.zeros([pieces-1,piecesY+1],dtype=float)
    ct = np.zeros([piecesY+1],dtype=float)
    
    if debug:
        fig, ax = plt.subplots(1,2,figsize=(18, 6))
    
    for i in range(0,(pieces-1)):

        y = newY[(newX > slices[i]) & (newX < slices[i+1])]
        #x = newX[(newX > slices[i]) & (newX < slices[i+1])]
        z = L[(newX > slices[i]) & (newX < slices[i+1])]
        
        xm[i] = (slices[i] + slices[i+1])/2
        
        if debug:
            iy, uy = get_shapeprofile(axi = y,light = z,fig = fig, ax = ax, debug = debug)
        else:
            iy, uy = get_shapeprofile(axi = y,light = z,fig = None, ax = None, debug = debug)
        
        ii = uy-np.round(np.min(newY))
        
        matrixY[i,ii.astype(int)] = iy
        
        ###
        perc = (1-(p/100))/2
        mm = uy[-1]-uy[0]
        idp = np.where((uy >= uy[0]+(mm*(perc))) & (uy <= uy[-1]-(mm*(perc))))
        ###
        
        zm[i] = np.mean(iy)
        zM[i] = np.max(iy)
        zs[i] = np.sum(iy[idp])
    
    if debug:
        fig1, ax1 = plt.subplots(2,2,figsize = (18, 18))
        ax1[1,1].plot(xm,zs,'x:')
        ax1[1,1].set_title("Sum of X Profile")
        ax1[0,1].plot(xm,zm,'o:')
        ax1[0,1].set_title("Mean of X Profile")
        ax1[1,0].plot(xm,zM,'s:')
        ax1[1,0].set_title("Max of X Profile")
    
    
    for jj in range(0,piecesY+1):
        aux = matrixY[:,jj]
        ct[jj] = np.size(aux[aux != 0])
    
    # - - - - - - - - - - - - - - - - - - - - -
    Yproj = np.mean(matrixY,axis = 0)
    errY = np.std(matrixY,axis = 0)/np.sqrt(ct)
    xy    = np.arange(0,piecesY+1)+np.min(newY)

    # weighted arithmetic mean (corrected - check the section below)
    mean  = sum(xy * Yproj) / sum(Yproj)
    sigma = np.sqrt(sum(Yproj * (xy - mean)**2) / sum(Yproj))
    sig   = np.std(Yproj)    
    
    def Gauss(x, a, x0, sigma):
        return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
    
    if debug:
        ax1[0,0].errorbar(xy, Yproj, yerr = errY, fmt = 'b+:', label='data')
    try:
        popt,pcov = curve_fit(Gauss, xy, Yproj, p0 = [max(Yproj), mean, sigma])
        if debug:
            ax1[0,0].plot(xy, Gauss(xy, *popt), 'r-',
               label='Gauss fit   \nAmpl     = %.1f\nMean    = %.1f\nSigma   = %.1f' %
               (popt[0], popt[1], popt[2]))
        g1     = Gaussian1D(np.max(Gauss(xy, *popt)), mean = popt[1], stddev = popt[2])
        #width  = popt[2]*sg
    except:
        print ("fit error")
    if debug:
        ax1[0,0].plot(np.linspace(popt[1]-(g1.fwhm/2),popt[1]+(g1.fwhm/2),10),
                      np.ones(10,dtype='float')*(g1.amplitude/2),
                      '-k',label = 'FWHM    = %.1f' % (g1.fwhm))
        ax1[0,0].set_xlabel('Y [pixel]')
        ax1[0,0].set_ylabel('average light profile [ph]')
        ax1[0,0].minorticks_on()
        ax1[0,0].legend()
        plt.show()
        fig1.hold
    
        
    ###### X information ############
    widthY = popt[2]
    widthX = np.max(xm)-np.min(xm)
    peakX  = np.max(zs)
    meanX  = np.mean(zs)
    if bp:
        return xm,zs
    else:
        return widthY, widthX, peakX, meanX


def get_shapeprofile(axi, light, fig, ax, debug = False):
    from scipy import stats
    import matplotlib.pyplot as plt
    
    uy = np.unique(np.round(axi))
    iy = np.zeros([np.size(uy),], dtype = float)
    my = np.zeros([np.size(uy),], dtype = float)

    for jj in range(0,np.size(uy)):
        ind = np.where(np.round(axi) == uy[jj])
        iy[jj]  = np.sum(light[ind])
        my[jj]  = stats.mode(light[ind])[0]

    if debug:
        #fig, ax = plt.subplots(1,2,figsize=(18, 6))
        ax[0].plot(uy,iy,'-x')
        ax[0].set_title("Sum of Y Profile")
        ax[1].plot(uy,my,'o')
        ax[1].set_title("Mode of Y Profile")
        fig.hold
    return iy,uy.astype(int)

def plottingCluster(df,colhead,cluN,x_resolution,y_resolution):
    import matplotlib.pyplot as plt
    from toolslib import colorbar
    
    index = df.index[cluN]
    Run   = df[colhead[0]][cluN]
    Nim   = df[colhead[1]][cluN]
    Xi    = df[colhead[3]][cluN]
    Yi    = df[colhead[4]][cluN]
    Lp    = df[colhead[5]][cluN]
    Lb    = df[colhead[6]][cluN]

    matrix = np.zeros([y_resolution,x_resolution],dtype=int)
    matrixb = np.zeros([y_resolution,x_resolution],dtype=int)
    
    matrix[Yi,Xi]=Lp
    matrixb[Yi,Xi]=Lb
    
    fig = plt.figure(figsize=(15,15))
    ax  = plt.gca()
    
    iax = ax.imshow(matrix,cmap="viridis", vmin=85,vmax=130)
    ax.set_ylim(np.max(Yi),np.min(Yi))
    ax.set_xlim(np.min(Xi),np.max(Xi))
    ax.set_title('%d - Run %d - # of Image %d' % (index, Run, Nim))
    colorbar(iax)
    plt.show(block=False)
    plt.close
    
def xstart(x):
    yo = np.concatenate([np.linspace(200,300,80), np.linspace(450,550,30), np.linspace(600,800,25), np.linspace(800,1100,15)])
    xo = np.linspace(330,0,np.size(yo))
    
    zo = np.polyfit(yo,xo, 6)
    fo = np.poly1d(zo)
    return fo(x)