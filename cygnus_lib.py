#
# Image analisys CYGNUS-RD Python Library
# G. Mazzitelli 2017 
# rev oct 2018 - swift direct access 
#

import numpy as np
import glob, os
import re
import sys


# ############################################################################################################################ #
# DATA Archive
# ############################################################################################################################ #

def logbookInfo(file, run):
    logbook = np.loadtxt(file, dtype=str, delimiter='\t')
    for i in range(1, len(logbook[:,0])):
        if logbook[:,0][i] == str(run):
            data = logbook[i]
    return logbook[0], data

def strigHMS2time(str):
    import datetime
    pt=datetime.datetime.strptime(str,'%H:%M:%S')
    sec=pt.second+pt.minute*60+pt.hour*3600
    return sec

def ScopeHeader(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    f.close()
    size_sco      = int(str.split(lines[1])[3])-1
    start_sec_sco = strigHMS2time(str.split(lines[3])[2])
    TA_start_sco  = lines[5]
    TA_end_sco    = lines[size_sco+5]
    T_start_sco   = float(str.split(TA_start_sco)[0])
    T_end_sco     = float(str.split(TA_end_sco)[0])

    return size_sco, start_sec_sco, T_start_sco, T_end_sco

def ScopeHeaderPath(path, ch):
    sch = 'C'+str(ch)+'wave'
    if path.find('BTF_2017-2') or dataSelection == 'LABOct2017':
        sch = 'C'+str(ch)+'Run'
    file_sco = path + sch + '%05d' % (0) +'.txt'
    size_sco, start_sec_sco, T_start_sco, T_end_sco = ScopeHeader(file_sco)

    return size_sco, start_sec_sco, T_start_sco, T_end_sco

def ReadScope(path, ch):
    file_in_dir=os.listdir(path)
    sch = 'C'+str(ch)+'wave'
    if path.find('BTF_2017-2')  or dataSelection == 'LABOct2017':
        sch = 'C'+str(ch)+'Run'
    if sch in str(file_in_dir):
        nfile = np.size(filter(lambda x: sch in x, file_in_dir))
        file_sco = path + sch + '%05d' % (0) +'.txt'
        size_sco, start_sec_sco, T_start_sco, T_end_sco = ScopeHeader(file_sco)
        t = np.empty((nfile, size_sco), dtype=np.double)
        a = np.empty((nfile, size_sco), dtype=np.double)
        for i in range(0, nfile):
            file_sco = path + sch + '%05d' % (i) +'.txt'
            t[i,], a[i,] = np.loadtxt(file_sco, delimiter=' ', skiprows=6, usecols=(0, 1), 
                                      unpack=True, dtype='double')
    return t, a

def ReadScopeTrace(path):

#    t, a = np.loadtxt(path, delimiter=' ', skiprows=6, usecols=(0, 1), 
#                       unpack=True, dtype='double')
    t, a = np.genfromtxt(path, delimiter=' ', skip_header=7, unpack=(0, 1))
    return t, a

def file2PathCygnus(dataSelection, fileNumber, ftype):
    #  return ftyp file path RUN, H5, TS, LOG, SCO in dataSelection 
    itype = ['RUN', 'H5', 'TS', 'LOG', 'SCO', 'TMP']
    base_path = 'Data/'+dataSelection+'/'
    RUN_path = 'Run%03d/' % np.int(fileNumber)
    H5_path   = base_path + 'Data_Camera/H5/'
    TS_path   = base_path + 'Data_Camera/TS/'
    LOG_path  = base_path + 'LOG/'
    SCO_path  = base_path + 'Data_Scope/'
    TMP_path  = base_path + 'TMP/'
    path = [RUN_path, H5_path, TS_path, LOG_path, SCO_path, TMP_path]
    return path[itype.index(ftype)]

def file2FullPathCygnus(dataSelection, fileNumber, ftype):
    itype = ['RUN', 'H5', 'TS', 'LOG', 'SCO', 'TMP']
    est = ['.HIS', '/', '', '.LOG', '/', '']
    return file2PathCygnus(dataSelection, fileNumber, ftype)+'Run%03d' % np.int(fileNumber) + est[itype.index(ftype)]

def imageFile2FullPathCygnus(dataSelection, fileNumber, traccia):
    return file2FullPathCygnus(dataSelection, fileNumber, 'H5')+'run%03d' % np.int(fileNumber)+('-%04d.h5' % traccia)

def scopeFile2FullPathCygnus(dataSelection, fileNumber, traccia, ch):
    sch = 'C'+str(ch)+'wave'
    if dataSelection == 'BTF_2017-2' or dataSelection == 'LABOct2017':
        sch = 'C'+str(ch)+'Run'
    return file2FullPathCygnus(dataSelection, fileNumber, 'SCO') + sch + '%05d' % (traccia) +'.txt'


# ############################################################################################################################ #
# swift storage read and write
# ############################################################################################################################ #

def swift_auth():
    # https://docs.openstack.org/python-swiftclient/latest/service-api.html
    # https://docs.openstack.org/python-swiftclient/

    import swiftclient
    from keystoneauth1 import session
    from keystoneauth1.identity import v3

    
    OS_REGION_NAME='lnf'
    OS_USER_DOMAIN_NAME='default'
    OS_PROJECT_DOMAIN_NAME='default'
    OS_PROJECT_NAME='cygnus-default'
    OS_IDENTITY_API_VERSION='3'
    OS_PASSWORD='bayRonPeOcan9Quiufvecfevesyailb7'
    OS_AUTH_TYPE='password'
    OS_AUTH_STRATEGY='keystone'
    OS_AUTH_URL='https://keystone.cloud.infn.it:5000/v3/'
    OS_USERNAME='cygnus'
    OS_STORAGE_URL='https://swift.cloud.infn.it:8080/v1/AUTH_1e60fe39fba04701aa5ffc0b97871ed8'


    _auth = v3.Password(
        user_domain_name    = OS_USER_DOMAIN_NAME,
        project_domain_name = OS_PROJECT_DOMAIN_NAME,
        project_name        = OS_PROJECT_NAME,
        username            = OS_USERNAME,
        password            = OS_PASSWORD,
        auth_url            = OS_AUTH_URL
    )
    _os_options={
        'region_name' : OS_REGION_NAME, 
        'object_storage_url': OS_STORAGE_URL
    }
    # Create session
    keystone_session = session.Session(auth = _auth)

    # Create swiftclient Connection
    swift = swiftclient.Connection(session      = keystone_session, 
                                    auth_version = OS_IDENTITY_API_VERSION,
                                    os_options   = _os_options
                                    )
    return swift

def swift_read_image_h5(file):
    # https://www.getdatajoy.com/learn/Read_and_Write_HDF5_from_Python#Reading_Data_from_HDF5_Files
    import numpy as np
    import h5py
    import os
    swift = swift_auth()
    obj_tuple = swift.get_object("Cygnus", file)
    tmpname = "./tmp." + str(os.getpid()) + ".h5"
    with open(tmpname, 'wb') as my_tmp:
            my_tmp.write(obj_tuple[1])
    image = read_image_h5(tmpname)
    try:
        os.remove(tmpname)
    except OSError:
        pass
    return image

def swift_listdir(dirname):
    swift = swift_auth()
    fileindir=[]
    for data in swift.get_container("Cygnus", full_listing=True)[1]:
        if dirname in str(data):
            fileindir.append(data['name'])
    return fileindir


# ############################################################################################################################ #
# TOOLS
# ############################################################################################################################ #

def read_image_h5(file):
    # https://www.getdatajoy.com/learn/Read_and_Write_HDF5_from_Python#Reading_Data_from_HDF5_Files
    import numpy as np
    import h5py
    with h5py.File(file,'r') as hf:
        data = hf.get('Image')
        np_data = np.array(data)
    return np_data

def write_image_h5(file, m1):
    import numpy as np
    import h5py
 
    with h5py.File(file, 'w') as hf:
        hf.create_dataset('Image', data=m1)
    return

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def smooth(y, box_pts):
    import numpy as np
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def OverTh2Array(ArrayIn, Th):
    # return x,y,status array find in ArrayIn when Threshold is passed
    OverTh    = False
    ThArr = []
    
    for i in range (0, len(ArrayIn)):
        if ArrayIn[i]>=Th and not OverTh:
            OverTh     = True
            ThArr.append([i, ArrayIn[i], OverTh])
        if ArrayIn[i]<Th and OverTh:
            OverTh     = False
            ThArr.append([i, ArrayIn[i], OverTh])
    return ThArr

def UnderTh2Array(ArrayIn, Th):
    UnderTh    = False
    ThArr = []

    for i in range (0, len(ArrayIn)):
        if ArrayIn[i]<=Th and not UnderTh:
            UnderTh     = True
            ThArr.append([i, ArrayIn[i], UnderTh])
        if ArrayIn[i]>Th and UnderTh:
            UnderTh     = False
            ThArr.append([i, ArrayIn[i], UnderTh])
    return  ThArr

# ############################################################################################################################ #
# Clustering
# ############################################################################################################################ #

def NNClustering(points, thC):
# poi: vettore dei punti X,Y nel piano
    
    import numpy as np
    C = np.zeros((len(points), 4))
    for i in range(0,  len(points)):
        C[i]=[i, 0, points[i,1], points[i,0]]
    
    NofC  = 0
    NeC   = 0
    nloop = 0
    while True:
        i = 1
        nordered = 0
        while (i < len(C)):
            j = 0
            while (j < i):
                sBreak = False
                if abs(C[j,2]-C[i,2])<thC and abs(C[j,3]-C[i,3])<thC:   # close point i to j
                    NofCj = C[j,0]
                    NeCj  = (C[C[:,0]==NofCj][:,1]).max()
                    NofCi = C[i,0]
                    NeCi  = (C[C[:,0]==NofCi][:,1]).max()
                    
                    if NofCi != NofCj:
                        if NeCi == 0:
                            C[i,0]= NofCj
                            C[i,1]= NeCj+1
                        else:
                            if NofCi>NofCj:
                                Ci = np.where(C[:,0]==NofCi)
                                for iCi in range(0, len(Ci)):
                                    C[Ci[iCi], 0] = NofCj
                                    C[Ci[iCi], 1] = NeCj+iCi+1
                            else: 
                                Cj = np.where(C[:,0]==NofCj)[0]
                                for iCj in range(0, len(Cj)):
                                    C[Cj[iCj], 0] = NofCi
                                    C[Cj[iCj], 1] = NeCi+iCj+1
                        sBreak = True
                        nordered += 1
                        break
                j += 1
            i +=1
        if nordered == 0:
            break
        nloop += 1

    sorted_C = np.array(sorted(C, key=lambda x:x[0]))
    return sorted_C

# HDBSCAN clustering. Uses HDBSCAN to generate clusters and returns an output matrix formatted as in NNClustering, 
# i.e. each row follows the format [cluster label, intracluster label, yval, xval]. The matrix has dimensions 
# points.shape[0] x 4, i.e. number of points x 4. 
def HDBSClustering(points, min_cluster_size):
    import hdbscan
#     blobs, labels = make_blobs(n_samples=points.shape[0], n_features=2) #centers = 3 is default; number of clusters
#     print 'blobs: ', blobs.size, blobs.shape, labels
#     print points, type(points), points.shape[0], points.shape[1]
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size)
    clusterer.fit(points)
    labels = clusterer.labels_
#     print('labels: ', labels, type(labels), labels.size, labels.shape, labels.shape[0])
#     print('number of labels/clusters: ', clusterer.labels_.max()+1)
    output = np.zeros((len(points), 4))
    intracluster_numbering = {} # dictionary to fill second column of output matrix, i.e. the numbering of a datum within a cluster
    for j in range(-1, clusterer.labels_.max()+1): #include clusterer.labels_.max()
        intracluster_numbering[j] = 0
    for i in range(len(points)):
        output[i]=[labels[i]+1, intracluster_numbering[labels[i]], points[i,1], points[i,0]]
        intracluster_numbering[labels[i]]+=1
    sorted_C = np.array(sorted(output, key=lambda x:x[0]))
    return sorted_C

def DBSClustering(points, eps, min_cluster_size):
    from sklearn.cluster import DBSCAN

    clusterer = DBSCAN(eps=eps,min_samples=min_cluster_size)
    clusterer.fit(points)
    labels = clusterer.labels_
    # print('labels: ', labels, type(labels), labels.size, labels.shape, labels.shape[0])
    # print('number of labels/clusters: ', clusterer.labels_.max()+1)
    output = np.zeros((len(points), 4))
    intracluster_numbering = {} # dictionary to fill second column of output matrix, i.e. the numbering of a datum within a cluster
    for j in range(-1, clusterer.labels_.max()+1): #include clusterer.labels_.max()
        intracluster_numbering[j] = 0
    for i in range(len(points)):
        output[i]=[labels[i]+1, intracluster_numbering[labels[i]], points[i,1], points[i,0]]
        intracluster_numbering[labels[i]]+=1
    sorted_C = np.array(sorted(output, key=lambda x:x[0]))
    return sorted_C


def IDBSClustering(points, iterative, vector_eps, vector_min_samples):
    from iDBSCAN import iDBSCAN

    clusterer = iDBSCAN(iterative = iterative, vector_eps = vector_eps, vector_min_samples = vector_min_samples)
    clusterer.fit(points)
    labels = clusterer.labels_
    # print('labels: ', labels, type(labels), labels.size, labels.shape, labels.shape[0])
    # print('number of labels/clusters: ', clusterer.labels_.max()+1)
    output = np.zeros((len(points), 4))
    intracluster_numbering = {} # dictionary to fill second column of output matrix, i.e. the numbering of a datum within a cluster
    for j in range(-1, clusterer.labels_.max()+1): #include clusterer.labels_.max()
        intracluster_numbering[j] = 0
    for i in range(len(points)):
        output[i]=[labels[i]+1, intracluster_numbering[labels[i]], points[i,1], points[i,0]]
        intracluster_numbering[labels[i]]+=1
    sorted_C = np.array(sorted(output, key=lambda x:x[0]))
    return sorted_C

def clusteringWithMethod(points, min_cluster_size, eps, iterative, vector_eps, vector_min_samples, Cmethod):
    thC = 2         # minimum distanca between points in a cluster (rebinne image)
    if Cmethod == 'hdbs': # nccs
        C = HDBSClustering(points, min_cluster_size)
    elif Cmethod == 'dbsc':
        C = DBSClustering(points, eps, min_cluster_size)
    elif Cmethod == 'idbsc':
        C = IDBSClustering(points, iterative, vector_eps, vector_min_samples)
    elif Cmethod == 'nccs':
        C = NNClustering(points, thC)
    else:
        C = []
    return C

def ClusteringInfo(C):
# ruturn NNClustering clusetr Info array, number of not clusterd, and size of lagre cluster 
    import numpy as np
    NC0  = 0
    NCL  = 0
    maxc = 0
    imax = 0
    info = []
    for i in range(0, len(C)):
        if C[i][1]>0:
            if C[i,1]==1:
                csize = np.where(C[:,0]==C[i,0])[0]
                if len(csize) > maxc:
                    maxc = len(csize)
                    imax = NCL
                info.append(csize)
                NCL +=1
            else:
                NC0 +=1 
    return maxc, imax, info 

# #####################

def PointDistMax(points):
# ruturn the max distance between to point along a line (only a line)
    import numpy as np
    dmax = np.sqrt((points[:,0].max()-points[:,0].min())**2+(points[:,1].max()-points[:,1].min())**2)
    return dmax

def PointDist(p1, p2):
    import numpy as np
    (x1, y1), (x2, y2) = p1, p2
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

# ############################################################################################################################ #
# filling 1d and 2d histograms, credits https://github.com/NichtJens/numpy-accumulative-histograms
# ############################################################################################################################ #

class Hist1D(object):

    def __init__(self, nbins, xlow, xhigh):
        self.nbins = nbins
        self.xlow  = xlow
        self.xhigh = xhigh

        self.range = (xlow, xhigh)

        self.hist, edges = np.histogram([], bins=nbins, range=self.range)
        self.bins = (edges[:-1] + edges[1:]) / 2.

    def fill(self, arr):
        hist, _ = np.histogram(arr, bins=self.nbins, range=self.range)
        self.hist += hist

    @property
    def data(self):
        return self.bins, self.hist


class Hist2D(object):

    def __init__(self, nxbins, xlow, xhigh, nybins, ylow, yhigh):
        self.nxbins = nxbins
        self.xhigh  = xhigh
        self.xlow   = xlow

        self.nybins = nybins
        self.yhigh  = yhigh
        self.ylow   = ylow

        self.nbins  = (nxbins, nybins)
        self.ranges = ((xlow, xhigh), (ylow, yhigh))

        self.hist, xedges, yedges = np.histogram2d([], [], bins=self.nbins, range=self.ranges)
        self.xbins = (xedges[:-1] + xedges[1:]) / 2.
        self.ybins = (yedges[:-1] + yedges[1:]) / 2.

    def fill(self, xarr, yarr):
        hist, _, _ = np.histogram2d(xarr, yarr, bins=self.nbins, range=self.ranges)
        self.hist += hist

    @property
    def data(self):
        return self.xbins, self.ybins, self.hist

# #####################

class fillXY:
    def __init__(self):
        self.n = 0
        self.y = np.array([])
        self.x = np.array([])
    def fill(self, x):
        self.n += 1
        self.x = np.append(self.x, self.n)
        self.y = np.append(self.y, x)
    @property
    def data(self):
        return self.x, self.y
    
    
# ############################################################################################################################ #
# Syile (ATLAS-style)
# ############################################################################################################################ #

def set_atlas_style(shape="medium"):
    import matplotlib.pyplot as plt
    """Set the plotting style to ATLAS-style and then point this function to 'None' so that it can only be called once. Called on canvas creation."""

    # Set figure layout
    if shape == "medium":
        plt.rcParams["figure.figsize"] = (10.0, 6.0)
    elif shape == "large":
        plt.rcParams["figure.figsize"] = (20.0, 20.0)
    elif shape == "xlarge":
        plt.rcParams["figure.figsize"] = (30.0, 30.0)
    elif shape == "long":
        plt.rcParams["figure.figsize"] = (20.0, 5.0)
    elif shape == "xlong":
        plt.rcParams["figure.figsize"] = (30.0, 10.0)
    elif shape == "square":
        plt.rcParams["figure.figsize"] = (6, 6)
    elif shape == "two":
        plt.rcParams['figure.figsize'] = (20.0, 10.0)
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["figure.subplot.bottom"] = 0.16
    plt.rcParams["figure.subplot.top"] = 0.95
    plt.rcParams["figure.subplot.left"] = 0.16
    plt.rcParams["figure.subplot.right"] = 0.95

    # Set font options
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Helvetica, helvetica, Nimbus Sans L, Mukti Narrow, FreeSans"  # alternatives if helvetica is unavailable
    plt.rcParams["font.cursive"] = "Apple Chancery, Textile, Zapf Chancery, Sand, Script MT, Felipa, cursive, Helvetica, helvetica"
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams["mathtext.default"] = "sf"
    plt.rcParams["mathtext.cal"] = "cursive"
    plt.rcParams["mathtext.bf"] = "sans:bold"
    plt.rcParams["mathtext.it"] = "sans:italic"
    plt.rcParams["mathtext.rm"] = "serif"
    plt.rcParams["mathtext.sf"] = "sans"
    plt.rcParams["mathtext.tt"] = "sans"

    # Set axes options
    plt.rcParams["axes.labelsize"] = 20
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["xtick.labelsize"] = 18
    plt.rcParams["xtick.major.size"] = 12
    plt.rcParams["xtick.minor.size"] = 6
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["ytick.labelsize"] = 18
    plt.rcParams["ytick.major.size"] = 14
    plt.rcParams["ytick.minor.size"] = 7

    # Set line options
    plt.rcParams["lines.markersize"] = 8
    plt.rcParams["lines.linewidth"] = 1

    # Set legend options
    plt.rcParams["legend.numpoints"] = 1
    plt.rcParams["legend.fontsize"] = 19
    plt.rcParams["legend.labelspacing"] = 0.3
    plt.rcParams["legend.frameon"] = True
    
    # set title dtyle
    plt.rcParams["axes.titlesize"] = 18
    # Disable calling this function again
    #set_atlas_style.func_code = (lambda: None).func_code