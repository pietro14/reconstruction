## IMPORTING LIBRARIES

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import ROOT
import math

import tools_lib as tl
import sys
sys.path.insert(1, '../')

from iDBSCAN import iDBSCAN
from root_numpy import hist2array
from matplotlib.pyplot import *
from scipy.spatial import distance
from sklearn import metrics
from IPython.display import set_matplotlib_formats
from scipy.ndimage import gaussian_filter, median_filter


def pedsub(img,pedarr):
    return img - pedarr
    
def zsfullres(img_sub,noisearr,nsigma=1):
    img_zs = np.where(img_sub > nsigma * noisearr, img_sub, 0)
    return img_zs

def noisearray(th2):
    noisearr = np.zeros( (th2.GetNbinsX(),th2.GetNbinsY()) )
    for ix in range(th2.GetNbinsX()):
        for iy in range(th2.GetNbinsY()):
            noisearr[ix][iy] = th2.GetBinError(ix+1,iy+1)
    return noisearr

def findedges(ybox,xbox,rescale):
    from skimage.measure import find_contours
    from numpy import zeros
    from scipy.ndimage import uniform_filter
    
    mm = zeros([rescale,rescale],dtype=int)
    mm[ybox,xbox]=10000
    mm = uniform_filter(mm, size=5)
    contours = find_contours(mm, 0.9)
    return contours

def noisereductor(edges,rescale):
    tpx = 10

    for i in range(1,rescale-2):
        for j in range(1,rescale-2):
            spx = edges[i,j]
            mpx = (np.sum(edges[i-1:i+2,j-1:j+2])-spx)/8.
            if np.abs(spx - mpx) > tpx :
                edges[i,j] = mpx
            if (mpx < 0.5):
                edges[i,j] = 0
    return edges

def clusters_neighborood(clusters,data,rescale,radius=3):
    Xtot_clusters = np.vstack(clusters)
    neighboroods = np.zeros([rescale,rescale],dtype=float)
    for pointInCluster in Xtot_clusters:
        ix,iy = pointInCluster[0],pointInCluster[1]
        for ixn in range(-radius,radius+1):
            for iyn in range(-radius,radius+1):
                x = max(0,min(ix+ixn,rescale-1))
                y = max(0,min(iy+iyn,rescale-1))
                neighboroods[x,y] = data[x,y]
    return neighboroods

def store_evolution_in(lst):
    """Returns a callback function to store the evolution of the level sets in
    the given list.
    """

    def _store(x):
        lst.append(np.copy(x))

    return _store

def supercluster(clustered_data):
    print "[Superclustering]"
    from skimage.segmentation import  inverse_gaussian_gradient,checkerboard_level_set, morphological_geodesic_active_contour
    
    gimage = inverse_gaussian_gradient(clustered_data)
    # Initial level set
    # this makes alternate squares active at the first iteration of 10 macro-pixels
    #init_ls = checkerboard_level_set(clustered_data.shape, 10)
    init_ls = np.zeros(clustered_data.shape, dtype=np.int8)
    init_ls[10:-10, 10:-10] = 1
    # List with intermediate results for plotting the evolution
    evolution = []
    callback = store_evolution_in(evolution)
    ls = morphological_geodesic_active_contour(gimage, 230, init_ls,
                                               smoothing=1, balloon=-1,
                                               threshold=0.69,
                                               iter_callback=callback)
    return ls

def supercluster_points(levels):
    from scipy import ndimage as ndi
    # fill the contours of the superclusters
    fill_contours = ndi.binary_fill_holes(levels)
    # remove the smallest superclusters
    label_objects, nb_labels = ndi.label(fill_contours)
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = sizes > 3
    mask_sizes[0] = 0
    contours_cleaned = mask_sizes[label_objects]
    labeled_pixels, nb_labels = ndi.label(contours_cleaned)

    superclusters = [None]*nb_labels

    for ix in range(512):
        for iy in range(512):
            lbl = labeled_pixels[ix,iy]
            if lbl>0: # super-clustered pixel
                if superclusters[lbl-1]: superclusters[lbl-1].append((ix,iy))
                else: superclusters[lbl-1] = [(ix,iy)]
    for isc in range(nb_labels):
        print "Supercluster # ",isc," has these pixels:"
        print superclusters[isc]
                
if __name__ == '__main__':

    ## Setting plotting parameters
    set_matplotlib_formats('png', 'pdf')
    sns.set_context('poster')
    sns.set_style('white')
    sns.set_color_codes()
    plot_kwds = {'alpha' : 0.5, 's' : 30, 'linewidths':0}
    cmapcolor = 'viridis' # or 'viridis'
    vmin      = 99
    vmax      = 125
    figsizeX  = 12 
    figsizeY  = 12
    
    ## Setting Debug Flags
    
    flag_full_image     = 0
    flag_rebin_image    = 0
    flag_edges_image    = 1
    flag_first_it       = 0
    flag_second_it      = 1
    flag_third_it       = 0
    flag_stats          = 1
    flag_all_it         = 1
    flag_plot_noise     = 1

    ## Setting environments variables
    
    rescale     = 512
    cimax       = 200
    nsigma      = 2.3        # numero di sigma sopra il piedistallo
    
    ## Setting i2DBSCAN parameters
    
    iterative     = 4
    tip           = '3D'               # 3D
    
    scale = 4

    if tip == '3D':
       vector_eps         =  [1,   2.5,  5.8,  4]    #[ 3,    5,  7,  9]   #[2.26, 3, 3.5, 4]     #[2, 3, 3.5, 4]      #FOR FNG    #[2.26, 3.5, 2.8, 6]
       vector_min_samples =  [1, 104, 30, 20]         #[30,  200,  100, 100] #[30,  150,  80, 40]   #[30,  55,  28, 13]    #[3,  55,  28, 13]            # [2, 30, 6, 2]
       #vector_eps         = list(np.array(vector_eps, dtype=float)/scale)
       #vector_min_samples = list(np.array(vector_min_samples, dtype=float)/scale)
    else:
       vector_eps         = [1, 2.9, 3.2, 4]
       vector_min_samples = [4,  18,  17, 7]

    cuts = [1500, 400]                # the cut on the length of the track for iteration 1 and 2
    cuts = list(np.array(cuts, dtype=float)/scale)

    ## File folder
    evt               = '00046'#'00016'         # Use always FIVE caracters
    numrun            = '00723'#'00724'         # Use always FIVE caracters
    formattype        = 'h5'            # If the root files comes from h5 conversion use = 'h5'
    filedir           = '../'   # Folder where the root file is placed
    peddir            = '../pedestals/'   # Folder where the pedestal file is placed
    expo              = '100'             # Exposure time of the pedestal in ms

    ## Loading image for analysis
    # otherwise use 'mid'
    if formattype == 'mid':
        imagename     = 'histograms_Run'
        picname       = 'pic_run%s_ev%d' % (numrun,int(evt))
    else:
        imagename     = 'histogram_Run'
        picname       = 'run%d_%s' % (int(numrun),evt)
    
    filename = '%s%s%s.root' % (filedir,imagename,numrun)

    print('Filename: %s' % (filename))
    print('Picname: %s' % (picname))
    
    tf2 = ROOT.TFile.Open(filename)
    imageth2 = tf2.Get(picname)
    image = hist2array(imageth2)
    tf2.Close()

    if flag_full_image == 1:
       fig = plt.figure(figsize=(figsizeX, figsizeY))
       plt.imshow(image,cmap=cmapcolor, vmin=vmin,vmax=vmax,origin='lower' )
       plt.title("Original Image")
       plt.colorbar()
       
       fig = plt.figure(figsize=(figsizeX, figsizeY))
       plt.imshow(image,cmap=cmapcolor, vmin=95,vmax=105,origin='lower' )
       plt.title("Original Image")
       plt.colorbar()

    print('[Image Loaded]')

    # Loading Pedestal files

    tf2ped  = ROOT.TFile.Open('%s/pedmap_ex%s_rebin1.root' % (peddir,expo))  # PEDESTAL MAP
    pedmap  = tf2ped.Get('pedmap').Clone()                               # Getting the pedmap
    pedmap.SetDirectory(None)
    m_image = hist2array(pedmap)
    
    s_image = noisearray(pedmap)
    
    tf2ped.Close()

    # pedrescale     = 2048
    # m_image = np.zeros((pedrescale,pedrescale),dtype=float)
    # s_image = np.zeros((pedrescale,pedrescale),dtype=float)
    # for x in range(1,pedrescale+1):
    #     for y in range(1,pedrescale+1):
    #         m_image[x-1,y-1] = pedmap.GetBinContent(x,y)
    #         s_image[x-1,y-1] = pedmap.GetBinError(x,y)

    #rebin_th_image   = np.round(m_image + nsigma*s_image)       # Threshold Image
    img_cimax = np.where(image < cimax, image, 0)
    img_fr_sub = pedsub(img_cimax,m_image)
    img_fr_zs  = zsfullres(img_fr_sub,s_image,nsigma=nsigma)


    print('[Pedestal Loaded]')

    # Subtracting Pedestal from Image
    rebin_image     = tl.rebin(image, (rescale, rescale))

    edges  = tl.rebin(img_fr_zs,(rescale, rescale))
    edges = noisereductor(edges,rescale)
    edcopy = edges.copy()
    edcopyMedian = median_filter(edcopy, size=4)
       
    points = np.array(np.nonzero(np.round(edcopyMedian))).astype(int).T
    lp = points.shape[0]
    
    
    if flag_rebin_image == 1:
        fig = plt.figure(figsize=(figsizeX, figsizeY))
        plt.imshow(rebin_image, cmap=cmapcolor, vmin=vmin, vmax=vmax, origin='lower' )
        plt.title("Rebinned Image")
    
    if flag_edges_image == 1:
        #     fig = plt.figure(figsize=(figsizeX, figsizeY))
        #     plt.imshow(edges, cmap=cmapcolor, vmin=0, vmax=1, origin='lower' )
        #     plt.title("Edges Image")
        f,ax = plt.subplots(1,2,figsize=(40,20))
        ax[0].imshow(edges, cmap=cmapcolor, vmin=0, vmax=1, origin='lower' )
        ax[0].set_title('Edges after pedestal subtraction')
        ax[1].imshow(edcopyMedian, cmap=cmapcolor, vmin=0, vmax=1, origin='lower' )
        ax[1].set_title('Edges after Filtering')


    print('[Edges Calculated]')

    ## Adding or not the third dimension

    X = points.copy()

    # if tip == '3D':
    #     lp = len(points)
    #     Xl = list(X.copy())
    #     for cor in X:
    #         for count in range(0,np.int(np.round(rebin_image[cor[0],cor[1]] - rebin_th_image[cor[0],cor[1]]) - 1)):
    #             Xl.append(cor)
    #     X = np.array(Xl)
    #     np.save('corImg.npz',X)
    #     print('[3D Method]')
    # else:
    #     print('[2D Method]')


    if tip=='3D':
        Xl = [(ix,iy) for ix,iy in points]          # Aux variable to simulate the Z-dimension
        X1 = np.array(Xl).copy()                    # variable to keep the 2D coordinates
        for ix,iy in points:                        # Looping over the non-empty coordinates
            nreplicas = int(edges[ix,iy])-1
            for count in range(nreplicas):                                # Looping over the number of 'photons' in that coordinate
                Xl.append((ix,iy))                              # add a coordinate repeatedly 
        X = np.array(Xl)                                        # Convert the list to an array
        print('[3D Method]')
    else:
        #     print('[2D Method]')
        X = points.copy()
        #X1 = X

    ## STARTING THE CLUSTERING
    
    ### Defining some environment variables for i2DBSCAN
    #------------------------------------------------------------------
    clusters = iDBSCAN(iterative = iterative, vector_eps = vector_eps, vector_min_samples = vector_min_samples, cuts = cuts, flag_plot_noise = flag_plot_noise ).fit(X)
    #------------------------------------------------------------------ 

    if tip == '3D':
        teste = clusters.tag_
        clusters.labels_ = clusters.labels_[range(0,lp)]
        clusters.tag_ = clusters.tag_[range(0,lp)]
        X = points.copy()
    print('[i2DBSCAN Calculated]')

    u,indices = np.unique(clusters.labels_,return_index = True)
    clu_it1 = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices]==1)[0])].tolist()]
    clu_it2 = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices]==2)[0])].tolist()]
    clu_it3 = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices]==3)[0])].tolist()]
    clu_anyit = clu_it1 + clu_it2
    custers_neighborood_it2 = clusters_neighborood(clu_anyit,edges,rescale)
    superclusters = supercluster(custers_neighborood_it2)
    
    if flag_stats == 1:
        print('[Statistics]')
        n_clusters_ = len(set(clusters.labels_)) - (1 if -1 in clusters.labels_ else 0)
        print("Total number of Clusters: %d" % (n_clusters_))
        u,indices = np.unique(clusters.labels_,return_index = True)
        print("Clusters found in iteration 1: %d" % (sum(clusters.tag_[indices] == 1)))
        print("Clusters found in iteration 2: %d" % (sum(clusters.tag_[indices] == 2)))
        print("Clusters found in iteration 3: %d" % (sum(clusters.tag_[indices] == 3)))

    if flag_first_it == 1:
        print('[Plotting 1st iteration]')
        u,indices = np.unique(clusters.labels_,return_index = True)
        clu = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices] == 1)[0])].tolist()]

        fig = plt.figure(figsize=(figsizeX, figsizeY))
        plt.imshow(rebin_image,cmap=cmapcolor, vmin=vmin,vmax=vmax,origin='lower' )
        plt.title("Clusters found in iteration 1")

        for j in range(0,np.shape(clu)[0]):

            ybox = clu[j][:,0]
            xbox = clu[j][:,1]

            if (len(ybox) > 0) and (len(xbox) > 0):
                contours = findedges(ybox,xbox,rescale)
                for n, contour in enumerate(contours):
                    plt.plot(contour[:, 1],contour[:, 0], '-r',linewidth=2.5)

    if flag_second_it == 1:
        print('[Plotting 2nd iteration]')
        u,indices = np.unique(clusters.labels_,return_index = True)
        clu = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices] == 2)[0])].tolist()]

        fig = plt.figure(figsize=(figsizeX, figsizeY))
        plt.imshow(rebin_image,cmap=cmapcolor, vmin=vmin,vmax=vmax,origin='lower' )
        plt.title("Clusters found in iteration 2")

        for j in range(0,np.shape(clu)[0]):

            ybox = clu[j][:,0]
            xbox = clu[j][:,1]

            if (len(ybox) > 0) and (len(xbox) > 0):
                contours = findedges(ybox,xbox,rescale)
                for n, contour in enumerate(contours):
                    plt.plot(contour[:, 1],contour[:, 0], '-b',linewidth=2.5)

    if flag_third_it == 1:
        print('[Plotting 3rd iteration]')
        u,indices = np.unique(clusters.labels_,return_index = True)
        clu = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices] == 3)[0])].tolist()]

        fig = plt.figure(figsize=(figsizeX, figsizeY))
        plt.imshow(rebin_image,cmap=cmapcolor, vmin=vmin,vmax=vmax,origin='lower' )
        plt.title("Clusters found in iteration 3")

        for j in range(0,np.shape(clu)[0]):

            ybox = clu[j][:,0]
            xbox = clu[j][:,1]

            if (len(ybox) > 0) and (len(xbox) > 0):
                contours = findedges(ybox,xbox,rescale)
                for n, contour in enumerate(contours):
                    plt.plot(contour[:, 1],contour[:, 0], '-y',linewidth=2.5)

    if flag_all_it == 1:
        print('[Plotting ALL iteration]')
        u,indices = np.unique(clusters.labels_,return_index = True)
        clu = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices] == 1)[0])].tolist()]

        fig = plt.figure(figsize=(figsizeX, figsizeY))
        plt.imshow(rebin_image,cmap=cmapcolor, vmin=vmin,vmax=vmax,origin='lower' )
        plt.title("Final Image")

        
        for j in range(0,np.shape(clu)[0]):

            ybox = clu[j][:,0]
            xbox = clu[j][:,1]

            if (len(ybox) > 0) and (len(xbox) > 0):
                contours = findedges(ybox,xbox,rescale)
                for n, contour in enumerate(contours):
                    line, = plt.plot(contour[:, 1],contour[:, 0], '-r',linewidth=2.5)
                if j == 0:
                    line.set_label('1st Iteration')

        clu = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices] == 2)[0])].tolist()]

        for j in range(0,np.shape(clu)[0]):

            ybox = clu[j][:,0]
            xbox = clu[j][:,1]

            if (len(ybox) > 0) and (len(xbox) > 0):
                contours = findedges(ybox,xbox,rescale)
                for n, contour in enumerate(contours):
                    line, = plt.plot(contour[:, 1],contour[:, 0], '-b',linewidth=2.5)
                if j == 0:
                    line.set_label('2nd Iteration')

        clu = [X[clusters.labels_ == i] for i in u[list(np.where(clusters.tag_[indices] == 3)[0])].tolist()]

        for j in range(0,np.shape(clu)[0]):

            ybox = clu[j][:,0]
            xbox = clu[j][:,1]

            if (len(ybox) > 0) and (len(xbox) > 0):
                contours = findedges(ybox,xbox,rescale)
                for n, contour in enumerate(contours):
                    line, = plt.plot(contour[:, 1],contour[:, 0], '-y',linewidth=2.5)
                if j == 0:
                    line.set_label('3rd Iteration')

        supercluster_contour = plt.contour(superclusters, [0.5], colors='g')
        supercluster_contour.collections[0].set_label('supercluster')
        
        plt.legend(loc='upper left')
        #plt.axis([175, 190, 275, 295])
        plt.savefig('allclusters.png', bbox_inches='tight', pad_inches=0)

        
        supercluster_points(superclusters)
