#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon April  20 11:15:28 2021

@author: atulprajapati
"""
import numpy as np
import root_numpy as rn
import ROOT
import cv2
import optparse
import os
import time
import sys
import math
import matplotlib.pyplot as plt
from skimage.morphology import skeletonize, thin, medial_axis,remove_small_holes,remove_small_objects
from scipy.ndimage import  median_filter
import mahotas as mh
from sklearn.decomposition import PCA
from statistics import stdev
from scipy import interpolate
from scipy.spatial.distance import pdist, squareform
from scipy.signal import find_peaks
from scipy import interpolate
from scipy.stats import skew

###############  Definitions ####################
def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops=dict(arrowstyle='->',
                    linewidth=2,
                    shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)
    
def getClusterMatrix(sc_xpixelcoord,sc_ypixelcoord,sc_zpixel):
    df = [sc_xpixelcoord,sc_ypixelcoord,sc_zpixel]
    hits = [None]*len(df[0])
    for ij in range(len(df[0])):
        hits[ij] = [x[ij] for x in df]
        
    xmin = int(min(df[0])); xmax = int(max(df[0]))
    ymin = int(min(df[1])); ymax = int(max(df[1]))
    
    data = np.zeros((int(xmax-xmin)+1,int(ymax-ymin)+1), dtype='float')
    
    for x,y,z in hits:
#            data[int(x-xmin-1),int(y-ymin-1)] = z    #filling the track with actual intensity of the pixel
        data[int(x-xmin),int(y-ymin)] = 255           # filling the track with white pixels just for the computation of skeleton, this avoids holes in the track
    return data
    
def endPoints(skel):
    endpoint1=np.array([[0, 0, 0],
                        [0, 1, 0],
                        [2, 1, 2]])
        
    endpoint2=np.array([[0, 0, 0],
                        [0, 1, 2],
                        [0, 2, 1]])
        
    endpoint3=np.array([[0, 0, 2],
                        [0, 1, 1],
                        [0, 0, 2]])
        
    endpoint4=np.array([[0, 2, 1],
                        [0, 1, 2],
                        [0, 0, 0]])
        
    endpoint5=np.array([[2, 1, 2],
                        [0, 1, 0],
                        [0, 0, 0]])
        
    endpoint6=np.array([[1, 2, 0],
                        [2, 1, 0],
                        [0, 0, 0]])
        
    endpoint7=np.array([[2, 0, 0],
                        [1, 1, 0],
                        [2, 0, 0]])
        
    endpoint8=np.array([[0, 0, 0],
                        [2, 1, 0],
                        [1, 2, 0]])
        
    ep1=mh.morph.hitmiss(skel,endpoint1)
    ep2=mh.morph.hitmiss(skel,endpoint2)
    ep3=mh.morph.hitmiss(skel,endpoint3)
    ep4=mh.morph.hitmiss(skel,endpoint4)
    ep5=mh.morph.hitmiss(skel,endpoint5)
    ep6=mh.morph.hitmiss(skel,endpoint6)
    ep7=mh.morph.hitmiss(skel,endpoint7)
    ep8=mh.morph.hitmiss(skel,endpoint8)
    ep = ep1+ep2+ep3+ep4+ep5+ep6+ep7+ep8
    return ep

def pruning(skeleton, size):   #remove iteratively end points "size" times from the skeleton
    for ip in range(size):
        endpoints = endPoints(skeleton)
        endpoints = np.logical_not(endpoints)
        skeleton = np.logical_and(skeleton,endpoints)
    return skeleton
    
def makeitthin(xpixels,ypixels,zpixels): #this finds the skeleton of the track by 2 methods, one by skeletonize and other by thinning
    
    img = getClusterMatrix(xpixels,ypixels,zpixels) # gives the img of the track
    _,th_skel = cv2.threshold(img.astype('uint8'),100,255,cv2.THRESH_BINARY) # returns the thresholded track
    
    blurred1 = median_filter(th_skel,6) #applying a median filter with kernel size of 6 pixels
    arr = blurred1 > 0                 # converting the img to boolean
    cleaned = remove_small_objects(arr, min_size=7)        # This removes the noisy pixels around the track
    cleaned = remove_small_holes(cleaned,7)                 # this removes the holes inside the track
    _,cln_img = cv2.threshold(cleaned.astype('uint8'),0, 255,cv2.THRESH_BINARY) # this returns the clean image

    img2 = cv2.cvtColor(cln_img, cv2.COLOR_GRAY2BGR) # converting it to 3 channel img for skeletonization
    blurred_skel = cv2.blur(img2,(5,5))   # bluring the image after cleaning to smoothen the edges
    blurred_thin = cv2.blur(cln_img,(5,5)) # bluring the image after cleaning to smoothen the edges

    skel_img = skeletonize(blurred_skel) #skeleton of the track
    thinned = thin(blurred_thin)        # skeleton by thinning
    gray_skel =cv2.cvtColor(skel_img, cv2.COLOR_BGR2GRAY) # converting RGB img to gray scale
    cts_skel = cv2.countNonZero(gray_skel)  # Counting the number of non zero pixels
    thin_gray = thinned.astype('uint8')
    cts_thin = cv2.countNonZero(thin_gray)
    _,thresh_thin = cv2.threshold(thin_gray,0, 255,cv2.THRESH_BINARY)
    
    skel_track_len = cts_skel
    thin_track_len = cts_thin
    
    if thin_track_len > 50:
        pruned_skel = pruning(gray_skel,8)          # Pruning if the track is longer than 50 pixels to remove the small branches
        pruned_thin = pruning(thin_gray,8)
        _,pruned_skel = cv2.threshold(pruned_skel.astype('uint8'),0, 255,cv2.THRESH_BINARY)
        _,pruned_thin = cv2.threshold(pruned_thin.astype('uint8'),0, 255,cv2.THRESH_BINARY)
        thresh_thin = pruned_thin
        cts_skel1 = cv2.countNonZero(pruned_skel)
        cts_thin1 = cv2.countNonZero(pruned_thin)
        skel_track = cts_skel1
        thin_track = cts_thin1
    else:
        skel_track = cts_skel          # Track length using skeletonization
        thin_track = cts_thin          # Track length using thinning
        
    ##### To find the min area rectangle ###
    (x_track,y_track) = np.nonzero(thresh_thin)
    arr_track = np.array([y_track,x_track]).T
    rect = cv2.minAreaRect(arr_track)
    box = cv2.boxPoints(rect)
    box = np.int0(box)
    (_, (width, height), angle_full) = cv2.minAreaRect(arr_track)

    length_full_rect = np.max([width,height])         # Length of the full track by using the major axis of the rectangle
    angle_full_track = angle_full                                                                       # Direction of the full track
    
    ##### To find the directionality of the track ####
    if len(arr_track) < 10:
        (_,_,angle) = cv2.minAreaRect(arr_track)
    else:
        (_,_,angle) = cv2.minAreaRect(arr_track[0:10,:])
    angle_10_pixels = angle                                              # Direction of the begining of track
    curlyness = length_full_rect/thin_track                        # Curlyness

    return skel_track ,thin_track,length_full_rect,angle_full_track,angle_10_pixels,curlyness,thresh_thin,img,img2,box
        

#######################################################################################################################################
if __name__ == "__main__":

    #CALL: python reco_track.py -I <<INPUT_FOLDER>> -E <<entries>> -P <<Plot>> -T <<textfile>> -O <<outfolder>>

    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")

    parser.add_option("-I", "--inputfile", dest="inputfile", help="specify the name of the input files")
    parser.add_option("-E", "--entries", dest="entries", default = 1,type = 'int', help="specify the number of entries")
    parser.add_option("-P", "--pt", dest="pt", default = 0, type = 'int', help="specify whether to plot or not")
    parser.add_option("-O", "--outputfolder", dest="outfolder", default=os.getcwd()+"/out", help="specify the output destination folder")
    parser.add_option("-T", "--textfile", dest="textfile", default = 0,type = 'int', help="to save the variables in text file")

  
    (opt, args) = parser.parse_args()
        
########################################################################################################################################



#### CODE EXECUTION ####

    t0=time.time()
        
    if not os.path.exists(opt.outfolder): #CREATING OUTPUT FOLDER
        os.makedirs(opt.outfolder)
            
    rootfile=ROOT.TFile.Open(opt.inputfile)
    tree=rootfile.Get('Events')            #GETTING NTUPLES
            
    infilename=opt.inputfile[:-5]
    outfile=ROOT.TFile(opt.outfolder+'/'+infilename+'.root', 'RECREATE')    # Creating an ouput root file
    
    
    ############### Defining the variables ################
    param_tree = ROOT.TTree("param_tree","param_tree")                      #creating a tree
    
    skel_track=np.empty((1),dtype="float32")
    thin_track=np.empty((1),dtype="float32")
    SDCD = np.empty((1),dtype="float32")
    CylThick = np.empty((1),dtype="float32")
    LAPA = np.empty((1),dtype="float32")
    MaxDen = np.empty((1),dtype="float32")
    eta = np.empty((1),dtype="float32")
    ChargeUnif = np.empty((1),dtype="float32")
    SC_size = np.empty((1),dtype="float32")
    SC_nhits = np.empty((1),dtype="float32")
    SC_integral = np.empty((1),dtype="float32")
    SC_length = np.empty((1),dtype="float32")
    SC_pathlength = np.empty((1),dtype="float32")
    SC_width = np.empty((1),dtype="float32")
    SC_xmean = np.empty((1),dtype="float32")
    SC_ymean = np.empty((1),dtype="float32")
    curlyness = np.empty((1),dtype="float32")
    angle_full_track = np.empty((1),dtype="float32")
    angle_10_pixels = np.empty((1),dtype="float32")
    length_full_rect = np.empty((1),dtype="float32")
    num_peaks = np.empty((1),dtype="float32")
    skew_track = np.empty((1),dtype="float32")
    skew_peaks = np.empty((1),dtype="float32")

    ######## Defing the branches of the tree ############
    param_tree.Branch('skel_track',skel_track,"skel_track/F")
    param_tree.Branch('thin_track',thin_track,"thin_track/F")
    param_tree.Branch('SDCD',SDCD,"SDCD/F")
    param_tree.Branch('CylThick',CylThick,"CylThick/F")
    param_tree.Branch('LAPA',LAPA,"LAPA/F")
    param_tree.Branch('MaxDen',MaxDen,"MaxDen/F")
    param_tree.Branch('eta',eta,"eta/F")
    param_tree.Branch('ChargeUnif',ChargeUnif,"ChargeUnif/F")
    param_tree.Branch('SC_size',SC_size,"SC_size/F")
    param_tree.Branch('SC_nhits',SC_nhits,"SC_nhits/F")
    param_tree.Branch('SC_integral',SC_integral,"SC_integral/F")
    param_tree.Branch('SC_length',SC_length,"SC_length/F")
    param_tree.Branch('SC_pathlength',SC_pathlength,"SC_pathlength/F")
    param_tree.Branch('SC_width',SC_width,"SC_width/F")
    param_tree.Branch('SC_xmean',SC_xmean,"SC_xmean/F")
    param_tree.Branch('SC_ymean',SC_ymean,"SC_ymean/F")
    param_tree.Branch('curlyness',curlyness,"curlyness/F")
    param_tree.Branch('angle_full_track',angle_full_track,"angle_full_track/F")
    param_tree.Branch('angle_10_pixels',angle_10_pixels,"angle_10_pixels/F")
    param_tree.Branch('length_full_rect',length_full_rect,"length_full_rect/F")
    param_tree.Branch('num_peaks',num_peaks,"num_peaks/F")
    param_tree.Branch('skew_track',skew_track,"skew_track/F")
    param_tree.Branch('skew_peaks',skew_peaks,"skew_peaks/F")
    ################################################################################
    start = 0
    if opt.entries==-1:
        totev=tree.GetEntries()
    else:
                    
        if opt.entries<=tree.GetEntries():
            totev=opt.entries
        else:
            totev=tree.GetEntries()
    if opt.pt == 1:
        start = opt.entries
        totev = start+1
    ################################################
    if opt.textfile == 1:
        txt_file = np.zeros((totev,23))     # Creating an output txt file
        vars_names = ['skel_track','thin_track','SDCD','CylThick','ChargeUnif','LAPA','MaxDen','eta','length_full_rect','angle_full_track','angle_10_pixels','curlyness','SC_size','SC_nhits','SC_integral','SC_length','SC_width','SC_pathlength','SC_xmean','SC_ymean','num_peaks','skew_track','skew_peaks']
    #######################################################
    for entry in range(start, totev): #RUNNING ON ENTRIES
        tree.GetEntry(entry)
        unique_scid = np.unique(tree.sc_ID)     # to select the number of super clusters in the given image
        scid = unique_scid[unique_scid >= 0]    # to remove the super clusters with scID == -1
        nsc = 0
        for id in scid:
            ##### thresholded pixels ####
            thresh_pixels = np.array([tree.sc_ID,tree.sc_xpixelcoord,tree.sc_ypixelcoord,tree.sc_zpixel])
            th_xpixels = thresh_pixels[1,thresh_pixels[0,:]==id].T  # selecting the coordinates of a given supercluster
            th_ypixels = thresh_pixels[2,thresh_pixels[0,:]==id].T
            th_zpixels = thresh_pixels[3,thresh_pixels[0,:]==id].T
            ####### all pixels ###########
            all_pixels = np.array([tree.sc_IDall,tree.sc_xallpixelcoord,tree.sc_yallpixelcoord,tree.sc_zallpixel])
            xpixels = all_pixels[1,all_pixels[0,:]==id].T
            ypixels = all_pixels[2,all_pixels[0,:]==id].T
            zpixels = all_pixels[3,all_pixels[0,:]==id].T
            
            SC_size[0] = tree.sc_size[nsc]                      # No. of pixels in a SC
            SC_nhits[0] = tree.sc_nhits[nsc]                    # No. of pixels above threshold in SC
            SC_integral[0] = tree.sc_integral[nsc]              # Integral of the SC
            SC_length[0] = tree.sc_length[nsc]                  # Major axis of the ellipse fitting the SC
            SC_pathlength[0] = tree.sc_pathlength[nsc]          # No. of slices* Radius of slice
            SC_width[0] = tree.sc_width[nsc]                    # Minor axis of the ellipse fitting the SC
            SC_xmean[0] = tree.sc_xmean[nsc]                    # Mean of the x coordinate of the SC
            SC_ymean[0] = tree.sc_ymean[nsc]                    # mean of the y coordinate of the SC

            ################# Skeletonization ######################
            skel_track[0], thin_track[0] ,length_full_rect[0],angle_full_track[0],angle_10_pixels[0],curlyness[0],thresh_thin,img,img2,box= makeitthin(xpixels,ypixels,zpixels)

            ############## PCA Analysis #############
            signal = thresh_pixels[1:3,thresh_pixels[0,:]==id].T # Track coordinates with given track id
            pca = PCA(n_components=2) # initializing 2 component PCA
            pca.fit(signal) # fitting the track to find the principal axis

            start_point = pca.mean_ - 5 * pca.components_[0]  # this is to draw the principal axis
            end_point = pca.mean_ + 5 * pca.components_[0]  # this is to draw the principal axis
            x = [start_point[0],end_point[0]]               # this is to draw the principal axis
            y = [start_point[1],end_point[1]]               # this is to draw the principal axis

            slope = np.empty(2,dtype = float) # empty vector to fill the slope of pricipal axis and the slope of axis perpendicular to Principal axis
            c = np.empty(2,dtype = float)  # this is the intercept
            for length, vector in zip(pca.explained_variance_, pca.components_):
                i = 0
                v = vector * 3 * np.sqrt(length)
                if v[0]==0 or v[1]==0:
                    slope[i]=0
                    c[i]= pca.mean_[1]
                else:
                    slope[i] = v[1]/v[0]
                    c[i]= pca.mean_[1]-(pca.mean_[0]*slope[i])
                i += 1
            ############# Discriminating Variables ###############
            X = signal
            x_min = abs(min(X.dot(pca.components_[0]))) # min projected value
            x_max = abs(max(X.dot(pca.components_[0]))) # max projected value
            
            LAPA[0] = abs(x_max-x_min)                                                              # Length along principal axis
            ###################### Number of peaks #################
            peak_hist = X.dot(pca.components_[0]);     # projection of track on the principal axis
            plot1 = plt.figure(1)
            (bin_val, num_bins, patches) = plt.hist(peak_hist, bins = int(LAPA[0]))  # to find the counts in each bin of the histogram
            plt.xlabel('Bins')
            plt.ylabel('Counts')
            tck, u= interpolate.splprep([num_bins[0:-1],bin_val], k=3,s=500)       # interpolating the bin counts to find a smooth curve (this helps in removing the sharp and small peaks which are generally fluctuations)
            interpolated_points = interpolate.splev(np.linspace(0,1,len(bin_val)), tck, der=0) # finding the interpolated points

            peaks, properties = find_peaks(interpolated_points[1],distance=6) # number of peaks in the insterpolated data (distance =10 means it ignores the peaks within 10pixels of 1 peak)
            num_peaks[0] = len(peaks)
#            print('Number of peaks is %d.'%(len(peaks)))
            skew_track[0] = abs(skew(interpolated_points[1]))  # here I compute the skewness of the track (NR tracks at lower energies are generally symmetrical)
            skew_peaks[0] = num_peaks[0] + skew_track[0]
            
            #################################################################
            data_mean = np.mean(X,axis =0)
            SDCD[0] = np.sqrt(np.sum((X[:,0]-data_mean[0])**2 +(X[:,1]-data_mean[1])**2)/len(X))    # Standard Deviation of Charge Distribution
            dist = abs(slope[0]*X[:,0]-X[:,1]+c[0])/np.sqrt(slope[0]**2 + 1)
            CylThick[0] = np.sum(dist**2)/len(X)                                                    # Cylindrical Thickness

            min_x = np.min(xpixels)
            max_x = np.max(xpixels)
            min_y = np.min(ypixels)
            max_y = np.max(ypixels)
            image_y = ROOT.TH2F('Original track'+str(entry+nsc),'original track'+str(entry+nsc),int(max_x-min_x),min_x,max_x,int(max_y-min_y),min_y,max_y)
            for n in range(len(th_xpixels)):
                image_y.Fill(int(th_xpixels[n]),int(th_ypixels[n]),th_zpixels[n])
        
            image_y.Rebin2D(2)  #Rebinning the image to find the maximum density exclusing the noisy pixels by rebinning
            MaxDen[0] = image_y.GetMaximum()                                                        # Maximum Density
            eta[0] = MaxDen[0]/thin_track[0]                                                        # Eta
            ############
            xbins = image_y.GetXaxis().GetNbins()
            ybins = image_y.GetYaxis().GetNbins()
            x_hist = np.zeros((xbins*ybins))
            y_hist = np.zeros((xbins*ybins))
            z_hist = np.zeros((xbins*ybins))
            n = 0
            for i in range(xbins):
                for j in range(ybins):
                    x_hist[n] = i
                    y_hist[n] = j
                    z_hist[n] = image_y.GetBinContent(i,j)
                    n +=1
            Xarr = np.array([x_hist,y_hist,z_hist])

            ################################
            Xarr = Xarr.T
            Xarr = Xarr[Xarr[:,2] != 0,:]
            ##### ChargeUnif ####
#                dist_vec = pdist(np.array([th_xpixels,th_ypixels]).T,'euclidean')  # This is distance computation without rebinning the track
            dist_vec = pdist(np.array([Xarr[:,0],Xarr[:,1]]).T,'euclidean')   # This is distance computation after rebinning the track
            dist_mat = squareform(dist_vec,'tomatrix')      # this converts distance matrix to sqaure form, from this we have information of distance of each point from all other points
            avg_dist = dist_mat.sum(axis = 0)/(len(Xarr)-1) #This is to find the  average distance of each point from all the other points (use len(th_xpixels) if computing the dist vec with th_pixels
            ChargeUnif[0] = stdev(avg_dist)                      # Charge Uniformity
            ##### Filling txt file ######
            if opt.textfile == 1:
                vars_mat = np.array([skel_track[0],thin_track[0],SDCD[0],CylThick[0],ChargeUnif[0],LAPA[0],MaxDen[0],eta[0],length_full_rect[0],angle_full_track[0],angle_10_pixels[0],curlyness[0],SC_size[0],SC_nhits[0],SC_integral[0],SC_length[0],SC_width[0],SC_pathlength[0],SC_xmean[0],SC_ymean[0],num_peaks[0],skew_track[0],skew_peaks[0]])
                txt_file[entry,:] = vars_mat

            #### Plotting the tracks and features #####
            if opt.pt == 1:
                print('Track length with thinning is %f pixels.'%(thin_track[0]))
                print('SDCD is %f, CylThick is %f, LAPA is %f and ChargeUnif is %f, MaxDen is %f , Eta is %f.'%(SDCD[0],CylThick[0],LAPA[0],ChargeUnif[0],MaxDen[0],eta[0]))
                cv2.imshow('Thinned Track',thresh_thin)
                #cv2.imwrite('Thinned_track_'+str(nsc)+'.png',thresh_thin)
                cv2.waitKey(-1)
                cv2.imshow('Original Track',img)
                cv2.waitKey(-1)
                cv2.imshow('Cleaned track',img2)
                cv2.waitKey(-1)
                cv2.drawContours(thresh_thin,[box],0,(255,0,0),0);
                cv2.imshow('Minimum Area Rectangle',thresh_thin)
                cv2.waitKey(-1)
                plot2 = plt.figure(2)
                for length, vector in zip(pca.explained_variance_, pca.components_):
                    v = vector * 2 * np.sqrt(length)
                    draw_vector(pca.mean_, pca.mean_ + v)
                plt.axis('equal');
                plt.scatter(signal[:, 0], signal[:, 1], alpha=0.2,label='original data')
                plt.plot(x,y,label='first principal axis',c='r')
                plt.legend()
                plt.show()
                
                y = X.dot(pca.components_[1])
                plot3 = plt.figure(3)
                plt.scatter(X.dot(pca.components_[0]), y, alpha=.3,label='Projected data onto first PCA component')
    #           plt.set(xlabel='Projected data onto first PCA component', ylabel='y')
                plt.legend()
                plt.show()
                
                plot4 = plt.figure(4)
                plt.plot(interpolated_points[1],label='Projected data onto first PCA component')
                plt.plot(peaks,interpolated_points[1][peaks],"s",label='peaks')
                plt.xlabel('Bins')
                plt.ylabel('Counts')
                plt.legend()
                plt.show()
            nsc += 1
            param_tree.Fill()
        print('Entry number %d.'%(entry))
    if opt.textfile == 1:
        np.savetxt(opt.outfolder+'/'+infilename+'.csv',txt_file,delimiter = ' ', header = ', '.join(vars_names))
    t1=time.time()
    print('Task Completed in %f s.'%(t1-t0))
    outfile.Write()
    outfile.Close()
                
