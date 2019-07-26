#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import ROOT,math,os,sys

from scipy.ndimage import gaussian_filter
from skimage import img_as_float
from skimage.morphology import reconstruction
from skimage import measure

from morphsnakes import(morphological_chan_vese,
                        morphological_geodesic_active_contour,
                        inverse_gaussian_gradient,
                        checkerboard_level_set)

from clusterTools import Cluster
from cameraChannel import cameraTools
import matplotlib.pyplot as plt

import debug_code.tools_lib as tl

class SnakesFactory:
    def __init__(self,img,img_fr,name,options):
        self.name = name
        self.options = options
        self.rebin = options.rebin
        ct = cameraTools()
        self.image = img
        self.imagelog = np.zeros((self.image.shape[0],self.image.shape[1]))
        for (x,y),value in np.ndenumerate(self.image):
            if value > 3.0/math.sqrt(self.rebin): # tresholding needed for tracking
                self.imagelog[x,y] = math.log(value)
        self.image_fr = img_fr
        self.contours = []
        
    def store_evolution_in(self,lst):
        """Returns a callback function to store the evolution of the level sets in
        the given list.
        """
        def _store(x):
            lst.append(np.copy(x))
        return _store
    
    def getContours(self,iterations,threshold=0.69):
        
        # Morphological GAC
        image = img_as_float(self.image)
        gimage = inverse_gaussian_gradient(image)

        # Initial level set
        init_ls = np.zeros(image.shape, dtype=np.int8)
        init_ls[10:-10, 10:-10] = 1
        # List with intermediate results for plotting the evolution
        evolution = []
        callback = self.store_evolution_in(evolution)
        ls = morphological_geodesic_active_contour(gimage, iterations, init_ls,
                                                   smoothing=1, balloon=-1,
                                                   threshold=threshold,
                                                   iter_callback=callback)
        # before returning the snakes, put them in the event
        self.contours = ls
        return ls

    def getClusters(self,maxDist=1000,minPoints=20,minPointsCore=10,plot=True):

        from sklearn.cluster import DBSCAN
        from iDBSCAN import iDBSCAN
        from sklearn import metrics
        from scipy.spatial import distance
        
        tip = '3D'
        
        #   IDBSCAN parameters  #
        
        scale              = 1
        iterative          = 4                         # number of iterations for the IDBSC
        if tip == '3D':
            vector_eps         = [2, 2.9, 3.5, 4]          #[2.26, 3.5, 2.8, 6]
            vector_min_samples = [1,  15,  28, 13]            # [2, 30, 6, 2]
        else:
            vector_eps         = [2, 2.9, 3.2, 4]
            vector_min_samples = [2,  18,  17, 7]
        
        vector_eps         = list(np.array(vector_eps, dtype=float)*scale)    
        cuts               = [0, 0]
        
        #-----------------------#
        
        # make the clustering with DBSCAN algo
        # this kills all macrobins with N photons < 1
        points = np.array(np.nonzero(np.round(self.image))).T.astype(int)
        lp = points.shape[0]
        
        Xl = [(ix,iy) for ix,iy in points]          # Aux variable to simulate the Z-dimension
        X1 = np.array(Xl).copy()                    # variable to keep the 2D coordinates
        for ix,iy in points:                        # Looping over the non-empty coordinates
            nreplicas = int(self.image[ix,iy])-1 if tip=='3D' else 1
            for count in range(nreplicas):                                # Looping over the number of 'photons' in that coordinate
                Xl.append((ix,iy))                              # add a coordinate repeatedly 
        X = np.array(Xl)                                        # Convert the list to an array
        
        # - - - - - - - - - - - - - -
        db = iDBSCAN(iterative = iterative, vector_eps = vector_eps, vector_min_samples = vector_min_samples, cuts = cuts).fit(X)
        # Returning to '2' dimensions
        if tip == '3D':
            db.labels_              = db.labels_[range(0,lp)]               # Returning theses variables to the length
            db.tag_                 = db.tag_[range(0,lp)]                  # of the 'real' edges, to exclude the fake repetitions.
        # - - - - - - - - - - - - - -
        
        labels = db.labels_
        
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

        clusters = []
        ##################### plot
        # the following is to preserve the square aspect ratio with all the camera pixels
        # plt.axes().set_aspect('equal','box')
        # plt.ylim(0,2040)
        # plt.xlim(0,2040)

        outname = self.options.plotDir
        if outname and not os.path.exists(outname):
            os.system("mkdir -p "+outname)
            os.system("cp utils/index.php "+outname)

        # Black removed and is used for noise instead.
        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each)
                  for each in np.linspace(0, 1, len(unique_labels))]
        #canv = ROOT.TCanvas('c1','',600,600)
        if plot:
            fig = plt.figure(figsize=(10, 10))
            plt.imshow(self.image.T,cmap='viridis', vmin=1, vmax=10, origin='lower' )
            
        for k, col in zip(unique_labels, colors):
            if k == -1:
                col = [0, 0, 0, 1]
                break # noise: the unclustered

            class_member_mask = (labels == k)
         
            #xy = X[class_member_mask & core_samples_mask]
            xy = X1[class_member_mask]
            
            x = xy[:, 0]; y = xy[:, 1]
                            
            # only add the cores to the clusters saved in the event
            if k>-1 and len(x)>1:
                cl = Cluster(xy,self.rebin,self.image_fr,debug=False)
                cl.iteration = db.tag_[labels == k][0]
                cl.xmax = max(x)
                cl.xmin = min(x)
                cl.ymax = max(y)
                cl.ymin = min(y)
                clusters.append(cl)
                if plot:
                    xri,yri = tl.getContours(x,y)
                    cline = {1:'r',2:'b',3:'y'}
                    plt.plot(xri,yri,'-{lcol}'.format(lcol=cline[cl.iteration]),linewidth=0.5)
                # if plot: cl.plotAxes(plot=plt,num_steps=100)
                # cl.calcProfiles(plot=None)
                # for dir in ['long','lat']:
                #     prof = cl.getProfile(dir)
                #     if prof and cl.widths[dir]>10: # plot the profiles only of sufficiently long snakes
                #         prof.Draw()
                #         for ext in ['png','pdf']:
                #             canv.SaveAs('{pdir}/{name}_snake{iclu}_{dir}profile.{ext}'.format(pdir=outname,name=self.name,iclu=k,dir=dir,ext=ext))

        if plot:
            for ext in ['png','pdf']:
                plt.savefig('{pdir}/{name}.{ext}'.format(pdir=outname,name=self.name,ext=ext))
            plt.gcf().clear()
            plt.close('all')
            
        return clusters

    def getTracks(self,plot=True):
        from skimage.transform import (hough_line, hough_line_peaks)
        # Classic straight-line Hough transform
        image = self.imagelog
        h, theta, d = hough_line(image)
        print("tracks found")
        
        tracks = []
        thr = 0.8 * np.amax(h)
        #######################   IMPLEMENT HERE THE SAVING OF THE TRACKS ############
        # loop over prominent tracks
        itrk = 0
        for _, angle, dist in zip(*hough_line_peaks(h, theta, d,threshold=thr)):
            print("Track # ",itrk)
            #points_along_trk = np.zeros((self.image.shape[1],self.image.shape[0]))
            points_along_trk = []
            for x in range(self.image.shape[1]):
                y = min(self.image.shape[0],max(0,int((dist - x * np.cos(angle)) / np.sin(angle))))
                #points_along_trk[x,y] = self.image[y,x]
                #print "adding point: %d,%d,%f" % (x,y,self.image[y,x])
                # add a halo fo +/- 20 pixels to calculate the lateral profile
                for iy in range(int(y)-5,int(y)+5):
                    if iy<0 or iy>=self.image.shape[0]: continue
                    points_along_trk.append((x,iy,self.image[iy,x]))
            xy = np.array(points_along_trk)
            trk = Cluster(xy,self.rebin)
            tracks.append(trk)
            itrk += 1
        ###################################
            
        if plot:
            # Generating figure
            from matplotlib import cm
            fig, ax = plt.subplots(2, 1, figsize=(18, 6))
            #ax = axes.ravel()

            ax[0].imshow(image, cmap=cm.gray)
            ax[0].set_title('Camera image')
            #ax[0].set_axis_off()            

            ax[1].imshow(image, cmap=cm.gray)
            for _, angle, dist in zip(*hough_line_peaks(h, theta, d,threshold=thr)):
                y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
                y1 = (dist - image.shape[1] * np.cos(angle)) / np.sin(angle)
                ax[1].plot((0, image.shape[1]), (y0, y1), '-r')
            ax[1].set_xlim((0, image.shape[1]))
            ax[1].set_ylim((image.shape[0], 0))
            #ax[1].set_axis_off()
            ax[1].set_title('Fitted tracks')

            plt.tight_layout()
            #plt.show()
            outname = self.options.plotDir
            if outname and not os.path.exists(outname):
                os.system("mkdir -p "+outname)
                os.system("cp ~/cernbox/www/Cygnus/index.php "+outname)
            for ext in ['pdf']:
                plt.savefig('{pdir}/{name}.{ext}'.format(pdir=outname,name=self.name,ext=ext))
            plt.gcf().clear()

        return tracks
        
    def plotClusterFullResolution(self,clusters):
        outname = self.options.plotDir
        for k,cl in enumerate(clusters):
            cl.plotFullResolution('{pdir}/{name}_cluster{iclu}'.format(pdir=outname,name=self.name,iclu=k))

    def calcProfiles(self,clusters):
        for k,cl in enumerate(clusters):
            cl.calcProfiles(plot=None)
                             
    def plotProfiles(self,clusters):
        outname = self.options.plotDir
        canv = ROOT.TCanvas('c1','',1200,600)
        for k,cl in enumerate(clusters):
            for dir in ['long','lat']:
                prof = cl.getProfile(dir)
                if prof and cl.widths[dir]>0.2: # plot the profiles only of sufficiently long snakes (>200 um)
                    prof.Draw("pe1")
                    for ext in ['pdf']:
                        canv.SaveAs('{pdir}/{name}_cluster{iclu}_{dir}profile.{ext}'.format(pdir=outname,name=self.name,iclu=k,dir=dir,ext=ext))
        
    def plotContours(self,contours):

        image = img_as_float(self.image)
        #fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        #ax = axes.flatten()
        fig, ax = plt.subplots()
        
        ax.imshow(image, cmap="gray")
        ax.set_axis_off()
        ax.contour(contours, [0.5], colors='r')
        #ax.set_title("Morphological GAC segmentation", fontsize=12)
        (run,event) = self.name.split('_')
        ax.set_title('Run={run}, Event={event}'.format(run=run,event=event), fontsize=12)
        
        # ax[1].imshow(contours, cmap="gray")
        # ax[1].set_axis_off()
        # contour = ax[1].contour(evolution[0], [0.5], colors='g')
        # contour.collections[0].set_label("Iteration 0")
        # contour = ax[1].contour(evolution[50], [0.5], colors='y')
        # contour.collections[0].set_label("Iteration 50")
        # contour = ax[1].contour(evolution[-1], [0.5], colors='r')
        # contour.collections[0].set_label("Iteration 100")
        # ax[1].legend(loc="upper right")
        # title = "Morphological GAC evolution"
        # ax[1].set_title(title, fontsize=12)
        
        fig.tight_layout()
        #plt.show()
        for ext in ['pdf']:
            plt.savefig('{name}.{ext}'.format(name=self.name,ext=ext))



class SnakesProducer:
    def __init__(self,sources,params,options):
        self.picture   = sources['picture']   if 'picture' in sources else None
        self.pictureHD = sources['pictureHD'] if 'pictureHD' in sources else None
        self.name      = sources['name']      if 'name' in sources else None
        self.algo      = sources['algo']      if 'algo' in sources else 'DBSCAN'
        
        self.snakeQualityLevel = params['snake_qual']   if 'snake_qual' in params else 3
        self.plot2D            = params['plot2D']       if 'plot2D' in params else False
        self.plotpy            = params['plotpy']       if 'plotpy' in params else False
        self.plotprofiles      = params['plotprofiles'] if 'plotprofiles' in params else False

        self.options = options

    def run(self):
        ret = []
        if any([x==None for x in (self.picture.any(),self.pictureHD.any(),self.name)]):
            return ret
        
        # Cluster reconstruction on 2D picture
        snfac = SnakesFactory(self.picture,self.pictureHD,self.name,self.options)

        # this plotting is only the pyplot representation.
        # Doesn't work on MacOS with multithreading for some reason... 
        if self.algo=='DBSCAN':
            snakes = snfac.getClusters(plot=self.plotpy)
        elif self.algo=='HOUGH':
            snakes = snfac.getTracks(plot=self.plotpy)            

        # print "Get light profiles..."
        snfac.calcProfiles(snakes)
        
        # snfac.calcProfiles(snakes) # this is for BTF
        
        # sort snakes by light integral
        snakes = sorted(snakes, key = lambda x: x.integral(), reverse=True)
        # and reject discharges (round)
        #snakes = [x for x in snakes if x.qualityLevel()>=self.snakeQualityLevel]
        
        # plotting
        if self.plot2D:       snfac.plotClusterFullResolution(snakes)
        if self.plotprofiles: snfac.plotProfiles(snakes)

        return snakes
