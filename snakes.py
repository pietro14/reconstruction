#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import ROOT,math,os

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

class SnakesFactory:
    def __init__(self,th2,name,options):
        self.name = name
        self.options = options
        self.rebin = options.rebin
        ct = cameraTools()
        self.data = ct.getData(th2)
        self.X = ct.getActiveCoords(th2)
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
        image = img_as_float(self.data)
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

    def getClusters(self,maxDist=500,minPoints=6,minPointsCore=3,plot=True):

        from sklearn.cluster import DBSCAN
        from sklearn import metrics
        from scipy.spatial import distance

        # make the clustering with DBSCAN algo
        X = self.X
        distance_matrix = distance.squareform(distance.pdist(X))
        db = DBSCAN(eps=maxDist, min_samples=minPoints,metric='euclidean',n_jobs=-1).fit(distance_matrix)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
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
            os.system("cp ~/cernbox/www/Cygnus/index.php "+outname)

        # Black removed and is used for noise instead.
        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each)
                  for each in np.linspace(0, 1, len(unique_labels))]
        #canv = ROOT.TCanvas('c1','',600,600)
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = [0, 0, 0, 1]

            class_member_mask = (labels == k)
         
            xy = X[class_member_mask & core_samples_mask]
            x = xy[:, 0]; y = xy[:, 1]
            if plot:
                plt.plot(x, y, 'o', markerfacecolor=tuple(col),
                         markeredgecolor='k', markersize=10)

            # only add the cores to the clusters saved in the event
            if k>-1 and len(xy)>minPointsCore:
                # print "Found cluster!"
                cl = Cluster(xy,self.rebin)
                clusters.append(cl)
                if plot: cl.plotAxes(plot=plt)
                # cl.calcProfiles(plot=None)
                # for dir in ['long','lat']:
                #     prof = cl.getProfile(dir)
                #     if prof and cl.widths[dir]>10: # plot the profiles only of sufficiently long snakes
                #         prof.Draw()
                #         for ext in ['png','pdf']:
                #             canv.SaveAs('{pdir}/{name}_snake{iclu}_{dir}profile.{ext}'.format(pdir=outname,name=self.name,iclu=k,dir=dir,ext=ext))

            # plot also the non-core hits
            # xy = X[class_member_mask & ~core_samples_mask]
            # if plot: plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
            #                   markeredgecolor='k', markersize=6)

        if plot:
            plt.title('Estimated number of snakes: %d' % n_clusters_)
            #plt.show()
            for ext in ['png','pdf']:
                plt.savefig('{pdir}/{name}.{ext}'.format(pdir=outname,name=self.name,ext=ext))
            plt.gcf().clear()

        return clusters

    def getTracks(self,plot=True):
        from skimage.transform import (hough_line, hough_line_peaks)
        # Classic straight-line Hough transform
        image = self.data
        h, theta, d = hough_line(image)

        if plot:
            # Generating figure
            from matplotlib import cm
            fig, axes = plt.subplots(1, 3, figsize=(15, 6))
            ax = axes.ravel()

            ax[0].imshow(image, cmap=cm.gray)
            ax[0].set_title('Input image')
            ax[0].set_axis_off()
            
            ax[1].imshow(np.log(1 + h),
                         extent=[np.rad2deg(theta[-1]), np.rad2deg(theta[0]), d[-1], d[0]],
                         cmap=cm.gray, aspect=1/1.5)
            ax[1].set_title('Hough transform')
            ax[1].set_xlabel('Angles (degrees)')
            ax[1].set_ylabel('Distance (pixels)')
            ax[1].axis('image')
            
            ax[2].imshow(image, cmap=cm.gray)
            thr = 0.7 * np.amax(h)
            for _, angle, dist in zip(*hough_line_peaks(h, theta, d,threshold=thr)):
                y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
                y1 = (dist - image.shape[1] * np.cos(angle)) / np.sin(angle)
                ax[2].plot((0, image.shape[1]), (y0, y1), '-r')
            ax[2].set_xlim((0, image.shape[1]))
            ax[2].set_ylim((image.shape[0], 0))
            ax[2].set_axis_off()
            ax[2].set_title('Detected lines')
            
            plt.tight_layout()
            plt.show()
            
        
    def plotClusterFullResolution(self,clusters,th2_fullres,pedmap_fullres):
        outname = self.options.plotDir
        for k,cl in enumerate(clusters):
            cl.plotFullResolution(th2_fullres,pedmap_fullres,'{pdir}/{name}_cluster{iclu}'.format(pdir=outname,name=self.name,iclu=k))

    def calcProfiles(self,clusters,th2_fullres=None,pedmap_fullres=None):
        for k,cl in enumerate(clusters):
            hits = cl.hitsFullResolution(th2_fullres,pedmap_fullres) if (th2_fullres and pedmap_fullres) else cl.hits
            cl.calcProfiles(hitscalc=hits,plot=None)
                             
    def plotProfiles(self,clusters):
        outname = self.options.plotDir
        canv = ROOT.TCanvas('c1','',600,600)
        for k,cl in enumerate(clusters):
            for dir in ['long','lat']:
                prof = cl.getProfile(dir)
                if prof and cl.widths[dir]>0.2: # plot the profiles only of sufficiently long snakes (>200 um)
                    prof.Draw()
                    line = ROOT.TLine(prof.GetXaxis().GetXmin(),0,prof.GetXaxis().GetXmax(),0)
                    line.SetLineWidth(2); line.SetLineColor(ROOT.kGray); line.SetLineStyle(ROOT.kDashed)
                    line.Draw("L")
                    for ext in ['png','pdf']:
                        canv.SaveAs('{pdir}/{name}_snake{iclu}_{dir}profile.{ext}'.format(pdir=outname,name=self.name,iclu=k,dir=dir,ext=ext))
        
    def plotContours(self,contours):

        image = img_as_float(self.data)
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
        for ext in ['png','pdf']:
            plt.savefig('{name}.{ext}'.format(name=self.name,ext=ext))



class SnakesProducer:
    def __init__(self,sources,params,options):
        self.picture   = sources['picture']   if 'picture' in sources else None
        self.pictureHD = sources['pictureHD'] if 'pictureHD' in sources else None
        self.pedmapHD  = sources['pedmapHD']  if 'pedmapHD' in sources else None
        self.name      = sources['name']      if 'name' in sources else None
        self.algo      = sources['algo']      if 'algo' in sources else 'DBSCAN'
        
        self.snakeQualityLevel = params['snake_qual']   if 'snake_qual' in params else 3
        self.plot2D            = params['plot2D']       if 'plot2D' in params else True
        self.plotpy            = params['plotpy']       if 'plotpy' in params else False
        self.plotprofiles      = params['plotprofiles'] if 'plotprofiles' in params else True

        self.options = options

    def run(self):
        ret = []
        if any([x==None for x in self.picture,self.pictureHD,self.pedmapHD,self.name]):
            return ret
        
        # Cluster reconstruction on 2D picture
        snfac = SnakesFactory(self.picture,self.name,self.options)

        # this plotting is only the pyplot representation.
        # Doesn't work on MacOS with multithreading for some reason... 
        if self.algo=='DBSCAN':
            snakes = snfac.getClusters(plot=self.plotpy)
        elif self.algo=='HOUGH':
            snakes = snfac.getTracks(plot=self.plotpy)            
            exit(0)
        snfac.calcProfiles(snakes,self.pictureHD,self.pedmapHD)

        # sort snakes by longitudinal width
        snakes = sorted(snakes, key = lambda x: x.widths['long'], reverse=True)
        # and reject discharges (round)
        snakes = filter(lambda x: x.qualityLevel()>=self.snakeQualityLevel, snakes)
        
        # plotting
        if self.plot2D:       snfac.plotClusterFullResolution(snakes,self.pictureHD,self.pedmapHD)
        if self.plotprofiles: snfac.plotProfiles(snakes)

        return snakes
