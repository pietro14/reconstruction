#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import ROOT

from scipy.ndimage import gaussian_filter
from skimage import img_as_float
from skimage.morphology import reconstruction

from morphsnakes import(morphological_chan_vese,
                        morphological_geodesic_active_contour,
                        inverse_gaussian_gradient,
                        checkerboard_level_set)
from skimage import measure
import matplotlib.pyplot as plt

class Cluster:
    def __init__(self,hits):
        self.hits = hits
    def integral(self):
        return sum([z for (x,y,z) in self.hits])
    def size(self):
        return len(self.hits)
    def dump(self):
        print self.hits
    
class SnakesFactory:
    def __init__(self,th2,name):
        self.data = self.getData(th2)
        self.X = self.getActiveCoords(th2)
        self.contours = []
        self.name = name
        
    def store_evolution_in(self,lst):
        """Returns a callback function to store the evolution of the level sets in
        the given list.
        """
        def _store(x):
            lst.append(np.copy(x))
        return _store

    def getData(self,th2):
        if not th2.InheritsFrom("TH2"):
            print "ERROR! The input object should be a TH2"
            return []
        x_bins = th2.GetNbinsX()
        y_bins = th2.GetNbinsY()
        bins = np.zeros((x_bins,y_bins))
        for y_bin in xrange(y_bins): 
            for x_bin in xrange(x_bins): 
                z = th2.GetBinContent(x_bin + 1,y_bin + 1)
                if z>0:
                    bins[y_bin,x_bin] = th2.GetBinContent(x_bin + 1,y_bin + 1)
        return bins

    def getActiveCoords(self,th2):
        ret = []
        if not th2.InheritsFrom("TH2"):
            print "ERROR! The input object should be a TH2"
            return ret
        x_bins = th2.GetNbinsX()
        y_bins = th2.GetNbinsY()
        for y_bin in xrange(y_bins): 
            for x_bin in xrange(x_bins): 
                z = th2.GetBinContent(x_bin + 1,y_bin + 1)
                if z>0:
                    ret.append((x_bin,y_bin,z))
        return np.array(ret)
    
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

    def plot_axis(self,center, dir, num_steps=100, step_size=0.5):
        line_x = []
        line_y = []
        for i in range(num_steps):
            dist_from_center = step_size * (i - num_steps / 2)
            point_on_line = center + dist_from_center * dir
            line_x.append(point_on_line[0])
            line_y.append(point_on_line[1])
        return (line_x, line_y)

    def getClusters(self,maxDist=100,minPoints=10,minPointsCore=3,plot=True):

        from sklearn.cluster import DBSCAN
        from sklearn import metrics
        from scipy.spatial import distance

        # make the clustering with DBSCAN algo
        X = self.X
        distance_matrix = distance.squareform(distance.pdist(X))
        db = DBSCAN(eps=maxDist, min_samples=minPoints,n_jobs=-1).fit(distance_matrix)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

        clusters = []
        ##################### plot
        import matplotlib.pyplot as plt

        # Black removed and is used for noise instead.
        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each)
                  for each in np.linspace(0, 1, len(unique_labels))]
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = [0, 0, 0, 1]

            class_member_mask = (labels == k)
         
            xy = X[class_member_mask & core_samples_mask]
            x = xy[:, 0]; y = xy[:, 1]
            if plot: plt.plot(x, y, 'o', markerfacecolor=tuple(col),
                              markeredgecolor='k', markersize=14)

            # only add the cores to the clusters saved in the event
            if len(xy)>minPointsCore:
                cl = Cluster(xy)
                clusters.append(cl)

                print "cluster size = ",len(xy)
                covmat = np.cov([x,y])
                print covmat

                eig_values, eig_vecs = np.linalg.eig(covmat)
                indexes = (np.argmax(eig_values),np.argmin(eig_values))
                eig_vec_vals = (eig_vecs[:, indexes[0]], eig_vecs[:, indexes[-1]])
                print "largest,smallest eignevalues = ",eig_vec_vals

                if plot:
                    mean_point=np.array([np.mean(x),np.mean(y)])
                    lines = [self.plot_axis(mean_point, ev) for ev in eig_vec_vals] 
                    for line in lines:
                        plt.plot(line[0], line[1], c="r")

                
            # plot also the non-core hits
            xy = X[class_member_mask & ~core_samples_mask]
            if plot: plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                              markeredgecolor='k', markersize=6)

        if plot:
            plt.title('Estimated number of clusters: %d' % n_clusters_)
            #plt.show()
            for ext in ['png','pdf']:
                plt.savefig('{name}.{ext}'.format(name=self.name,ext=ext))
            plt.gcf().clear()

            
        return clusters
        
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
            
