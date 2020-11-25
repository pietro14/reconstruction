#!/usr/bin/env python
import numpy as np
from skimage.measure import LineModelND, ransac
from scipy import spatial
from clusterTools import Cluster

def array_row_intersection(a,b):
   tmp=np.prod(np.swapaxes(a[:,:,None],1,2)==b,axis=2)
   return a[np.sum(np.cumsum(tmp,axis=0)*tmp==1,axis=1).astype(bool)]

class ClusterMatcher:
    def __init__(self,params):
        self.min_length = 100
        self.npixx = 2304
        self.min_pix_intercept = 10
        
    def fitCluster(self,hits,plotting=False):
        #print (hits)
        x = np.array([h[0] for h in hits])
        y = np.array([h[1] for h in hits])
        data = np.column_stack([x, y])
        
        # fit line using all data
        model = LineModelND()
        model.estimate(data)

        # robustly fit line only using inlier data with RANSAC algorithm
        model_robust, inliers = ransac(data, LineModelND, min_samples=4,
                                       residual_threshold=1, max_trials=1000)
        outliers = inliers == False

        line_x = np.arange(0, self.npixx)
        line_y = model.predict_y(line_x)
        line_y_robust = model_robust.predict_y(line_x)

        if plotting:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.plot(data[inliers, 0], data[inliers, 1], '.b', alpha=0.6,
                    label='Inlier data')
            ax.plot(data[outliers, 0], data[outliers, 1], '.r', alpha=0.6,
                    label='Outlier data')
            ax.plot(line_x, line_y, '-k', label='Line model from all data')
            ax.plot(line_x, line_y_robust, '-b', label='Robust line model')
            ax.set_ylim(0,self.npixx)
            ax.legend(loc='lower left')
            plt.show()

        extrap_xy = np.column_stack([line_x, line_y]).astype(int)
        return extrap_xy

    def matchCluster(self,killer,targets):
        if killer.shapes['long_width'] > self.min_length:
            killer_x = np.array([h[0] for h in killer.hits])
            killer_y = np.array([h[1] for h in killer.hits])
            killer_data = np.column_stack([killer_x, killer_y])
            extrap_xy = self.fitCluster(killer.hits)
            #print ("KILLER extrap")
            #print(extrap_xy)
            for i,clu in enumerate(targets):
                hits = clu.hits
                x = np.array([h[0] for h in hits])
                y = np.array([h[1] for h in hits])
                data = np.column_stack([x, y])
                #print ("data cluster = ",i)
                #print (data)
                intersect = array_row_intersection(extrap_xy,data)
                #print ('INTERSECT = ')
                #print (intersect)
                
                if len(intersect)>self.min_pix_intercept:
                    # This solution is optimal when data is very large
                    tree = spatial.cKDTree(data)
                    mindist, minid = tree.query(killer_data)
                    #print("Min dist = ",min(mindist))
                    clu.minDist = min(mindist)

        
if __name__ == '__main__':

    hits = np.load('debug_code/sclu_4.npy')
    fitter = ClusterMatcher({})
    fitter.fitCluster(hits)

    # test the cluster matcher
    killer = Cluster(hits,1,None,None,'lime')
    killer.shapes['long_width'] = 300 # cluster shapes are not calculated when running standalone. Just set to a value to pass the threshold
    targets = []
    for icl in range(7):
        if icl==4:
            continue
        cluhits = np.load('debug_code/sclu_{i}.npy'.format(i=icl))
        cl = Cluster(cluhits,1,None,None,'lime')
        targets.append(cl)
    fitter.matchCluster(killer,targets)

    # check if a cluster has been killed
    for i,cl in enumerate(targets):
        print ("Cluster i = ",i,"has min distance from a killer = ",cl.minDist)
