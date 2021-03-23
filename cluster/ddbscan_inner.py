from __future__ import division

import numpy as np
from sklearn.linear_model import RANSACRegressor
from sklearn.metrics import mean_squared_error
from operator import itemgetter
import time,math

import warnings
warnings.simplefilter('ignore', np.RankWarning)

class PolynomialRegression(object):
    def __init__(self, degree=3, coeffs=None):
        self.degree = degree
        self.coeffs = coeffs

    def fit(self, X, y):
        self.coeffs = np.polyfit(X.ravel(), y, self.degree)

    def get_params(self, deep=False):
        return {'coeffs': self.coeffs}

    def set_params(self, coeffs=None, random_state=None):
        self.coeffs = coeffs

    def predict(self, X):
        poly_eqn = np.poly1d(self.coeffs)
        y_hat = poly_eqn(X.ravel())
        return y_hat

    def score(self, X, y):
        return mean_squared_error(y, self.predict(X))

def ransac_polyfit(x,y,order,t,n=0.8,k=100,f=0.9):

    print("\t\t*** doing polyfit with order ",order)
    besterr = np.inf
    bestfit = np.array([None])
    bestfitderi = np.array([None])
    for kk in range(k):
        maybeinliers = np.random.randint(len(x), size=int(n*len(x)))
        maybemodel = np.polyfit(x[maybeinliers], y[maybeinliers], order)
        polyderi = []
        for i in range(order):
            polyderi.append(maybemodel[i]*(order-i))
        res_th = t / np.cos(np.arctan(np.polyval(np.array(polyderi),x)))

        alsoinliers = np.abs(np.polyval(maybemodel, x)-y) < res_th
        if sum(alsoinliers) > len(x)*f:
            bettermodel = np.polyfit(x[alsoinliers], y[alsoinliers], order)
            polyderi = []
            for i in range(order):
                polyderi.append(maybemodel[i]*(order-i))
            thiserr = np.sum(np.abs(np.polyval(bettermodel, x[alsoinliers])-y[alsoinliers]))
            if thiserr < besterr:
                bestfit = bettermodel
                besterr = thiserr
                bestfitderi = np.array(polyderi)

    print("\t\t=== polyfit DONE")                     
    return bestfit, bestfitderi

#Parameters of the new ransac function:
# x, y - x and y coordinates 
# order - Order of the polynomial
# t - Thickness of the track
# n - Random fraction of the data used to find the polyfit
# k - Number of tries 
# f - Accuracy of the RANSAC to consider the fit a good one

def ddbscaninner(data, is_core, neighborhoods, neighborhoods2, labels, min_samples, dir_radius, dir_min_accuracy, dir_minsamples, dir_thickness, time_threshold, max_attempts, isolation_radius):
    #Definitions
    
    #Beginning of the algorithm - DBSCAN check part
    label_num = 0
    stack = []
    clu_stra = []
    acc = []
    length = []
    clu_labels = []
    
    print("Start DBSCAN seeding...")
    #Loop 
    for i in range(labels.shape[0]):
        if labels[i] != -1 or not is_core[i]:
            continue
        while True:
            if labels[i] == -1:
                labels[i] = label_num
                if is_core[i]:     #Only core points are expanded
                    neighb = neighborhoods[i]
                    for i in range(neighb.shape[0]):
                        v = neighb[i]
                        if labels[v] == -1:
                            stack.append(v)
            if len(stack) == 0:
                break
            i = stack[len(stack)-1]
            del(stack[len(stack)-1])


        #Ransac part
        print("test cluster n. ",i)
        if sum(labels==label_num) > dir_minsamples:
            print("==> cluster ",i," has ",sum(labels==label_num)," samples")
            x = data[labels==label_num][:,0]
            y = data[labels==label_num][:,1]
            ransac = RANSACRegressor(PolynomialRegression(degree=3),
                                     residual_threshold=3 * np.std(y),
                                     random_state=0)
            ransac.fit(np.expand_dims(x, axis=1), y)
            
            accuracy = sum(ransac.inlier_mask_)/len(y)

            if accuracy > dir_min_accuracy:
                clu_stra.append(label_num)
                acc.append(accuracy)
                length.append(sum(labels==label_num))
                clu_labels.append(label_num)
        label_num += 1
    
    #End of DBSCAN loop - check if directional part is viable
    print("Clusters found in DBSCAN: %d" %(len(set(labels)) - (1 if -1 in labels else 0)))

    if len(clu_stra) == 0:
        #If no cluster has a good fit model, the output will be the same of the DBSCAN
        la_aux = np.copy(labels)
        labels = np.zeros([la_aux.shape[0],2], dtype=np.intp)
        labels[:,0] = la_aux
        return labels
    else:
        #If any cluster has a good fit model, it'll be marked from the worst fitted cluster to the best, each of them respecting the accuracy threshold
        auxiliar_points = []
        vet_aux = np.zeros([len(clu_stra),4])
        vet_aux[:,0] = np.asarray(clu_stra)
        vet_aux[:,1] = np.asarray(acc)
        vet_aux[:,2] = np.asarray(length)
        vet_aux = np.asarray(sorted(vet_aux,key=itemgetter(1),reverse=1))
        # in case there is more than 1 cluster with 0 energy around, sort them by accuracy
        if (sum(vet_aux[:,1]==1) > 1):
            l1 = sum(vet_aux[:,1]==1)
            vet_aux[0:l1,:] = np.asarray(sorted(vet_aux[0:l1,:],key=itemgetter(2),reverse=1))
        for u in range(len(clu_stra)):
            lt = (labels==vet_aux[u][0])*is_core
            auxiliar_points.append(np.where(lt)[0][0])
            print("The point %d has been assigned as part of a good fit" %(np.where(lt)[0][0]))
        
        #Now the clusterization will begin from zero with directionality enabled for the clusters that have a good fit model
        label_num = 0
        labels = np.full(data.shape[0], -1, dtype=np.intp)
        stack = []
        for i in auxiliar_points:
            print("Auxiliar point ",i)
            if labels[i] != -1 or not is_core[i]:
                continue
            while True:
                if labels[i] == -1:
                    labels[i] = label_num
                    if is_core[i]:
                        neighb = neighborhoods[i]
                        for i in range(neighb.shape[0]):
                            v = neighb[i]
                            if labels[v] == -1:
                                stack.append(v)

                if len(stack) == 0:
                    break
                i = stack[len(stack)-1]
                del(stack[len(stack)-1])

            print("An attempt cluster fround with ",sum(labels==label_num)," requested ",dir_minsamples, " min samples.  Dir search now...")
            #Now that the provisional cluster has been found, directional search begins
            if sum(labels==label_num) > dir_minsamples:
                #Taking unique points to use on the ransac
                clu_coordinates = [tuple(row) for row in data[labels==label_num]] 
                uniques = np.unique(clu_coordinates,axis=0)
                x = uniques[:,0]
                y = uniques[:,1]
                center_i = (np.average(np.unique(x)),np.average(np.unique(y)))
                print ("cluster with center ",center_i," has samples > dir_min: = ",sum(labels==label_num))
                
                #RANSAC fit
                fit_model, fit_deri = ransac_polyfit(x,y,order=1, t = dir_thickness)
                counter = 1
                #Adding new points to the cluster (If the fit_model output is None, then no model was found)
                if sum(fit_model == None) == 0:
                    control = 1
                    pts1 = 0
                    t1 = time.time()
                    while True:
                        print ("***** attempt n ",counter)
                        #Filling stack list with possible new points to be added (start point)
                        pts0 = pts1
                        moment_lab = np.where(labels==label_num)[0]
                        stack = []
                        for j in moment_lab:
                            neig2 = neighborhoods2[j]
                            for k in neig2:
                                if labels[k] != label_num and labels[k]!=-1:
                                #if (labels[k] != label_num) or 1:
                                    stack.append(k)
                        stack = np.unique(stack).tolist()
                        if len(stack) == 0:
                            break
                        
                        res_th = dir_thickness / np.cos(np.arctan(np.polyval(fit_deri,data[:,0])))
                        inliers_bool = np.abs(np.polyval(fit_model, data[:,0])-data[:,1]) < res_th
                        inliers = np.where(inliers_bool)[0]
                            

                        i = stack[len(stack)-1]
                        #Adding the inliers points from stack list and filling stack with more possible points
                        while True:
                            if i in inliers and (labels[i] != label_num):
                            #if i in inliers:
                                labels[i] = label_num
                                if is_core[i]:
                                    neig2 = neighborhoods2[i]
                                    for i in range(neig2.shape[0]):
                                        v = neig2[i]
                                        if labels[v] != label_num:
                                            stack.append(v)

                            if len(stack) == 0:
                                break
                            i = stack[len(stack)-1]
                            del(stack[len(stack)-1])

                        #Checking current cluster for possible fit model update

                        clu_coordinates = [tuple(row) for row in data[labels==label_num]] 
                        uniques = np.unique(clu_coordinates,axis=0)
                        x = uniques[:,0]
                        y = uniques[:,1]


                        #Updating the ransac model
                        if control == 1:
                            fit_model, fit_deri = ransac_polyfit(x,y,order=1, t = dir_thickness)
                        else:
                            fit_model, fit_deri = ransac_polyfit(x,y,order=3, t = dir_thickness)
                        pts1 = sum(labels==label_num)
                        #Stop criteria - time
                        t2 = time.time()
                        print ("INFO: fit time = ",t2-t1)
                        if (t2 - t1) > time_threshold:
                            break
                        #Stop criteria - max attempts of improvements
                        if counter > max_attempts:
                            break
                        #Stop criteria - When there is no more point to be added or if the fit is not good anymore
                        counter = counter + 1
                        if (pts1 == pts0) or (sum(fit_model == None) != 0):
                            print ("last control = ",control)
                            if control == 0:
                                print('The cluster %d' %(label_num) + ' needed %d attempts' %(counter))
                                break
                            else:
                                fit_model, fit_deri = ransac_polyfit(x,y,order=3, t = dir_thickness)
                                control = 0
                                if sum(fit_model == None) != 0:
                                    break
                
                
            #label_num += 1
            if sum(labels==label_num) > min_samples:
                label_num += 1
            else:
                labels[labels==label_num] = len(data)
            
        # Now that the clusters with good fit models were found, the rest of the data will be clustered with the standard DBSCAN logic
        clustered_data = np.unique(data[labels!=-1],axis=0)
        for i in range(labels.shape[0]):
            if labels[i] != -1 or not is_core[i]:
                continue
            while True:
                if labels[i] == -1:                    
                    labels[i] = label_num
                    if is_core[i]:     #Only core points are expanded
                        neighb = neighborhoods[i]
                        for i in range(neighb.shape[0]):
                            v = neighb[i]
                            isolated = True
                            for cd in clustered_data:
                                dist = math.dist((cd[0],cd[1]),(data[v][0],data[v][1]))
                                if dist<isolation_radius:
                                    #print ("excluding point ",data[v], "close to ",cd," with dist = ",dist)
                                    isolated = False
                                    break
                            if labels[v] == -1 and isolated:
                                stack.append(v)

                if len(stack) == 0:
                    break
                i = stack[len(stack)-1]
                del(stack[len(stack)-1])
        
            #label_num += 1
            if sum(labels==label_num) > min_samples:
                label_num += 1
            else:
                labels[labels==label_num] = len(data)
        
        #False clusters remotion
        labels[labels==len(data)] = -1
        
        #Clusterization has finished, now the clusters found with ransac fit model will be marked
        la_aux = np.copy(labels)
        labels = np.zeros([la_aux.shape[0],2], dtype=np.intp)
        labels[:,0] = la_aux
        labels[auxiliar_points,1] = 1


        return labels
        
            
            
        
        
        
        
        
            
