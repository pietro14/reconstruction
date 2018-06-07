#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import ROOT

class Cluster:
    def __init__(self,hits,rebin):
        self.hits = hits
        self.x = hits[:, 0]; self.y = hits[:, 1]
        self.mean_point = np.array([np.mean(self.x),np.mean(self.y)])
        self.EVs = self.eigenvectors()
        self.profiles = {}
        self.rebin = rebin
    def integral(self):
        return sum([z for (x,y,z) in self.hits])
    def size(self):
        return len(self.hits)
    def dump(self):
        print self.hits
    def eigenvectors(self):
        covmat = np.cov([self.x,self.y])
        eig_values, eig_vecs = np.linalg.eig(covmat)
        indexes = (np.argmax(eig_values),np.argmin(eig_values))
        eig_vec_vals = (eig_vecs[:, indexes[0]], eig_vecs[:, indexes[-1]])
        return eig_vec_vals

    def plotAxes(self,plot):
        def plot_line(center, dir, num_steps=100, step_size=0.5):
            line_x = []
            line_y = []
            for i in range(num_steps):
                dist_from_center = step_size * (i - num_steps / 2)
                point_on_line = center + dist_from_center * dir
                line_x.append(point_on_line[0])
                line_y.append(point_on_line[1])
            return (line_x, line_y)
        eigen_vectors = self.EVs
        lines = [plot_line(self.mean_point, ev) for ev in eigen_vectors]
        for line in lines:
            plot.plot(line[0], line[1], c="r")

    def calcProfiles(self,plot=None):
        # if they have been attached to the cluster, do not recompute them
        if len(self.profiles)>0:
            return

        def rotate_around_point(hit, dir, pivot):
            x,y = hit[:-1]
            ox, oy = pivot
            cos,sin = dir        
            qx = ox + cos * (x - ox) + sin * (y - oy)
            qy = oy - sin * (x - ox) + cos * (y - oy)
            return qx, qy

        # rotate the hits of the cluster along the major axis
        rot_hits=[]
        for h in self.hits:
            rx,ry = rotate_around_point(h,self.EVs[0],self.mean_point)
            rh_major_axis = (rx,ry,h[-1])
            rot_hits.append(rh_major_axis)
        if plot!=None:
            rx = [h[0] for h in rot_hits]; ry = [h[1] for h in rot_hits]; 
            plot.plot(rx, ry, color='green', marker='^',markersize=3)

        # now comput e the length along major axis, long profile, etc
        rxmin = min([h[0] for h in rot_hits]); rxmax = max([h[0] for h in rot_hits])
        rymin = min([h[1] for h in rot_hits]); rymax = max([h[1] for h in rot_hits])

        bwidth=3
        nbinsx = int((rxmax-rxmin)/float(bwidth))
        nbinsy = int((rymax-rymin)/float(bwidth))
        if nbinsx>1:
            longprof = ROOT.TProfile('longprof','longitudinal profile',nbinsx,0,self.rebin*(rxmax-rxmin))
            longprof.SetDirectory(None)
        else: longprof = None
        if nbinsy>1:
            latprof = ROOT.TProfile('latprof','lateral profile',nbinsy,0,self.rebin*(rymax-rymin))
            latprof.SetDirectory(None)
        else: latprof = None

        for h in rot_hits:
            x,y,z=h[0],h[1],h[2]
            if longprof: longprof.Fill(self.rebin*(x-rxmin),z)
            if latprof: latprof.Fill(self.rebin*(y-rymin),z)

        profiles = [longprof,latprof]
        for p in profiles:
            if p:
                p.GetXaxis().SetTitle('depth (pixels)')
                p.GetYaxis().SetTitle('average counts')
            
        self.profiles['long'] = longprof
        self.profiles['lat'] = latprof
        
    def getProfile(self,name='long'):
        if len(self.profiles)==0:
            self.calcProfiles()
        return self.profiles[name] if name in self.profiles else None
    
