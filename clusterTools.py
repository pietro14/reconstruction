#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math,itertools
import ROOT
from array import array
from cameraChannel import cameraGeometry

import utilities
utilities = utilities.utils()

class Cluster:
    def __init__(self,hits,rebin):
        self.hits = hits
        self.rebin = rebin
        self.x = hits[:, 0]; self.y = hits[:, 1]
        self.mean_point = np.array([np.mean(self.x),np.mean(self.y)])
        self.EVs = self.eigenvectors()
        self.widths = {}
        self.profiles = {}
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
        def plot_line(center, dir, num_steps=400, step_size=0.5):
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

    def calcProfiles(self,hitscalc=[],plot=None):
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
        # this is in case one wants to make the profile with a different resolution wrt the clustering
        hits = hitscalc if len(hitscalc)>0 else self.hits
        for h in hits:
            rx,ry = rotate_around_point(h,self.EVs[0],self.mean_point)
            rh_major_axis = (rx,ry,h[-1])
            rot_hits.append(rh_major_axis)
        if plot!=None:
            rx = [h[0] for h in rot_hits]; ry = [h[1] for h in rot_hits]; 
            plot.plot(rx, ry, color='green', marker='^',markersize=3)

        # now compute the length along major axis, long profile, etc
        rxmin = min([h[0] for h in rot_hits]); rxmax = max([h[0] for h in rot_hits])
        rymin = min([h[1] for h in rot_hits]); rymax = max([h[1] for h in rot_hits])

        xedg = utilities.dynamicProfileBins(rot_hits,'x',relError=0.3)
        yedg = utilities.dynamicProfileBins(rot_hits,'y',relError=0.3)
        geo = cameraGeometry()
        xedg = [(x-int(rxmin))*geo.pixelwidth for x in xedg]
        yedg = [(y-int(rymin))*geo.pixelwidth for y in yedg]

        length=(rxmax-rxmin)*geo.pixelwidth; width=(rymax-rymin)*geo.pixelwidth
        if len(xedg)>1:
            longprof = ROOT.TProfile('longprof','longitudinal profile',len(xedg)-1,array('f',xedg),'i')
            longprof.SetDirectory(None)
        else: longprof = None
        if len(yedg)>1:
            latprof = ROOT.TProfile('latprof','lateral profile',len(yedg)-1,array('f',yedg),'i')
            latprof.SetDirectory(None)
        else: latprof = None
        
        for h in rot_hits:
            x,y,z=h[0],h[1],h[2]
            if longprof: longprof.Fill((x-rxmin)*geo.pixelwidth,z)
            if latprof: latprof.Fill((y-rymin)*geo.pixelwidth,z)

        profiles = [longprof,latprof]
        for p in profiles:
            if p:
                p.GetXaxis().SetTitle('X (mm)')
                p.GetYaxis().SetTitle('Average photons per bin')
                self.applyProfileStyle(p)
                
        # now set the cluster shapes and profiles
        self.profiles['long'] = longprof
        self.profiles['lat'] = latprof
        self.widths['long'] = length
        self.widths['lat'] = width
        
    def getProfile(self,name='long'):
        if len(self.profiles)==0:
            self.calcProfiles()
        return self.profiles[name] if name in self.profiles else None
    
    def applyProfileStyle(self,prof):
        prof.SetMarkerStyle(ROOT.kFullCircle)
        prof.SetMarkerSize(1)
        prof.SetMarkerColor(ROOT.kBlack)
        prof.SetLineColor(ROOT.kBlack)
        prof.SetMinimum(0)
                
    def hitsFullResolution(self,th2_fullres,pedmap_fullres):
        if hasattr(self,'hits_fr'):
            return self.hits_fr
        else:
            retdict={} # need dict not to duplicate hits after rotation (non integers x,y)
            margin = 30 # in pixels
            for h in self.hits:
                xfull = range(int(h[0]-margin),int(h[0]+margin))
                yfull = range(int(h[1]-margin),int(h[1]+margin))
                fullres = []
                for x in xfull:
                    for y in yfull:
                       xbfull = th2_fullres.GetXaxis().FindBin(x)
                       ybfull = th2_fullres.GetYaxis().FindBin(y)
                       ped = pedmap_fullres.GetBinContent(xbfull,ybfull)
                       z = max(th2_fullres.GetBinContent(xbfull,ybfull)-ped,0)
                       fullres.append((x,y,z))
                for hfr in fullres:
                    x = hfr[0]; y=hfr[1]
                    retdict[(x,y)]=hfr[2]
            ret=[]
            for k,v in retdict.iteritems():
                ret.append((k[0],k[1],v))
            self.hits_fr = np.array(ret)
            return self.hits_fr
    
    def plotFullResolution(self,th2_fullres,pedmap_fullres,name,option='colz'):
        hits_fr = self.hitsFullResolution(th2_fullres,pedmap_fullres)
        border = 30
        xmin,xmax = (min(hits_fr[:,0])-border, max(hits_fr[:,0])+border)
        ymin,ymax = (min(hits_fr[:,1])-border, max(hits_fr[:,1])+border)
        zmax = max(hits_fr[:,2])
        nbinsx = int(xmax-xmin)
        nbinsy = int(ymax-ymin)
        snake_fr = ROOT.TH2D(name,'',nbinsx,xmin,xmax,nbinsy,ymin,ymax)
        for (x,y,zrebin) in hits_fr:
            xb = snake_fr.GetXaxis().FindBin(x)
            yb = snake_fr.GetYaxis().FindBin(y)
            xb_fr = th2_fullres.GetXaxis().FindBin(x)
            yb_fr = th2_fullres.GetYaxis().FindBin(y)
            ped = pedmap_fullres.GetBinContent(xb_fr,yb_fr)
            z = max(th2_fullres.GetBinContent(xb_fr,yb_fr) - ped,0)
            #print "Filling xb,yb =",xb," ",yb," xb_fr,yb_fr",xb_fr," ",yb_fr," with ",x," "," ",y," ",z,"   zrebin = ",zrebin
            snake_fr.SetBinContent(xb,yb,z)
            
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)

        cFR = ROOT.TCanvas("cfr","",600,600)
        snake_fr.GetXaxis().SetTitle('x (pixels)')
        snake_fr.GetYaxis().SetTitle('y (pixels)')
        snake_fr.GetZaxis().SetTitle('counts')
        #snake_fr.GetZaxis().SetRangeUser(0,(zmax*1.05))
        snake_fr.Draw(option)
        for ext in ['png','pdf']:
            cFR.SaveAs('{name}.{ext}'.format(name=name,ext=ext))
        
