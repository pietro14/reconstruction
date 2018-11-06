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
        self.shapes = {}
        
    def integral(self):
        return sum([z for (x,y,z) in self.hits])

    def getSize(self,name='long'):
        if len(self.profiles)==0:
            self.calcProfiles()
        if name in self.widths: return self.widths[name]
        else:
            print "ERROR! You can only get 'long' or 'lat' sizes!"
            return -999

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

        # rotate the hits of the cluster along the major axis
        rot_hits=[]
        # this is in case one wants to make the profile with a different resolution wrt the clustering
        hits = hitscalc if len(hitscalc)>0 else self.hits
        for h in hits:
            rx,ry = utilities.rotate_around_point(h,self.EVs[0],self.mean_point)
            rh_major_axis = (rx,ry,h[-1])
            rot_hits.append(rh_major_axis)
        if plot!=None:
            rx = [h[0] for h in rot_hits]; ry = [h[1] for h in rot_hits]; 
            plot.plot(rx, ry, color='green', marker='^',markersize=3)

        # now compute the length along major axis, long profile, etc
        rxmin = min([h[0] for h in rot_hits]); rxmax = max([h[0] for h in rot_hits])
        rymin = min([h[1] for h in rot_hits]); rymax = max([h[1] for h in rot_hits])
        xedg = utilities.dynamicProfileBins(rot_hits,'x',relError=0.1)
        yedg = utilities.dynamicProfileBins(rot_hits,'y',relError=0.5)
        # rxmin = 1000; rxmax = 1500 # central 5 cm of the track
        # xedg = range(int(rxmin),int(rxmax)) # full binning for BTF tracks
        # yedg = range(int(rymin),int(rymax)) # full binning for BTF tracks
        geo = cameraGeometry()
        xedg = [(x-int(rxmin))*geo.pixelwidth for x in xedg]
        yedg = [(y-int(rymin))*geo.pixelwidth for y in yedg]

        length=(rxmax-rxmin)*geo.pixelwidth; width=(rymax-rymin)*geo.pixelwidth
        if len(xedg)>1:
            longprof = ROOT.TH1F('longprof','longitudinal profile',len(xedg)-1,array('f',xedg),'i')
            longprof.SetDirectory(None)
        else: longprof = None
        if len(yedg)>1:
            latprof = ROOT.TH1F('latprof','lateral profile',len(yedg)-1,array('f',yedg),'i')
            latprof.SetDirectory(None)
        else: latprof = None
        
        for h in rot_hits:
            x,y,z=h[0],h[1],h[2]
            if longprof: longprof.Fill((x-rxmin)*geo.pixelwidth,z)
            if latprof: latprof.Fill((y-rymin)*geo.pixelwidth,z)

        profiles = [longprof,latprof]
        titles = ['longitudinal','transverse']
        for ip,p in enumerate(profiles):
            if p:
                p.GetXaxis().SetTitle('X_{%s} (mm)' % titles[ip])
                p.GetYaxis().SetTitle('Average photons per bin')
                self.applyProfileStyle(p)
                
        # now set the cluster shapes and profiles
        self.profiles['long'] = longprof
        self.profiles['lat'] = latprof
        # those are not used, since they include the "margins" at 0
        # just used as the starting values in clusterShapes()
        self.widths['long'] = length
        self.widths['lat'] = width

        # get the peaks inside the profile
        for direction in ['lat','long']:
            self.clusterShapes(direction)
        
    def getProfile(self,name='long'):
        if len(self.profiles)==0:
            self.calcProfiles()
        return self.profiles[name] if name in self.profiles else None

    def clusterShapes(self,name='long'):
        # ensure the cluster profiles are ready
        if name not in ['lat','long']:
            print "ERROR! Requested profile along the ",name," direction. Should be either 'long' or 'lat'. Exiting clusterShapes()."
            return
        self.getProfile(name)

        from waveform import PeakFinder,simplePeak

        # find first the length/width with intersection of the base of the large peak
        # threshold = 3
        # min_distance_peaks = 5 # number of bins of the profile, to be converted in mm later... TO DO
        # prominence = 2 # noise seems <1
        # width = 10  # find only 1 big peak
        # pf = PeakFinder(self.profiles[name])        
        # pf.findPeaks(threshold,min_distance_peaks,prominence,width)
        # self.widths[name] = pf.getFWHMs()[0] if len(pf.getFWHMs()) else 0 # first should be the only big peak
        self.shapes['%s_width' % name] = self.widths[name]
        
        # find the peaks and store their properties
        # thresholds on the light. Should be configurable...
        threshold = 3
        min_distance_peaks = 3 # number of bins of the profile, to be converted in mm later... TO DO
        prominence = 2 # noise seems <1
        width = 1 # minimal width of the signal
        pf = PeakFinder(self.profiles[name])        
        pf.findPeaks(threshold,min_distance_peaks,prominence,width)

        amplitudes = pf.getAmplitudes()
        prominences = pf.getProminences()
        fwhms = pf.getFWHMs()
        peakPositions = pf.getPeakTimes()
        
        peaksInProfile = [simplePeak(amplitudes[i],prominences[i],peakPositions[i],fwhms[i]) for i in xrange(len(amplitudes))]
        peaksInProfile = sorted(peaksInProfile, key = lambda x: x.mean, reverse=True)

        self.shapes[name+'_fullrms']          = self.profiles[name].GetRMS()
        if len(peaksInProfile):
            mainPeak = peaksInProfile[0]
            self.shapes[name+'_p0amplitude']  = mainPeak.amplitude
            self.shapes[name+'_p0prominence'] = mainPeak.prominence
            self.shapes[name+'_p0mean']       = mainPeak.mean
            self.shapes[name+'_p0fwhm']       = mainPeak.fwhm
        else:
            self.shapes[name+'_p0amplitude']  = -999
            self.shapes[name+'_p0prominence'] = -999
            self.shapes[name+'_p0mean']       = -999
            self.shapes[name+'_p0fwhm']       = -999
            
    def applyProfileStyle(self,prof):
        prof.SetMarkerStyle(ROOT.kFullCircle)
        prof.SetMarkerSize(0.5)
        prof.SetMarkerColor(ROOT.kBlack)
        prof.SetLineColor(ROOT.kGray)
        prof.SetLineWidth(1)
        prof.SetMinimum(-1.0)
        
    def hitsFullResolution(self,th2_fullres,pedmap_fullres,zs=False):
        if hasattr(self,'hits_fr'):
            return self.hits_fr
        else:
            retdict={} # need dict not to duplicate hits after rotation (non integers x,y)
            #latmargin = 15 # in pixels
            #longmargin = 100 # in pixels
            latmargin = 30 # in pixels
            longmargin = 5 # in pixels
            for h in self.hits:
                rx,ry = utilities.rotate_around_point(h,self.EVs[0],self.mean_point)
                rxfull = range(int(rx-longmargin),int(rx+longmargin))
                ryfull = range(int(ry-latmargin),int(ry+latmargin))
                fullres = []
                for rxf in rxfull:
                    for ryf in ryfull:
                        rhf = (rxf,ryf,-1)
                        xfull,yfull = utilities.rotate_around_point(rhf,self.EVs[0],self.mean_point,inverse=True)
                        # these for are to ensure that one includes all bins after rounding/rotation
                        for xfullint in range(int(xfull-1),int(xfull+1)):
                            for yfullint in range(int(yfull-1),int(yfull+1)):
                                xbfull = th2_fullres.GetXaxis().FindBin(xfullint)
                                ybfull = th2_fullres.GetYaxis().FindBin(yfullint)
                                ped = pedmap_fullres.GetBinContent(xbfull,ybfull)
                                noise = pedmap_fullres.GetBinError(xbfull,ybfull)
                                z = th2_fullres.GetBinContent(xbfull,ybfull)-ped
                                if zs and z<0.5*noise:
                                    continue
                                fullres.append((xfullint,yfullint,z))
                for hfr in fullres:
                    x = hfr[0]; y=hfr[1]
                    retdict[(x,y)]=hfr[2]
            ret=[]
            for k,v in retdict.iteritems():
                ret.append((k[0],k[1],v))
            self.hits_fr = np.array(ret)
            return self.hits_fr
    
    def plotFullResolution(self,th2_fullres,pedmap_fullres,name,option='colz'):
        # REPETITION!!!! IMPROVE...
        if self.rebin==1:
            hits_fr = self.hits
        else:
            hits_fr = self.hitsFullResolution(th2_fullres,pedmap_fullres)
        border = 15
        xmin,xmax = (min(hits_fr[:,0])-border, max(hits_fr[:,0])+border)
        ymin,ymax = (min(hits_fr[:,1])-border, max(hits_fr[:,1])+border)
        zmax = max(hits_fr[:,2])
        nbinsx = int(xmax-xmin)
        nbinsy = int(ymax-ymin)
        snake_fr = ROOT.TH2D(name,'',nbinsx,xmin,xmax,nbinsy,ymin,ymax)
        for (x,y,z) in hits_fr:
            xb = snake_fr.GetXaxis().FindBin(x)
            yb = snake_fr.GetYaxis().FindBin(y)
            snake_fr.SetBinContent(xb,yb,z)
            
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)

        cFR = ROOT.TCanvas("cfr","",600,600)
        snake_fr.GetXaxis().SetTitle('x (pixels)')
        snake_fr.GetYaxis().SetTitle('y (pixels)')
        snake_fr.GetZaxis().SetTitle('counts')
        # just for the 2D plotting, cut at 1.5 (mean of the RMS of all the pixels)
        #snake_fr.GetZaxis().SetRangeUser(.0,(zmax*1.05))
        snake_fr.GetZaxis().SetRangeUser(0,10)
        snake_fr.Draw(option)
        print "cluster integral = ",snake_fr.Integral()
        #cFR.SetRightMargin(0.2); cFR.SetLeftMargin(0.1); cFR.SetBottomMargin(0.1);
        cFR.SetBottomMargin(0.3); cFR.SetLeftMargin(0.2); cFR.SetRightMargin(0.2); 
        for ext in ['png','pdf']:
            cFR.SaveAs('{name}.{ext}'.format(name=name,ext=ext))


    def qualityLevel(self):
        # result: 1=loose, 2=medium, 3=tight
        kGood = 0
        # sanity (they are not sparse points clustered)
        if self.shapes['lat_p0fwhm']>-999: kGood += 1
        # sphericity
        if self.shapes['lat_width']>0 and self.shapes['long_width']/self.shapes['lat_width']>2.0: kGood += 1
        # uniformity in long/lateral shape
        if self.shapes['lat_fullrms']>0 and self.shapes['long_fullrms']/self.shapes['lat_fullrms']>2.0: kGood += 1
        return kGood
