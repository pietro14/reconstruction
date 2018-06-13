#!/usr/bin/env python

import os,math
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.SetBatch(True)
from channelsTools import cameraChannel

from snakes import SnakesFactory

class analysis:

    def __init__(self,rfile,options):
        self.xmax = 2048
        self.rebin = options.rebin        
        self.rfile = rfile
        self.options = options
        self.pedfile_name = '{base}_ped_rebin{rb}.root'.format(base=os.path.splitext(self.rfile)[0],rb=self.rebin)
        if not os.path.exists(self.pedfile_name):
            print "WARNING: pedestal file ",self.pedfile_name, " not existing. First calculate them..."
            self.calcPedestal(options.numPedEvents)
        self.pedfile_fullres_name = '{base}_ped_rebin1.root'.format(base=os.path.splitext(self.rfile)[0])
        if not os.path.exists(self.pedfile_fullres_name):
            print "WARNING: pedestal file with full resolution ",self.pedfile_fullres_name, " not existing. First calculate them..."
            self.calcPedestal(options.numPedEvents,1)
        print "Pulling pedestals..."
        # first the one for clustering with rebin
        pedrf = ROOT.TFile.Open(self.pedfile_name)
        self.pedmap = pedrf.Get('pedmap').Clone()
        self.pedmap.SetDirectory(None)
        self.pedmean = pedrf.Get('pedmean').GetMean()
        self.pedrms = pedrf.Get('pedrms').GetMean()
        pedrf.Close()
        # then the full resolution one
        pedrf_fr = ROOT.TFile.Open(self.pedfile_fullres_name)
        self.pedmap_fr = pedrf_fr.Get('pedmap').Clone()
        self.pedmap_fr.SetDirectory(None)
        pedrf_fr.Close()
        
    def calcPedestal(self,maxImages=-1,alternativeRebin=-1):
        nx=ny=self.xmax
        rebin = self.rebin if alternativeRebin<0 else alternativeRebin
        nx=int(nx/rebin); ny=int(ny/rebin); 
        pedfile = ROOT.TFile.Open(self.pedfile_name,'recreate')
        pedmap = ROOT.TProfile2D('pedmap','pedmap',nx,0,self.xmax,ny,0,self.xmax,'s')
        tf = ROOT.TFile.Open(self.rfile)
        for i,e in enumerate(tf.GetListOfKeys()):
            if maxImages>-1 and i<len(tf.GetListOfKeys())-maxImages: continue
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print "Calc pedestal with event: ",name
            obj.RebinX(rebin); obj.RebinY(rebin); 
            for ix in xrange(nx+1):
                for iy in xrange(ny+1):
                    x = obj.GetXaxis().GetBinCenter(ix+1)
                    y = obj.GetYaxis().GetBinCenter(iy+1)
                    pedmap.Fill(x,y,obj.GetBinContent(ix+1,iy+1)/float(math.pow(self.rebin,2)))

        tf.Close()
        pedfile.cd()
        pedmap.Write()
        pedmean = ROOT.TH1D('pedmean','pedestal mean',500,97,103)
        pedrms = ROOT.TH1D('pedrms','pedestal RMS',500,0,5)
        for ix in xrange(nx):
            for iy in xrange(ny):
               pedmean.Fill(pedmap.GetBinContent(ix,iy)) 
               pedrms.Fill(pedmap.GetBinError(ix,iy)) 
        pedmean.Write()
        pedrms.Write()
        pedfile.Close()
        print "Pedestal calculated and saved into ",self.pedfile_name

    def reconstruct(self):
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        tf = ROOT.TFile.Open(self.rfile)
        c1 = ROOT.TCanvas('c1','',600,600)
        cc = cameraChannel()
        # loop over events (pictures)
        for ie,e in enumerate(tf.GetListOfKeys()) :
            if ie==options.maxEntries: break
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print "Processing histogram: ",name
            # keep the full resolution 2D image
            #h2zs_fullres = cc.zs(obj,self.pedmap_fr,1,plot=False)
            h2_fullres = obj.Clone(obj.GetName()+'_fr')
            
            # rebin to make the cluster finding fast
            obj.RebinX(self.rebin); obj.RebinY(self.rebin)
            obj.Scale(1./float(math.pow(self.rebin,2)))
            h2zs = cc.zs(obj,self.pedmap)
            print "Analyzing its contours..."
            snfac = SnakesFactory(h2zs,name,options)
            snakes = snfac.getClusters()
            snfac.plotClusterFullResolution(snakes,h2_fullres,self.pedmap_fr)
            snfac.plotProfiles(snakes,h2_fullres,self.pedmap_fr)
            
            #snakes = snfac.getContours(iterations=100)
            #snfac.plotContours(snakes,fill=True)
            #snfac.filledSnakes(snakes)
            
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog h5file1,...,h5fileN [opts] ')
    parser.add_option('-r', '--rebin', dest='rebin', default=4, type='int', help='Rebin factor (same in x and y)')
    parser.add_option(      '--numPedEvents', dest='numPedEvents', default=-1, type='float', help='Use the last n events to calculate the pedestal. Default is all events')

    parser.add_option(      '--max-entries', dest='maxEntries', default=-1, type='float', help='Process only the first n entries')
    parser.add_option(      '--pdir', dest='plotDir', default='./', type='string', help='Directory where to put the plots')

    (options, args) = parser.parse_args()

    inputf = args[0]
    ana = analysis(inputf,options)
    print "Will save plots to ",options.plotDir
    ana.reconstruct()
