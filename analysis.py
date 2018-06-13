#!/usr/bin/env python

import os,math
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.SetBatch(True)


from snakes import SnakesFactory

class analysis:

    def __init__(self,rfile,options):
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
        
    def zs(self,th2,pedmap,nsigma=5,plot=False):
        nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
        xmin,xmax=(th2.GetXaxis().GetXmin(),th2.GetXaxis().GetXmax())
        ymin,ymax=(th2.GetYaxis().GetXmin(),th2.GetYaxis().GetXmax())
        th2_zs = ROOT.TH2D(th2.GetName()+'_zs',th2.GetName()+'_zs',nx,xmin,xmax,ny,ymin,ymax)
        th2_zs.SetDirectory(None)
        for ix in xrange(1,nx+1):
            for iy in xrange(1,ny+1):
                if not self.isGoodChannel(pedmap,ix,iy): continue
                ped = pedmap.GetBinContent(ix,iy)
                noise = pedmap.GetBinError(ix,iy)
                z = max(th2.GetBinContent(ix,iy)-ped,0)
                if z>nsigma*noise:
                    th2_zs.SetBinContent(ix,iy,z)
                    #print "x,y,z=",ix," ",iy," ",z,"   noise = ",noise
        #th2_zs.GetZaxis().SetRangeUser(0,1)
        if plot:
            canv = ROOT.TCanvas('zs','',600,600)
            th2_zs.Draw('colz')
            for ext in ['png','pdf']:
                canv.SaveAs('{name}.{ext}'.format(name=th2.GetName()+'_zs',ext=ext))
        return th2_zs

    def calcPedestal(self,maxImages=-1,alternativeRebin=-1):
        nx=ny=2048
        rebin = self.rebin if alternativeRebin<0 else alternativeRebin
        nx=int(nx/rebin); ny=int(ny/rebin); 
        pedfile = ROOT.TFile.Open(self.pedfile_name,'recreate')
        pedmap = ROOT.TProfile2D('pedmap','pedmap',nx,0,nx,ny,0,ny,'s')
        tf = ROOT.TFile.Open(self.rfile)
        for i,e in enumerate(tf.GetListOfKeys()):
            if maxImages>-1 and i<len(tf.GetListOfKeys())-maxImages: continue
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print "Calc pedestal with event: ",name
            obj.RebinX(rebin); obj.RebinY(rebin); 
            for ix in xrange(nx):
                for iy in xrange(ny):
                    pedmap.Fill(ix,iy,obj.GetBinContent(ix+1,iy+1)/float(math.pow(self.rebin,2)))

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

    def isGoodChannel(self,pedmap,ix,iy):
        pedval = pedmap.GetBinContent(ix,iy)
        pedrms = pedmap.GetBinError(ix,iy)
        if pedval > 110: return False
        if pedrms < 0.2: return False
        if pedrms > 5: return False
        return True

    def reconstruct(self):
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        tf = ROOT.TFile.Open(self.rfile)
        c1 = ROOT.TCanvas('c1','',600,600)
        # loop over events (pictures)
        for ie,e in enumerate(tf.GetListOfKeys()) :
            if ie==options.maxEntries: break
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print "Processing histogram: ",name
            # keep the full resolution 2D image
            h2zs_fullres = self.zs(obj,self.pedmap_fr,1,plot=False)
            
            # rebin to make the cluster finding fast
            obj.RebinX(self.rebin); obj.RebinY(self.rebin)
            obj.Scale(1./float(math.pow(self.rebin,2)))
            h2zs = self.zs(obj,self.pedmap)
            print "Analyzing its contours..."
            snfac = SnakesFactory(h2zs,name,options)
            snakes = snfac.getClusters()
            snfac.plotClusterFullResolution(snakes,h2zs_fullres)
            snfac.plotProfiles(snakes,h2zs_fullres)
            
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
