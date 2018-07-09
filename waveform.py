#!/usr/bin/env python

import os,math,sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
from scipy.signal import find_peaks

class PeakFinder:
    def __init__(self,tgraph,xmin=-1,xmax=-1):
        self.importTGraph(tgraph,xmin,xmax)
        self.name = tgraph.GetName()
        
    def importTGraph(self,tgraph,xmin,xmax):
        # transform to positive signals
        y = np.array([-y for y in tgraph.GetY()])
        x = np.array(tgraph.GetX())
        xmax = xmax if xmax>0 else x[-1]
        ix = np.array([i for i,v in enumerate(x) if v>xmin and v<xmax])
        self.x = np.array(x[ix])
        self.y = np.array(y[ix])
        self.binsize = self.x[1]-self.x[0]
        
    def findPeaks(self,thr,mindist,prominence=1,width=5):
        peaks, properties = find_peaks(self.y, distance=mindist, height=thr, prominence=prominence,width=width)
        self.peaks = peaks
        self.properties = properties
        return peaks

    def plotpy(self,pdir='./'):
        import matplotlib.pyplot as plt

        # plot data and the found peaks
        plt.plot(self.x,self.y)
        plt.plot(self.getPeakTimes(), self.y[self.peaks], "x")
        plt.plot(self.x, np.zeros_like(self.y), "--", color="gray")

        # plot some properties
        plt.vlines(x=self.getPeakTimes(), ymin=self.y[self.peaks] - self.getProminences(),
                   ymax = self.y[self.peaks], color = "C1")
        plt.hlines(y=self.getHMs(), xmin=self.getPeakBoundaries('left'),
                   xmax=self.getPeakBoundaries('right'), color = "C1")
        
        plt.xlabel('Time (ns)')
        plt.ylabel('amplitude (mV)')
        for ext in ['png','pdf']:
            plt.savefig('{pdir}/{name}.{ext}'.format(pdir=pdir,name=self.name,ext=ext))
        plt.gcf().clear()

    def getPeakBoundaries(self,side):
        if side=='left': return np.array([self.x[int(x)] for x in self.properties["left_ips"]])
        return np.array([self.x[int(x)] for x in self.properties["right_ips"]])
        
    def getFWHMs(self):
        return self.properties["widths"]
        
    def getHMs(self):
        return self.properties["width_heights"]

    def getPeakTimes(self):
        return self.x[self.peaks]

    def getProminences(self):
        return self.properties["prominences"]

    def getAmplitudes(self):
        return self.properties["peak_heights"]


from cameraChannel import cameraGeometry
class PMTSignal:
    def __init__(self,tgraph,clusters,options):
        self.waveform = tgraph
        self.clusters = clusters
        self.options = options
        
    def plotNice(self):
        sig_width = 150 #ns
        sig_min = 6150 # at least at FNG with DAQ

        canv = ROOT.TCanvas("cfr","",600,600)
        canv.SetLeftMargin(0.20)
        canv.SetBottomMargin(0.15)
        self.waveform.Draw('AL')
        self.waveform.GetXaxis().SetRangeUser(sig_min,sig_min+sig_width)
        self.waveform.GetXaxis().SetTitle('Time (ns)')
        self.waveform.GetYaxis().SetTitle('Amplitude (mV)')

        maxwidth = 0
        if len(self.clusters): maxwidth = max([cl.widths['long'] for cl in self.clusters]) # mm
        title = 'N clusters = {nclu}, max length = {maxl:.1f}mm'.format(nclu=len(self.clusters), maxl=maxwidth)
        self.waveform.SetTitle(title)

        for ext in ['png','pdf']:
            canv.SaveAs('{od}/{name}.{ext}'.format(od=self.options.plotDir,name=self.waveform.GetName(),ext=ext))
        

if __name__ == '__main__':

    inputf = sys.argv[1]
    print "testing ",inputf 

    tf = ROOT.TFile(inputf)
    # sampling was 5 GHz (5/ns). Separate peaks of at least 1ns

    threshold = 10 # min threshold for a signal
    min_distance_peaks = 10 # number of samples (10 samples = 2ns)
    prominence = 2 # noise seems ~1 mV
    width = 5 # minimal width of the signal

    # single peak example
    gr = tf.Get('wfm_run00070_ev310')
    pf = PeakFinder(gr,6160,6300)
    pf.findPeaks(threshold,min_distance_peaks,prominence)
    pf.plotpy()
    
    # two clear peaks example
    gr = tf.Get('wfm_run00070_ev39')
    pf = PeakFinder(gr,6160,6300)
    pf.findPeaks(threshold,min_distance_peaks,prominence)
    pf.plotpy()

    # mess
    gr = tf.Get('wfm_run00070_ev9')
    pf = PeakFinder(gr,6160,6300)
    pf.findPeaks(threshold,min_distance_peaks,prominence)
    pf.plotpy()
