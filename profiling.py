import os,math,sys,ctypes
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
from scipy.signal import find_peaks,peak_widths

class simplePeak:
    def __init__(self,ampli,prominence,mean,fwhm):
        self.amplitude = ampli
        self.prominence = prominence
        self.mean = mean
        self.fwhm = fwhm
    def __repr__(self):
        return "(Ampli={ampli:.2f}, Prom={prom:.2f}, Mean={mean:.2f}, FWHM={fwhm:.2f})".format(ampli=self.amplitude,prom=self.prominence,mean=self.mean,fwhm=self.fwhm)

class PeakFinder:
    def __init__(self,graph,xmin=None,xmax=None,rebin=None,negative=True):
        if graph.InheritsFrom('TH1'):
            self.importTH1(graph,xmin,xmax,rebin,negative)
        self.name = graph.GetName()
        self.xmin = xmin; self.xmax=xmax

    def importTH1(self,th1,xmin,xmax,rebin,negative=True):
        if rebin:
            if th1.InheritsFrom('TProfile'):
                print("WARNING! Rebinning for TProfile not implemented yet!")
            else:
                th1.Rebin(rebin)
        ysign = -1 if negative else 1
        x    = np.array([th1.GetXaxis().GetBinCenter(b) for b in range(1,th1.GetNbinsX()+1)])
        y    = np.array([ysign*th1.GetBinContent(b) for b in range(1,th1.GetNbinsX()+1)])
        yerr = np.array([th1.GetBinError(b) for b in range(1,th1.GetNbinsX()+1)])
        self.setData(x,y,xmin,xmax,yerr)
        
    def setData(self,x,y,xmin,xmax,yerr=np.array([])):
        xmax = xmax if xmax!=None else x[-1]
        xmin = xmin if xmin!=None else x[0]
        ix = np.array([i for i,v in enumerate(x) if v>xmin and v<xmax])
        if len(ix):
            self.x = np.array(x[ix])
            self.y = np.array(y[ix])
            if len(yerr)==len(ix):
                self.yerr = np.array(yerr[ix])
            else:
                self.yerr = np.zeros(len(ix))
        else:
            self.x = self.y = self.yerr = np.array([])
        self.binsize = self.x[1]-self.x[0] if len(self.x)>1 else 0
        
    def findPeaks(self,thr,mindist,prominence=1,width=5):
        peaks, properties = find_peaks(self.y, distance=mindist, height=thr, prominence=prominence,width=width)
        self.peaks = peaks
        self.properties = properties
        self.setTot(thr)
        return peaks

    def plotpy(self,pdir='./',xlabel='Time (ns)',ylabel='amplitude (mV)'):
        import matplotlib.pyplot as plt
        ## enable TeX 
        from matplotlib import rc
        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        rc('text', usetex=True)

        # plot data and the found peaks
        plt.errorbar(self.x, self.y, self.yerr, ls='', ecolor='lightgrey',elinewidth=1, marker='o', mfc = 'black', ms=3, mew=0)
        plt.plot(self.getPeakTimes(), self.y[self.peaks], "x")
        plt.plot(self.x, np.zeros_like(self.y), "--", color="gray")

        # plot some properties
        plt.vlines(x=self.getPeakTimes(), ymin=self.y[self.peaks] - self.getProminences(),
                   ymax = self.y[self.peaks], color = "C1")
        plt.hlines(y=self.getHMs(), xmin=self.getPeakBoundaries('left'),
                   xmax=self.getPeakBoundaries('right'), color = "C1")        
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        for ext in ['pdf','png']:
            plt.savefig('{pdir}/{name}.{ext}'.format(pdir=pdir,name=self.name,ext=ext))
        plt.gcf().clear()

    def getPeakBoundaries(self,side):
        if side=='left': return np.array([self.x[int(x)] for x in self.properties["left_ips"]])
        return np.array([self.x[int(x)] for x in self.properties["right_ips"]])

    def getFWHMs(self):
        return self.properties["widths"]

    def getFullWidths(self):
        self.widths_full = peak_widths(self.y, self.peaks, rel_height=1)        
        return self.widths_full

    def getTimes(self,side='rise'):
        if side=='rise': index=2
        elif side=='fall': index=3
        else:
            print("ERROR! Side should be either rise or fall. Exiting.")
            return []
        if not hasattr(self,'widths_full'):
            self.getFullWidths()
        # intersection points are interpolated
        times = np.array([self.x[int(x)] for x in self.widths_full[index]])
        return times

    def getHMs(self):
        return self.properties["width_heights"]

    def getPeakTimes(self):
        return self.x[self.peaks]

    def getProminences(self):
        return self.properties["prominences"]

    def getAmplitudes(self):
        return self.properties["peak_heights"]

    def setTot(self,threshold=0):
        x0,x1=(-1,-1)
        # for robustness, look for the rise from the left and for the fall from the right
        for i,y in enumerate(self.y):
            if y>threshold:
                x0 = self.x[i]
                break
        for i in range(1,len(self.y)):
            ip = -i
            if self.x[ip]<=x0:
                x1 = x0
                break
            y = self.y[ip]
            if y>threshold:
                x1=self.x[ip]
                break
        if self.xmin is not None and self.xmax is not None:
            self.x0 = np.nanmax(np.array([x0,self.xmin]))
            self.x1 = np.nanmin(np.array([x1,self.xmax]))
        else:
            self.x0 = -999
            self.x1 = -999

    def getTot(self):
        return self.x1-self.x0

    def getIntegral(self):
        # range of x with y over threshold
        ix = np.array([i for i,v in enumerate(self.x) if v>self.x0 and v<self.x1])
        if len(ix)==0:
            return 0
        else:
            return sum(self.y[ix])