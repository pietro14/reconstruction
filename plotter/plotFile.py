#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
from copy import *
import ROOT

class PlotFile:
    def __init__(self,fileName,options):
        self._options = options
        self._plots = []
        defaults = {}
        infile = open(fileName,'r')
        for line in infile:
            if re.match("\s*#.*", line) or len(line.strip())==0: continue
            while line.strip()[-1] == "\\":
                line = line.strip()[:-1] + infile.next()
            extra = {}
            if ";" in line:
                (line,more) = line.split(";")[:2]
                more = more.replace("\\,",";")
                for setting in [f.strip().replace(";",",") for f in more.split(',')]:
                    if "=" in setting:
                        # in following line, if setting has more than one '=' an error will occur, e.g., if you have XTitle="muon isolation (#DeltaR=0.4)"                                                                                                                      
                        # therefore, split only on the first occurrence of '='                                                                                                                                                                                                  
                        (key,val) = [f.strip() for f in setting.split("=",1)]
                        extra[key] = eval(val)
                    else: extra[setting] = True
            line = re.sub("#.*","",line)
            field = [f.strip().replace(";",":") for f in line.replace("::",";;").replace("\\:",";").split(':')]
            if len(field) == 1 and field[0] == "*":
                if len(self._plots): raise RuntimeError, "PlotFile defaults ('*') can be specified only before all plots"
                print "Setting the following defaults for all plots: "
                for k,v in extra.iteritems():
                    print "\t%s: %r" % (k,v)
                    defaults[k] = v
                continue
            else:
                for k,v in defaults.iteritems():
                    if k not in extra: extra[k] = v
            if len(field) <= 2: continue
            if len(options.plotselect):
                skipMe = True
                for p0 in options.plotselect:
                    for p in p0.split(","):
                        if re.match(p+"$", field[0]): skipMe = False
                if skipMe: continue
            if len(options.plotexclude):
                skipMe = False
                for p0 in options.plotexclude:
                    for p in p0.split(","):
                        if re.match(p+"$", field[0]): skipMe = True
                if skipMe: continue
            if options.globalRebin: extra['rebinFactor'] = options.globalRebin
            self._plots.append(PlotSpec(field[0],field[1],field[2],extra))
    def plots(self):
        return self._plots[:]
    
class PlotSpec:
    def __init__(self,name,expr,bins,opts):
        self.name = name
        self.expr = expr
        self.bins = bins
        self.opts = opts
        self.logs = {}
    def hasOption(self,name):
        return (name in self.opts)
    def getOption(self,name,default=None):
        return self.opts[name] if (name in self.opts) else default
    def setOption(self,name,value):
        self.opts[name] = value
    def setLog(self,name,value):
        self.logs[name] = value
    def allLogs(self):
        return self.logs.iteritems()

def stylePlot(plot,spec,getOption):
    ## Sample specific-options, from self
    if getOption('FillColor',None) != None:
        plot.SetFillColor(getOption('FillColor',0))
        plot.SetFillStyle(getOption('FillStyle',1001))
    else:
        plot.SetFillStyle(0)
        plot.SetLineWidth(getOption('LineWidth',1))
    plot.SetLineColor(getOption('LineColor',1))
    plot.SetLineStyle(getOption('LineStyle',1))
    plot.SetMarkerColor(getOption('MarkerColor',1))
    plot.SetMarkerStyle(getOption('MarkerStyle',20))
    plot.SetMarkerSize(getOption('MarkerSize',1.1))
    ## Plot specific-options, from spec
    if "TH3" not in plot.ClassName():
        plot.GetYaxis().SetTitle(spec.getOption('YTitle',"Events"))
        plot.GetXaxis().SetTitle(spec.getOption('XTitle',spec.name))
        plot.GetXaxis().SetNdivisions(spec.getOption('XNDiv',510))
        plot.GetXaxis().SetMoreLogLabels(True)

def makeHistFromBinsAndSpec(name,expr,bins,plotspec):
    profile1D      = plotspec.getOption('Profile1D',False) if plotspec != None else False
    profile2D      = plotspec.getOption('Profile2D',False) if plotspec != None else False
    nvars = expr.replace("::","--").count(":")+1
    if nvars == 1 or (nvars == 2 and profile1D):
        if bins[0] == "[":
            edges = [ float(f) for f in bins[1:-1].split(",") ]
            if profile1D:
                histo = ROOT.TProfile(name,name,len(edges)-1,array('f',edges))
            else:
                histo = ROOT.TH1D(name,name,len(edges)-1,array('f',edges))
        else:
            (nb,xmin,xmax) = bins.split(",")
            if profile1D:
                histo = ROOT.TProfile(name,name,int(nb),float(xmin),float(xmax))
            else:
                histo = ROOT.TH1D(name,name,int(nb),float(xmin),float(xmax))
    elif nvars == 2 or (nvars == 3 and profile2D):
        if bins[0] == "[":
            xbins, ybins = bins.split("*")
            xedges = [ float(f) for f in xbins[1:-1].split(",") ]
            yedges = [ float(f) for f in ybins[1:-1].split(",") ]
            if profile2D:
                histo = ROOT.TProfile2D(name,name,len(xedges)-1,array('d',xedges),len(yedges)-1,array('d',yedges))
            else:
                histo = ROOT.TH2D(name,name,len(xedges)-1,array('f',xedges),len(yedges)-1,array('f',yedges))
        else:
            (nbx,xmin,xmax,nby,ymin,ymax) = bins.split(",")
            if profile2D:
                histo = ROOT.TProfile2D(name,name,int(nbx),float(xmin),float(xmax),int(nby),float(ymin),float(ymax))
            else:
                histo = ROOT.TH2D(name,name,int(nbx),float(xmin),float(xmax),int(nby),float(ymin),float(ymax))
    elif nvars == 3:
        ez,ey,ex = [ e.replace("--","::") for e in expr.replace("::","--").split(":") ]
        if bins[0] == "[":
            xbins, ybins, zbins = bins.split("*")
            xedges = [ float(f) for f in xbins[1:-1].split(",") ]
            yedges = [ float(f) for f in ybins[1:-1].split(",") ]
            zedges = [ float(f) for f in zbins[1:-1].split(",") ]
            histo = ROOT.TH3D(name,name,len(xedges)-1,array('f',xedges),len(yedges)-1,array('f',yedges),len(zedges)-1,array('f',zedges))
        else:
            (nbx,xmin,xmax,nby,ymin,ymax,nbz,zmin,zmax) = bins.split(",")
            histo = ROOT.TH3D(name,name,int(nbx),float(xmin),float(xmax),int(nby),float(ymin),float(ymax),int(nbz),float(zmin),float(zmax))
        histo.GetXaxis().SetTitle(ex)
        histo.GetYaxis().SetTitle(ey)
        histo.GetZaxis().SetTitle(ez)
    else:
        raise RuntimeError, "Can't make a plot with %d dimensions" % nvars
    histo.Sumw2()
    return histo
