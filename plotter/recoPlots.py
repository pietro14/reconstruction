#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch(True)
import os, sys
from plotFile import *

class PlotMaker:
    def __init__(self,tree,tdir,options):
        self._options = options
        self._dir = tdir
        self._tree = tree
        self._cut = options.cut
        ROOT.gROOT.ProcessLine(".x ~/cpp/tdrstyle.cc")
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)

    def run(self,cuts,plots,makeCanvas=True):
        dir = self._dir
        dir.cd()
        pspecs = plots.plots()
        for pspec in pspecs:
            print "    plot: ",pspec.name
            self.printOnePlot(pspec,cuts,makeCanvas=makeCanvas,
                              outputDir=dir, printDir=self._options.printDir)

    def getPlot(self,plotspec,cut,fsplit=None,closeTreeAfter=False):
        ret = makeHistFromBinsAndSpec(plotspec.name,plotspec.expr,plotspec.bins,plotspec)
        drawOpt = "goff"
        if "TProfile" in ret.ClassName(): drawOpt += " PROF";
        self._tree.Draw("%s>>%s" % (plotspec.expr,plotspec.name), cut, drawOpt)
        ret.SetDirectory(0)
        # fold overflow
        if ret.ClassName() in [ "TH1F", "TH1D" ] :
            n = ret.GetNbinsX()
            if plotspec.getOption('IncludeOverflows',True) and ("TProfile" not in ret.ClassName()):
                ret.SetBinContent(1,ret.GetBinContent(0)+ret.GetBinContent(1))
                ret.SetBinContent(n,ret.GetBinContent(n+1)+ret.GetBinContent(n))
                ret.SetBinError(1,hypot(ret.GetBinError(0),ret.GetBinError(1)))
                ret.SetBinError(n,hypot(ret.GetBinError(n+1),ret.GetBinError(n)))
                ret.SetBinContent(0,0)
                ret.SetBinContent(n+1,0)
                ret.SetBinContent(0,0)
                ret.SetBinContent(n+1,0)
            if plotspec.getOption('IncludeOverflow',False) and ("TProfile" not in ret.ClassName()):
                ret.SetBinContent(n,ret.GetBinContent(n+1)+ret.GetBinContent(n))
                ret.SetBinError(n,hypot(ret.GetBinError(n+1),ret.GetBinError(n)))
                ret.SetBinContent(n+1,0)
                ret.SetBinContent(n+1,0)
            if plotspec.getOption('IncludeUnderflow',False) and ("TProfile" not in ret.ClassName()):
                ret.SetBinContent(1,ret.GetBinContent(0)+ret.GetBinContent(1))
                ret.SetBinError(1,hypot(ret.GetBinError(0),ret.GetBinError(1)))
                ret.SetBinContent(0,0)
                ret.SetBinContent(0,0)
            rebin = plotspec.getOption('rebinFactor',0)
            if plotspec.bins[0] != "[" and rebin > 1 and n > 5:
                while n % rebin != 0: rebin -= 1
                if rebin != 1: ret.Rebin(rebin)
            if plotspec.getOption('Density',False):
                for b in xrange(1,n+1):
                    ret.SetBinContent( b, ret.GetBinContent(b) / ret.GetXaxis().GetBinWidth(b) )
                    ret.SetBinError(   b, ret.GetBinError(b) / ret.GetXaxis().GetBinWidth(b) )
                    
        stylePlot(ret,plotspec,plotspec.getOption)
        return ret.Clone(plotspec.name)
        
    def printOnePlot(self,pspec,cuts,makeCanvas=True,outputName=None,outputDir=None,printDir=None):
        options = self._options
        if printDir == None: printDir=self._options.printDir
        if outputDir == None: outputDir = self._dir
        if outputName == None: outputName = pspec.name
        if ROOT.gROOT.FindObject("dummy") != None: ROOT.gROOT.FindObject("dummy").Delete()

        savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kWarning;
                
        # define aspect ratio
        doWide = True if self._options.wideplot or pspec.getOption("Wide",False) else False
        plotformat = (1200,600) if doWide else (600,600)
        sf = 20./plotformat[0]
        ROOT.gStyle.SetPadLeftMargin(600.*0.18/plotformat[0])

        histo = self.getPlot(pspec,cuts)
        ytitle = "Events"
        histo.GetXaxis().SetTitleFont(42)
        histo.GetXaxis().SetTitleSize(0.05)
        histo.GetXaxis().SetTitleOffset(0.9)
        histo.GetXaxis().SetLabelFont(42)
        histo.GetXaxis().SetLabelSize(0.04)
        histo.GetXaxis().SetLabelOffset(0.007)
        histo.GetYaxis().SetTitleFont(42)
        histo.GetYaxis().SetTitleSize(0.05)
        histo.GetYaxis().SetTitleOffset(0.90 if doWide else 1.7)
        histo.GetYaxis().SetLabelFont(42)
        histo.GetYaxis().SetLabelSize(0.04)
        histo.GetYaxis().SetLabelOffset(0.007)
        histo.GetYaxis().SetTitle(pspec.getOption('YTitle',ytitle))
        histo.GetXaxis().SetTitle(pspec.getOption('XTitle',outputName))
        histo.GetXaxis().SetNdivisions(pspec.getOption('XNDiv',510))

        # create canvas

        height = plotformat[1]
        c1 = ROOT.TCanvas(outputName+"_canvas", outputName, plotformat[0], height)
        c1.SetTopMargin(c1.GetTopMargin()*1.2);
        topsize = 0.06*600./height
        c1.Draw()
        c1.SetLogy(True if pspec.hasOption('Logy') else False)
        c1.SetLogz(True if pspec.hasOption('Logz') else False)
        c1.SetWindowSize(plotformat[0] + (plotformat[0] - c1.GetWw()), plotformat[1] + (plotformat[1] - c1.GetWh()));
        if "TH2" in histo.ClassName() or "TProfile2D" in histo.ClassName():
            c1.SetRightMargin(0.20)
            histo.SetContour(100)
            ROOT.gStyle.SetPaintTextFormat(pspec.getOption("PaintTextFormat","g"))
            histo.SetMarkerSize(pspec.getOption("MarkerSize",1))
            if pspec.hasOption('ZMin') and pspec.hasOption('ZMax'):
                histo.GetZaxis().SetRangeUser(pspec.getOption('ZMin',1.0), pspec.getOption('ZMax',1.0))
            # histo.SetMarkerStyle(mca.getProcessOption(p,'MarkerStyle',1))
            # histo.SetMarkerColor(mca.getProcessOption(p,'FillColor',ROOT.kBlack))
            histo.Draw(pspec.getOption("PlotMode","COLZ"))
            # Z axis setting ######################
            histo.GetZaxis().SetTitle(pspec.getOption('ZTitle',ytitle)) # use same content of default ytitle defined above, Events or Events/XX
            histo.GetZaxis().SetTitleFont(42)
            histo.GetZaxis().SetTitleSize(0.055)
            histo.GetZaxis().SetTitleOffset(0.90 if doWide else 1.2)
            histo.GetZaxis().SetLabelFont(42)
            histo.GetZaxis().SetLabelSize(0.05)
            histo.GetZaxis().SetLabelOffset(0.007)
        else:
            histo.Draw("HIST")
        for ext in ['png','pdf']:
            c1.SaveAs("%s/%s.%s" % (printDir, outputName, ext))
        ROOT.gErrorIgnoreLevel = savErrorLevel
        c1.Close()
        
def addPlotterOptions(parser):
    parser.add_option("-n",     "--name",      dest="name",                    default='plot',  help="name");
    parser.add_option("-c",     "--cut",       dest="cut",      type="string", default="1",   help="cut")
    parser.add_option("--pdir", "--print-dir", dest="printDir", type="string", default="plots", help="print out plots in this directory");
    parser.add_option("--select-plot", "--sP", dest="plotselect", action="append", default=[], help="Select only these plots out of the full file")
    parser.add_option("--exclude-plot", "--xP", dest="plotexclude", action="append", default=[], help="Exclude these plots from the full file")
    parser.add_option("--rebin", dest="globalRebin", type="int", default="0", help="Rebin all plots by this factor")
    parser.add_option("--wide", dest="wideplot", action="store_true", default=False, help="Draw a wide canvas")

    
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tree.root cutfile.txt plotfile.txt")
    addPlotterOptions(parser)
    (options, args) = parser.parse_args()

    if len(args)<3:
        print "You need to give at least tree.root cutfile.txt plotfile.txt. Exiting."
        sys.exit(0)

    inputf   = args[0]
    cutfile  = args[1] 
    plotfile = args[2]
    
    tf = ROOT.TFile(inputf)
    tree = tf.Get('Events')

    fileo = open(cutfile,'r')
    cuts = eval(fileo.read())
    cut = ' && '.join([c for k,c in cuts.iteritems()])
    
    if options.cut:
        cut += ' && '+options.cut
    print "Full cut applied = ",cut
        
    if not os.path.exists(options.printDir):
        os.system("mkdir -p "+options.printDir)
        os.system('cp ../utils/index.php '+options.printDir)
        
    plots = PlotFile(plotfile,options)

    outname = 'fngplots.root'
    outfile  = ROOT.TFile(outname,"recreate")
    plotter = PlotMaker(tree,outfile,options)
    plotter.run(cut,plots)
    outfile.Close()
