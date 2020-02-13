import ROOT as rt
rt.gROOT.LoadMacro('./libCpp/histFitter.C+')

from ROOT import histFitter

import re
import math


def histFitterNominal( sample, fitWorkspaceParam  ):

    fitWorkspaceFunc = [
        "Gaussian::sigPdf(x,meanS,sigmaS)",
        ]

    fitWorkspace = []
    fitWorkspace.extend(fitWorkspaceParam)
    fitWorkspace.extend(fitWorkspaceFunc)
    
    ## init fitter
    print "sample.path = ",sample.path
    infile = rt.TFile.Open( sample.path, "read")
    hData = infile.Get( sample.hist )
    fitter = histFitter( hData, 'hist_data' )
    infile.Close()

    ## setup
    fitter.useMinos()
    rootfile = rt.TFile.Open(sample.fitOutput, 'update')
    fitter.setOutputFile( rootfile )

    ## PDF for Bkg from template
    fileBkg  = rt.TFile.Open(sample.pathBkg,'read')
    hBkg = fileBkg.Get( sample.histBkg )
    fitter.setBkgPDF(hBkg)
    fileBkg.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in fitWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = 'boh'
    fitter.fits(title)
    rootfile.Close()

    ## give info on efficiency for a set of cuts
    cuts = [16,17,18,19,20]
    for c in cuts:
        print "signal Effi for x>",c," = ",fitter.efficiency(c,'sigPdf')
        print "bkg Effi for x>",c," = ",fitter.efficiency(c,'bkgPdf')
    

def histPlotter( filename, plotDir):
    print 'opening ', filename
    rootfile = rt.TFile.Open(filename,"read")

    print '  get canvas: hist_data_Canv'
    c = rootfile.Get( 'hist_data_Canv' )
    c.SaveAs( '%s/%s.pdf' % (plotDir,'hist_data') )

    rootfile.Close()
    
