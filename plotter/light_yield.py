#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch(True)
import sys

def getLYGraph(runs,color,rebin=1):
    ret = ROOT.TGraphErrors(len(runs))
    ret.SetTitle("")
    r = 0
    for run,dist in runs.iteritems():
        fname = 'reco_run%s.root' % run
        print "Using file ",fname
        f = ROOT.TFile.Open(fname)
        tree = f.Get('Events')
        print "Has ",tree.GetEntries()," entries."

        if ROOT.gROOT.FindObject("ly") != None: ROOT.gROOT.FindObject("ly").Delete()
        histo = ROOT.TH1F("ly","ly",20,0,1000); histo.Sumw2()
        
        tree.Draw("track_integral/track_length>>ly",cut)

        mean = histo.GetMean()*rebin
        rms = histo.GetRMS()*rebin
        nev = histo.GetEntries()

        print "run %s: ly = %d +/- %d" % (run,mean,rms)
        print "from %d tracks" % nev
        print "distance = %d cm" %dist
        
        ret.SetPoint(r,dist,mean)
        ret.SetPointError(r,0,rms)
        r += 1
        f.Close()

    ret.SetMarkerStyle(ROOT.kOpenCircle)
    ret.SetMarkerColor(color)
    ret.SetMarkerSize(3)
    ret.SetLineColor(color)
    ret.GetXaxis().SetTitle("Distance (cm)") 
    ret.GetYaxis().SetTitle("Light yield (Photons/mm)") 
    
    return ret

if __name__ == '__main__':

    # use large Times-Roman fonts
    ROOT.gStyle.SetTitleFont(132,"xyz");  # set the all 3 axes title font
    ROOT.gStyle.SetTitleFont(132," ");    # set the pad title font
    ROOT.gStyle.SetTitleSize(0.06,"xyz"); # set the 3 axes title size
    ROOT.gStyle.SetTitleSize(0.06," ");   # set the pad title size
    ROOT.gStyle.SetLabelFont(132,"xyz");
    ROOT.gStyle.SetLabelSize(0.05,"xyz");
    ROOT.gStyle.SetTextFont(132);
    ROOT.gStyle.SetTextSize(0.08);
    ROOT.gStyle.SetStatFont(132);

    cut = 'nTrack<4'

    btf_runs = {"058" : -7.0,
                "065" :  0.0,
                "075" :  8.0}
    cosmics_runs = {"121" : -7.0,
                    "123" : 8.0}
    
    btf_result = getLYGraph(btf_runs,ROOT.kAzure-6)
    cosmics_result = getLYGraph(cosmics_runs,ROOT.kRed+1,4)

    c = ROOT.TCanvas('c','light yield',600,600)
    c.SetBottomMargin(0.3); c.SetLeftMargin(0.2); c.SetRightMargin(0.2); 
    btf_result.Draw("AP")
    btf_result.GetYaxis().SetRangeUser(0,400)
    cosmics_result.Draw("P")

    leg = ROOT.TLegend(0.55,0.65,0.76,0.82);
    leg.AddEntry(btf_result,"BTF 2017")
    leg.AddEntry(cosmics_result,"cosmics 2018")
    leg.Draw()
    
    for ext in ['pdf','png']:
        c.SaveAs("light_yield.%s" % ext)
    
    
