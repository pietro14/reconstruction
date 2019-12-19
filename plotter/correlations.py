#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch(True)

filein = "reco.root"

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kRainBow)

f = ROOT.TFile.Open(filein)
c1 = ROOT.TCanvas("c1","",600,600)

npeaksVtlength = ROOT.TH2D("npeaksVtlength_prevEvent","",15,0,70,15,1,15)
npeaks2Vtlength = ROOT.TH2D("npeaksVtlength_thisEvent","",15,0,70,15,1,15)

nPeakPrev=0
for ie,event in enumerate(f.Events):
    if ie>0:
        if event.nTrack>0:
            npeaksVtlength.Fill(event.track_length[0],nPeakPrev)
            npeaks2Vtlength.Fill(event.track_length[0],event.nPeak)
            #print "ev = ",event.event,"np current = ",event.nPeak," np prev = ",nPeakPrev, "t length = ",event.track_length[0]
    nPeakPrev=event.nPeak
        
th2s = [npeaksVtlength,npeaks2Vtlength]
for h2 in th2s:
    h2.GetXaxis().SetTitle('track length (mm)')
    h2.GetYaxis().SetTitle('number of PMT peaks')
    h2.Draw('colz')
    h2.GetZaxis().SetRangeUser(0,20)
    print h2.GetName(),"  correlation = ",h2.GetCorrelationFactor()
    for ext in ['pdf','png']:
        c1.SaveAs("{name}.{ext}".format(name=h2.GetName(),ext=ext))


