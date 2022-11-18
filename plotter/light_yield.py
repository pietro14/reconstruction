#!/usr/bin/env python

import ROOT
ROOT.gROOT.SetBatch(True)
import os,sys,csv

def getLYGraph(runs,color):
    ret = ROOT.TGraphErrors()
    ret.SetTitle("")
    r = 0
    print(runs)
    for run,time in runs.items():
        fname = '/Users/emanuele/data/cygnus/RECO/lngs_Nov22/reco_run%05d_3D.root' % run
        if not os.path.isfile(fname):
            print ("run ",run," not present, skip...")
            continue
        date = time.replace('"','').split(" ")[1]
        tim = time.replace('"','').split(" ")[2]
        year,month,day=date.split("-")
        hour,minute,sec=tim.split(":")
        print ("run = ",run,"time = ",time)
        tdate = ROOT.TDatime(int(year),int(month),int(day),int(hour),int(minute),int(sec))
        
        f = ROOT.TFile.Open(fname)
        tree = f.Get('Events')

        if ROOT.gROOT.FindObject("ly") != None: ROOT.gROOT.FindObject("ly").Delete()
        histo = ROOT.TH1F("ly","ly",25,2e3,1e4); histo.Sumw2()
        
        tree.Draw("sc_integral>>ly",cut)

        mean = histo.GetMean()
        meanerr = histo.GetMeanError()
        nev = histo.GetEntries()

        ret.SetPoint(r,tdate.Convert(),mean)
        ret.SetPointError(r,0,meanerr)
        r += 1
        f.Close()

    ret.SetMarkerStyle(ROOT.kOpenCircle)
    ret.SetMarkerColor(color)
    ret.SetMarkerSize(0.5)
    ret.SetLineColor(color)
    ret.GetXaxis().SetTitleFont(42)
    ret.GetXaxis().SetLabelFont(42)
    ret.GetXaxis().SetLabelSize(0.03)
    ret.GetYaxis().SetTitleFont(42)
    ret.GetYaxis().SetLabelFont(42)
    ret.GetXaxis().SetTimeDisplay(1);
    ret.GetXaxis().SetTimeFormat("%H:%M");
    ret.GetXaxis().SetTitle("time") 
    ret.GetYaxis().SetTitle("Light yield (counts)") 
    
    
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

    cut = 'sc_rms>6 && sc_width/sc_length>0.8 && sc_integral<1e4 && TMath::Hypot(sc_xmean-2304./2.,sc_ymean-2304./2.)<900 && sc_integral>2e3'

    runs_night = {}
    with open('runs_day.txt',"r") as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
        next(csvreader)
        for row in list(csvreader):
            runs_night[int(row[0])] = row[1]
    #print(("dic = ",runs_night))
    
    night_result = getLYGraph(runs_night,ROOT.kAzure-6)
    # cosmics_result = getLYGraph(cosmics_runs,ROOT.kRed+1,4)

    c = ROOT.TCanvas('c','light yield',600,600)
    c.SetBottomMargin(0.2); c.SetLeftMargin(0.2); c.SetRightMargin(0.1); 
    night_result.Draw("AP")
    # btf_result.GetYaxis().SetRangeUser(0,400)
    # cosmics_result.Draw("P")

    # leg = ROOT.TLegend(0.55,0.65,0.76,0.82);
    # leg.AddEntry(btf_result,"BTF 2017")
    # leg.AddEntry(cosmics_result,"cosmics 2018")
    # leg.Draw()
    
    for ext in ['pdf','png','root']:
        c.SaveAs("light_yield.%s" % ext)
    
    
