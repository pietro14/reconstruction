import ROOT
ROOT.gROOT.SetBatch(True)
import os,sys,csv,re

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L functions.cc+");

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

    runmin = 8882#8882
    runmax = 9084#9250

    timeAx = True
    makePlot = True
    if makePlot:
     
        inputfile = "trees-lngs/merged_feruns_8882_9250.root"
        base_selection = "sc_integral<1e4 && sc_integral > 700 && 0.152*sc_length < 50 && sc_rms>6 &&  sc_tgausssigma/sc_lgausssigma>0.6 &&  sc_tgausssigma/sc_lgausssigma<1.1 && R(sc_xmean,sc_ymean)<800"
     
        tfo = ROOT.TFile.Open("ly-history-graph.root","recreate")
     
        ret = ROOT.TGraphErrors()
        ret.SetTitle("")
        ret.SetName("history")
     
        
        tf = ROOT.TFile.Open(inputfile)
        tree = tf.Get("Events")
        histo = ROOT.TH1F("histo","",500,1000,10000)
        
        ir = 0
        with open('../pedestals/runlog_LNGS.csv',"r") as csvfile:
            csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
            next(csvreader)
            for row in list(csvreader):
                time=None
                for ifield in range(len(row)):
                    if re.match(r"\d+\-\d+\-\d+.*",row[ifield].lstrip()):
                        time=row[ifield].lstrip()
                        break
                run = int(row[0])
                if "S003:DATA:Fe" not in row[1] or run<runmin or run>runmax: continue
     
                date = time.replace('"','').split(" ")[0]
                tim = time.replace('"','').split(" ")[1]
                year,month,day=date.split("-")
                hour,minute,sec=tim.split(":")
                print ("run = ",run,"time = ",time)
                tdate = ROOT.TDatime(int(year),int(month),int(day),int(hour),int(minute),int(sec))
                tdate.Print()
                
                histo.Reset()
                sel = "({base})*(run=={run})".format(base=base_selection,run=run)
                tree.Draw("sc_integral>>histo",sel)
                mean = histo.GetMean()
                rms = histo.GetRMS()
                print ("run = ",run,"    mean = ",mean," rms = ",rms)
     
                f = ROOT.TF1('f','gaus',mean-rms,mean+rms) # there is a tail, so +/- 1sigma is sufficient to fit the core
                f.SetParameter(1,mean);
                f.SetParLimits(1,mean-2*rms,mean+2*rms);
                f.SetParameter(2,rms/2.); # there is a tail
                fitr_xmin = mean-rms
                fitr_xmax = mean+rms
                fitRe = histo.Fit(f,'S','',fitr_xmin,fitr_xmax)
                rMean  = f.GetParameter(1)
                rMeanError  = f.GetParError(1)
                rSigma = f.GetParameter(2)
                rSigmaError = f.GetParError(2)

                x = tdate.Convert() if timeAx else run
                ret.SetPoint(ir, x, rMean)
                ret.SetPointError(ir, 0, rMeanError)
                print ("r = ",ir)
                ir += 1
     
        tf.Close()
     
                
        ret.SetMarkerStyle(ROOT.kOpenCircle)
        ret.SetMarkerColor(ROOT.kRed+1)
        ret.SetMarkerSize(0.5)
        ret.SetLineColor(ROOT.kRed+1)
        ret.GetXaxis().SetTitleFont(42)
        ret.GetXaxis().SetLabelOffset(0.03)
        ret.GetXaxis().SetLabelFont(42)
        ret.GetXaxis().SetLabelSize(0.03)
        ret.GetXaxis().SetRangeUser(runmin,runmax)
        ret.GetYaxis().SetTitleFont(42)
        ret.GetYaxis().SetLabelFont(42)
        if timeAx:
            ret.GetXaxis().SetTimeDisplay(1);
            ret.GetXaxis().SetTimeFormat("#splitline{%d\/%m}{%H:%M}");
            ret.GetXaxis().SetTitle("date")
        else:
            ret.GetXaxis().SetTitle("run")
        ret.GetYaxis().SetTitle("LY (counts)") 
     
        c = ROOT.TCanvas('c','',600,600)
        c.SetBottomMargin(0.2); c.SetLeftMargin(0.2); c.SetRightMargin(0.1); 
        ret.Draw("AP")
        
        for ext in ['pdf','png','root']:
            c.SaveAs("light-yield-history-range-%d_%d.%s" % (runmin,runmax,ext))
     
        tfo.cd()
        ret.Write()
        tfo.Close()
     
        print ("Bye")

    else:
        tfgr = ROOT.TFile.Open("ly-history-graph.root")
        gr = tfgr.Get("history")

        c = ROOT.TCanvas('c','',600,600)
        c.SetBottomMargin(0.2); c.SetLeftMargin(0.2); c.SetRightMargin(0.1); 
        gr.Draw("AP")

        gr.Fit('pol2')
        for ext in ['pdf','png','root']:
            c.SaveAs("light-yield-history-range-%d_%d-fit.%s" % (runmin,runmax,ext))
        
        
        
