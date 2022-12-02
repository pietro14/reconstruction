import ROOT
ROOT.gROOT.SetBatch(True)
import os,sys,csv,re


def getGraph(name,title,runs):

    ret = ROOT.TGraphErrors()
    ret.SetTitle("")
    r = 0
    print(runs)
    for run,time in runs.items():
        fname = "../pedestals/pedmap_run%d_rebin1.root" % run
        if not os.path.isfile(fname):
            print ("run ",run," not present, skip...")
            continue
        date = time.replace('"','').split(" ")[0]
        tim = time.replace('"','').split(" ")[1]
        year,month,day=date.split("-")
        hour,minute,sec=tim.split(":")
        print ("run = ",run,"time = ",time)
        tdate = ROOT.TDatime(int(year),int(month),int(day),int(hour),int(minute),int(sec))
        
        f = ROOT.TFile.Open(fname)

        histo = f.Get(name)
        
        mean = histo.GetMean()
        meanerr = histo.GetMeanError()

        ret.SetPoint(r,tdate.Convert(),mean)
        ret.SetPointError(r,0,meanerr)
        r += 1
        f.Close()

    ret.SetMarkerStyle(ROOT.kOpenCircle)
    ret.SetMarkerColor(ROOT.kRed+1)
    ret.SetMarkerSize(0.5)
    ret.SetLineColor(ROOT.kRed+1)
    ret.GetXaxis().SetTitleFont(42)
    ret.GetXaxis().SetLabelFont(42)
    ret.GetXaxis().SetLabelSize(0.03)
    ret.GetYaxis().SetTitleFont(42)
    ret.GetYaxis().SetLabelFont(42)
    ret.GetXaxis().SetTimeDisplay(1);
    ret.GetXaxis().SetTimeFormat("%d\/%m");
    ret.GetXaxis().SetTitle("time") 
    ret.GetYaxis().SetTitle(title) 

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

    runs = {}
    with open('../pedestals/runlog_LNGS_SummerPedestals_1.csv',"r") as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
        next(csvreader)
        for row in list(csvreader):
            date=None
            for ifield in range(len(row)):
                if re.match(r"\d+\-\d+\-\d+.*",row[ifield].lstrip()):
                    date=row[ifield].lstrip()
                    break
            runs[int(row[0])] = date

    variables = {'pedrms':    "pedestal RMS",
                 'pedmean':   "pedestal mean"}
    
    for var,title in variables.items():
        result = getGraph(var,title,runs)

        c = ROOT.TCanvas('c','',600,600)
        c.SetBottomMargin(0.2); c.SetLeftMargin(0.2); c.SetRightMargin(0.1); 
        result.Draw("AP")

        for ext in ['pdf','png','root']:
            c.SaveAs("%s.%s" % (var,ext))
