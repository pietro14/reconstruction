import ROOT
ROOT.gROOT.SetBatch(True)
import os,sys,csv,re
import pandas as pd
import numpy as np

def doLegend(histos,labels,styles,corner="TR",textSize=0.035,legWidth=0.18,legBorder=False,nColumns=1):
    nentries = len(histos)
    (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .89)
    if corner == "TR":
        (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .89)
    elif corner == "TC":
        (x1,y1,x2,y2) = (.5, .7 - textSize*max(nentries-3,0), .5+legWidth, .89)
    elif corner == "TL":
        (x1,y1,x2,y2) = (.2, .7 - textSize*max(nentries-3,0), .2+legWidth, .89)
    elif corner == "BR":
        (x1,y1,x2,y2) = (.85-legWidth, .15 + textSize*max(nentries-3,0), .90, .25)
    elif corner == "BC":
        (x1,y1,x2,y2) = (.5, .15 + textSize*max(nentries-3,0), .5+legWidth, .25)
    elif corner == "BL":
        (x1,y1,x2,y2) = (.2, .23 + textSize*max(nentries-3,0), .33+legWidth, .35)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetNColumns(nColumns)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)  # should make the legend semitransparent (second number is 0 for fully transparent, 1 for full opaque)
    #leg.SetFillStyle(0) # transparent legend, so it will not cover plots (markers of legend entries will cover it unless one changes the histogram FillStyle, but this has other effects on color, so better not touching the FillStyle)
    leg.SetShadowColor(0)
    if not legBorder:
        leg.SetLineColor(0)
        leg.SetBorderSize(0)  # remove border  (otherwise it is drawn with a white line, visible if it overlaps with plots
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    for (plot,label,style) in zip(histos,labels,styles): leg.AddEntry(plot,label,style)
    leg.Draw()
    ## assign it to a global variable so it's not deleted
    global legend_
    legend_ = leg
    return leg

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L functions.cc+");

if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog workdir runs [options] ')
    parser.add_option('-f', '--root-input', dest='inputTree',  default="trees-lngs/merged_feruns_8882_9857.root", help="ROOT file with the Events TTree");
    parser.add_option('-p', '--plot-only',  dest='plotOnly',  action='store_true', default=False, help="don'trun the analysis from the trees, use the results file");
    parser.add_option('-i', '--input-table', dest='inputTable',  default=None, help="Pickle file with the saved panda DataFrame with analysis results");
    parser.add_option('-r', '--run-range',  dest='runRange',  default=[8882,9857], nargs=2, help="minimum and maximum");
    parser.add_option('-a', '--analysis',  dest='analysis',  default='history', help="Type of analysis (default LY history)");
    parser.add_option('-v', '--variable',  dest='variable',  default='fitm', help="variable to be plotted (in case of zscans analysis)");
    (options, args) = parser.parse_args()

    
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

    runmin = int(options.runRange[0])
    runmax = int(options.runRange[1])
    pklfile = "light-yield-history-range-%d_%d.pkl" % (runmin,runmax)

    timeAx = True
    doFit = True
    HV1 = 420

    results = pd.DataFrame(columns=['run','date','vgem1','vgem2','vgem3','nev','nclu','mean','meanerr','rms','fitm','fitmerr','fits','fitserr',    'meanr','meanrerr','rmsr','fitmr','fitmrerr','fitsr','fitsrerr'])
    
    if not options.plotOnly:
     
        base_selection = "sc_integral<2e4 && sc_integral > 1500 && 0.152*sc_length < 50 && sc_rms>8 && sc_tgausssigma/sc_lgausssigma>0.6 &&  sc_tgausssigma/sc_lgausssigma<1.1 && sc_tgausssigma*0.152>0.3 && R(sc_xmean,sc_ymean)<900"
     
        tfo = ROOT.TFile.Open("ly-history-graph.root","recreate")
     
        ret = ROOT.TGraphErrors()
        ret.SetTitle("")
        ret.SetName("history")
     
        
        tf = ROOT.TFile.Open(options.inputTree)
        tree = tf.Get("Events")
        tree.AddFriend("Friends",options.inputTree.replace(".root","_Friend.root"))
        histo = ROOT.TH1F("histo","",500,1000,20000)
        histor = ROOT.TH1F("histor","",500,1000,20000)
        dummy = ROOT.TH1F("dummy","",1,0,1)
        
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
                vgem1 = int(row[-14])
                vgem2 = int(row[-15])
                vgem3 = int(row[-16])
                print ("HV1,2,3 = {v1}, {v2}, {v3}".format(v1=vgem1,v2=vgem2,v3=vgem3))
                #if vgem1 != HV1 or vgem2 != HV1 or vgem3 != HV1: continue
                
                date = time.replace('"','').split(" ")[0]
                tim = time.replace('"','').split(" ")[1]
                year,month,day=date.split("-")
                hour,minute,sec=tim.split(":")
                print ("run = ",run,"time = ",time)
                tdate = ROOT.TDatime(int(year),int(month),int(day),int(hour),int(minute),int(sec))
                tdate.Print()
                
                histo.Reset()
                histor.Reset()
                dummy.Reset()
                sel = "({base})*(run=={run})".format(base=base_selection,run=run)
                runsel = "(run=={run})".format(run=run)
                tree.Draw("sc_integral>>histo",sel)
                tree.Draw("sc_qregr_integral>>histor",sel)
                tree.Draw("0.5>>dummy",runsel)
                mean = histo.GetMean(); meanr = histor.GetMean()
                meanErr = histo.GetMeanError(); meanrErr = histor.GetMeanError()
                rms = histo.GetRMS(); rmsr = histor.GetRMS()
                print ("run = %-20s, mean(raw) = %d(%d) rms(raw) = %d(%d)   nevents = %d"%(str(run),meanr,mean,rmsr,rms,dummy.Integral()))

                rMean = rMeanError = rSigma = rSigmaError = -999
                if doFit:
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

                    # fit the regressed integral distribution
                    f.SetParameter(1,meanr);
                    f.SetParLimits(1,meanr-2*rmsr,meanr+2*rmsr);
                    f.SetParameter(2,rmsr/2.); # there is a tail
                    fitr_xmin = meanr-rmsr
                    fitr_xmax = meanr+rmsr
                    fitRe = histor.Fit(f,'S','',fitr_xmin,fitr_xmax)
                    rMeanr  = f.GetParameter(1)
                    rMeanrError  = f.GetParError(1)
                    rSigmar = f.GetParameter(2)
                    rSigmarError = f.GetParError(2)
                    
                    x = tdate.Convert() if timeAx else run
                    ret.SetPoint(ir, x, rMean)
                    ret.SetPointError(ir, 0, rMeanError)

                else:
                    x = tdate.Convert() if timeAx else run
                    ret.SetPoint(ir, x, mean)
                    ret.SetPointError(ir, 0, meanErr)

                entry = pd.DataFrame.from_dict({'run':[run],'date':[time],'vgem1':[vgem1],'vgem2':[vgem2],'vgem3':[vgem3],
                                                'nev':[dummy.Integral()],'nclu':[histo.Integral()],'mean':[mean],'meanerr':[meanErr],'rms':[rms],
                                                'fitm':[rMean],'fitmerr':[rMeanError],'fits':[rSigma],'fitserr':[rSigmaError],
                                                'meanr':[meanr],'meanrerr':[meanrErr],'rmsr':[rmsr],'fitmr':[rMeanr],'fitmrerr':[rMeanrError],'fitsr':[rSigmar],'fitsrerr':[rSigmarError]})
                results = pd.concat([results,entry], ignore_index=True)
                
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
        results.to_pickle("light-yield-history-range-%d_%d.pkl" % (runmin,runmax))
        
        print ("Bye")

        
        
        
    else:
        if options.inputTable: pklfile=options.inputTable
        print ("Retrieving the data from panda df in: ",pklfile)
        df = pd.read_pickle(pklfile)
        if options.analysis == 'history':
            ret = ROOT.TGraphErrors()
            ret.SetTitle("")
            ret.SetName("history")

            HV1=420
            data = df[(df['vgem1']==HV1)&(df['vgem2']==HV1)&(df['vgem3']==HV1)&(df['run']>=runmin)&(df['run']<=runmax)]
            #data = df[(df['vgem2']==HV1)&(df['vgem3']==HV1)]
            t = data['date'].values
            r = data['run'].values
            if options.variable=='mean':
                y = data['mean'].values
                ye = data['meanerr'].values
            else:
                y = data['fitm'].values
                ye = data['fitmerr'].values                
            for i in range(len(data)):
                date = t[i].replace('"','').split(" ")[0]
                tim = t[i].replace('"','').split(" ")[1]
                year,month,day=date.split("-")
                hour,minute,sec=tim.split(":")
                print ("run = ",r[i],"time = ",t[i])
                tdate = ROOT.TDatime(int(year),int(month),int(day),int(hour),int(minute),int(sec))
                tdate.Print()

                x = tdate.Convert() if timeAx else r[i]
                ret.SetPoint(i, x, y[i])
                ret.SetPointError(i, 0, ye[i])
                
            ret.SetMarkerStyle(ROOT.kOpenCircle)
            ret.SetMarkerColor(ROOT.kRed+1)
            ret.SetMarkerSize(0.5)
            ret.SetLineColor(ROOT.kRed+1)
            ret.GetXaxis().SetTitleFont(42)
            ret.GetXaxis().SetLabelOffset(0.03)
            ret.GetXaxis().SetLabelFont(42)
            ret.GetXaxis().SetLabelSize(0.03)
            ret.GetYaxis().SetTitleFont(42)
            ret.GetYaxis().SetLabelFont(42)
            if timeAx:
                ret.GetXaxis().SetTimeDisplay(1);
                ret.GetXaxis().SetTimeFormat("#splitline{%d\/%m}{%H:%M}");
                ret.GetXaxis().SetTimeOffset(0,"gmt")
                ret.GetXaxis().SetTitle("date")
            else:
                ret.GetXaxis().SetTitle("run")
            ret.GetYaxis().SetTitle("LY (counts)") 
         
            c = ROOT.TCanvas('c','',600,600)
            c.SetBottomMargin(0.2); c.SetLeftMargin(0.2); c.SetRightMargin(0.1); 
            ret.Draw("AP")
            
            for ext in ['pdf','png','root']:
                c.SaveAs("light-yield-history-range-%d_%d.%s" % (runmin,runmax,ext))
         
        if options.analysis == 'zscans':
            zpos = {'3/4': (9364,9372), '10/11': (9378,9386), '17/18': (9390,9440), '24/25': (9745,9753), '32/33': (9735,9742)}
            zmap = {'3/4': 5, '10/11': 15, '17/18': 25, '24/25': 36, '32/33': 48}

            graphs = {}
            retpd = pd.DataFrame(columns=['vgem','z','ly','lyerr'])
            
            for z,rrange in zpos.items():
                rmin,rmax=rrange
                data = df[(df['run']>=rmin)&(df['run']<=rmax)]
                print ("zpos = ",z)
                print (data)

                x = data['vgem1'].unique()
                if (options.variable in ['fitm','mean','fitmr','meanr']):
                    y = [np.mean(data[data['vgem1']==v][options.variable].values) for v in x]
                    print ("y = ",y)
                    ye = [np.mean(data[data['vgem1']==v][options.variable+'err'].values) for v in x]
                elif (options.variable in ['fits','rms','fitsr','rmsr']):
                    y =  [np.mean(data[data['vgem1']==v][options.variable].values)/np.mean(data[data['vgem1']==v]['fitm'].values) for v in x]
                    ye = [np.mean(data[data['vgem1']==v]['fitserr' if options.variable in ['fits','rms'] else 'fitsrerr'].values)/np.mean(data[data['vgem1']==v]['fitm'].values) for v in x]                    
                elif options.variable=='nclu':
                    y = [np.sum(data[data['vgem1']==v]['nclu'].values)/np.sum(data[data['vgem1']==v]['nev'].values) for v in x]
                    print (y)
                    ye = [np.sqrt(np.mean(data[data['vgem1']==v]['nclu'].values))/np.mean(data[data['vgem1']==v]['nev'].values) for v in x]
                else:
                    ye = [0 for v in x]

                
                ret = ROOT.TGraphErrors()
                ret.SetTitle("")
                ret.SetName("hvscan zpos = %s" % z)
                for i in range(len(x)):
                    ret.SetPoint(i, x[i], y[i])
                    ret.SetPointError(i, 0, ye[i])
                    entry = pd.DataFrame.from_dict({'vgem':[x[i]],'z':[zmap[z]],'ly':[y[i]],'lyerr':[ye[i]]})
                    retpd = pd.concat([retpd,entry], ignore_index=True)
                    
                ret.SetMarkerStyle(ROOT.kOpenCircle)
                ret.SetMarkerColor(ROOT.kRed+1)
                ret.SetMarkerSize(0.8)
                ret.SetLineColor(ROOT.kRed+1)
                ret.GetXaxis().SetTitleFont(42)
                ret.GetXaxis().SetLabelOffset(0.03)
                ret.GetXaxis().SetLabelFont(42)
                ret.GetXaxis().SetLabelSize(0.04)
                ret.GetXaxis().SetRangeUser(runmin,runmax)
                ret.GetYaxis().SetTitleFont(42)
                ret.GetYaxis().SetLabelFont(42)
                ret.GetYaxis().SetLabelSize(0.04)
                ret.GetXaxis().SetTitle("V_{GEM1} (V)")
                if (options.variable in ['fitm','mean']):
                    ret.GetYaxis().SetTitle("LY (counts)")
                elif (options.variable=='nclu'):
                    ret.GetYaxis().SetTitle("cluster multiplicity")
                elif (options.variable in ['fits','rms']):
                    ret.GetYaxis().SetTitle("energy resolution")
                else:
                    ret.GetYaxis().SetTitle("arbitrary units")
                graphs[z] = ret
                
            c = ROOT.TCanvas('c','',600,600)
            c.SetBottomMargin(0.2); c.SetLeftMargin(0.2); c.SetRightMargin(0.1);
            ig=0
            colors = {'3/4':ROOT.kRed+1,'10/11':ROOT.kBlack,'17/18':ROOT.kSpring+2, '24/25':ROOT.kAzure+2, '32/33':ROOT.kViolet+2}
            for zp,gr in graphs.items():
                gr.SetMarkerColor(colors[zp])
                gr.SetLineColor(colors[zp])
                if (options.variable in ['fitm','mean','fitmr','meanr']):
                    gr.GetYaxis().SetRangeUser(0,15e3)
                elif (options.variable in ['nclu']):
                    gr.GetYaxis().SetRangeUser(0,3)
                elif (options.variable in ['fits','rms','fitsr','rmsr']):
                    gr.GetYaxis().SetRangeUser(0,0.5)
                else:
                    pass
                if ig==0: gr.Draw("AP")
                else: gr.Draw("P")
                ig+=1

            responses = [gr for gr in graphs.values()]
            titles = ["z = %d cm" % zmap[k] for k in graphs.keys()]
            styles = ['pl' for k in graphs.keys()]
            legend = doLegend(responses,titles,styles,corner="TL")
            
            for ext in ['pdf','png','root']:
                c.SaveAs("%s-zhvscans.%s" % (options.variable,ext))
            retpd.to_pickle("%s-zhvscans.pkl" % options.variable)
