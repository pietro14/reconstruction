import os, math, optparse, ROOT
from array import array
from simple_plot import getCanvas, doLegend
from mcPlots import doTinyCmsPrelim

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

def formatHisto(histo):
    histo.GetXaxis().SetTitleFont(42)
    histo.GetXaxis().SetTitleSize(0.05)
    histo.GetXaxis().SetTitleOffset(1.1)
    histo.GetXaxis().SetLabelFont(42)
    histo.GetXaxis().SetLabelSize(0.05)
    histo.GetXaxis().SetLabelOffset(0.007)
    histo.GetYaxis().SetTitleFont(42)
    histo.GetYaxis().SetTitleSize(0.05)
    histo.GetYaxis().SetTitleOffset(1.2)
    histo.GetYaxis().SetLabelFont(42)
    histo.GetYaxis().SetLabelSize(0.05)
    histo.GetYaxis().SetLabelOffset(0.007)


def fitFe(rfile,calib=True):
    tf = ROOT.TFile.Open(rfile)
    histo_sig = tf.Get("energyfe_diff" if calib else "integralfe_diff")
    histo_sig.SetMarkerColor(ROOT.kBlack)
    histo_sig.SetLineColor(ROOT.kBlack)
    histo_sig.GetXaxis().SetTitle('energy (keV)' if calib else 'I_{SC} (counts)')
    histo_sig.GetXaxis().SetTitleSize(0.05)
    histo_sig.GetYaxis().SetTitle('superclusters (bkg subtracted)')
    
    c = getCanvas()
    histo_sig.Draw("pe 1")

    par = array( 'd', 6*[0.] )

    xmin1,xmax1 = (3,9) if calib else (1600,3400)
    xmin2,xmax2 = (11,16) if calib else (4400,6000)    
    g1 = ROOT.TF1("g1","gaus",xmin1,xmax1);
    g2 = ROOT.TF1("g2","gaus",xmin2,xmax2);
    total = ROOT.TF1('total','gaus(0)+gaus(3)',xmin1,xmax2)
    total.SetLineColor(ROOT.kBlue+1)
    histo_sig.Fit('g1','RS')
    histo_sig.Fit('g2','R+S')
    mean  = g1.GetParameter(1); mErr = g1.GetParError(1)
    sigma = g1.GetParameter(2); mSigma = g1.GetParError(2)
    lat = ROOT.TLatex()
    unit = 'keV' if calib else 'counts'
    ndigits = 2 if calib else 0
    lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.03)
    lat.DrawLatex(0.55, 0.70, "m_{{1}} = {m:.{nd}f} #pm {em:.{nd}f} {unit}".format(m=mean,em=mErr,unit=unit,nd=ndigits))
    lat.DrawLatex(0.55, 0.65, "#sigma_{{1}} = {s:.{nd}f} #pm {es:.{nd}f} {unit}".format(s=sigma,es=mSigma,unit=unit,nd=ndigits))

    mean  = g2.GetParameter(1); mErr = g2.GetParError(1)
    sigma = g2.GetParameter(2); mSigma = g2.GetParError(2)
    lat.DrawLatex(0.55, 0.50, "m_{{2}} = {m:.{nd}f} #pm {em:.{nd}f} {unit}".format(m=mean,em=mErr,unit=unit,nd=ndigits))
    lat.DrawLatex(0.55, 0.45, "#sigma_{{2}} = {s:.{nd}f} #pm {es:.{nd}f} {unit}".format(s=sigma,es=mSigma,unit=unit,nd=ndigits))



    par1 = g1.GetParameters()
    par2 = g2.GetParameters()

    par[0], par[1], par[2] = par1[0], par1[1], par1[2]
    par[3], par[4], par[5] = par2[0], par2[1], par2[2]
    total.SetParameters( par )
    histo_sig.Fit( total, 'R+' )
    
    c.SaveAs('fe_diff_simplefit.pdf')

def fitDensity(rfile,plotdir):
    tf = ROOT.TFile.Open(rfile)
    histo1 = tf.Get("density")
    histo1.SetMarkerColor(ROOT.kBlack)
    histo1.SetLineColor(ROOT.kBlack)
    histo1.GetXaxis().SetTitle('#delta (photons/pixel)')
    histo1.GetXaxis().SetTitleSize(0.05)
    histo1.GetYaxis().SetTitle('superclusters')
    histo2 = tf.Get("cosm_density")
    histo3 = tf.Get("fe_density")

    c = getCanvas()

    xmin1,xmax1 = (15,25)
    xmin2,xmax2 = (15,25)
    xmin3,xmax3= (15,24)

    if histo3.Integral():
        histo3.SetMarkerColor(ROOT.kGray+2)
        histo3.SetMarkerStyle(ROOT.kFullCircle)
        histo3.SetLineColor(ROOT.kGray+2)
        histo3.Scale(1./histo3.Integral())
        histo3.Sumw2()
        histo3.Draw("pe 1")
        g3 = ROOT.TF1("g3","gaus",xmin3,xmax3);
        g3.SetLineColor(ROOT.kGray+2)
        histo3.Fit('g3','R+S')

    histo1.Sumw2()
    histo1.SetMaximum(0.5)
    histo1.GetYaxis().SetRangeUser(0.,0.5)
    histo1.Scale(1./histo1.Integral())
    histo1.Draw("pe 1 same")

    g1 = ROOT.TF1("g1","gaus",xmin1,xmax1);
    g1.SetLineColor(ROOT.kBlack)
    histo1.Fit('g1','R+S')

    histo2.SetMarkerColor(ROOT.kTeal+4)
    histo2.SetMarkerStyle(ROOT.kFullCircle)
    histo2.SetLineColor(ROOT.kTeal+4)
    histo2.Sumw2()
    histo2.Scale(1./histo2.Integral())
    histo2.Draw("pe 1 same")
    g2 = ROOT.TF1("g2","gaus",xmin2,xmax2);
    g2.SetLineColor(ROOT.kTeal+4)
    histo2.Fit('g2','R+S')


    mean  = g1.GetParameter(1); mErr = g1.GetParError(1)
    sigma = g1.GetParameter(2); mSigma = g1.GetParError(2)
    lat = ROOT.TLatex()
    unit = ''
    ndigits = 2
    lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.03)
    lat.DrawLatex(0.55, 0.80, "m_{{AmBe}} = {m:.{nd}f} #pm {em:.{nd}f} {unit}".format(m=mean,em=mErr,unit=unit,nd=ndigits))
    lat.DrawLatex(0.55, 0.75, "#sigma_{{AmBe}} = {s:.{nd}f} #pm {es:.{nd}f} {unit}".format(s=sigma,es=mSigma,unit=unit,nd=ndigits))

    mean  = g2.GetParameter(1); mErr = g2.GetParError(1)
    sigma = g2.GetParameter(2); mSigma = g2.GetParError(2)
    lat.DrawLatex(0.55, 0.60, "m_{{no source}} = {m:.{nd}f} #pm {em:.{nd}f} {unit}".format(m=mean,em=mErr,unit=unit,nd=ndigits))
    lat.DrawLatex(0.55, 0.55, "#sigma_{{no source}} = {s:.{nd}f} #pm {es:.{nd}f} {unit}".format(s=sigma,es=mSigma,unit=unit,nd=ndigits))

    if histo3.Integral():
        mean  = g3.GetParameter(1); mErr = g3.GetParError(1)
        sigma = g3.GetParameter(2); mSigma = g3.GetParError(2)
        lat.DrawLatex(0.55, 0.40, "m_{{Fe}} = {m:.{nd}f} #pm {em:.{nd}f} {unit}".format(m=mean,em=mErr,unit=unit,nd=ndigits))
        lat.DrawLatex(0.55, 0.35, "#sigma_{{Fe}} = {s:.{nd}f} #pm {es:.{nd}f} {unit}".format(s=sigma,es=mSigma,unit=unit,nd=ndigits))

    for ext in ['png','pdf']:
        c.SaveAs('{pdir}/cosm_scalefit.{ext}'.format(pdir=plotdir,ext=ext))


def makeEff(f1,histo1,f2,histo2,f3=None,histo3=None,plotdir="./",xmax=31):
    # numerator: selected events
    tf1 = ROOT.TFile.Open(f1)
    hpass = tf1.Get(histo1)
    hpass.GetYaxis().SetTitle('efficiency')
    hpass.GetXaxis().SetRangeUser(0,xmax)

    # used to calculate the average efficiency up to a given energy
    averageXmax = 20 # keV
    upperBin = hpass.GetXaxis().FindFixBin(averageXmax)
    passInt = hpass.Integral(0,upperBin+1)
    
    # denominator: all events
    tf2 = ROOT.TFile.Open(f2)
    htotal = tf2.Get(histo2)
    htotal.GetYaxis().SetTitle('signal efficiency ( #varepsilon_{S}^{total} )')
    htotal.GetXaxis().SetRangeUser(0,xmax)

    totalInt = htotal.Integral(1,upperBin+1)

    eff = passInt/totalInt
    errEff = math.sqrt(eff*(1-eff)/totalInt)
    print("Average efficiency from 0 to ",averageXmax," keV = ",eff," +/- ",errEff)
    
    ## default is 68% CL
    col = ROOT.kRed+1
    teffi = ROOT.TEfficiency(hpass,htotal)
    teffi.SetStatisticOption(ROOT.TEfficiency.kFCP);
    teffi.SetFillColorAlpha(col,0.1);
    teffi.SetMarkerStyle(ROOT.kFullSquare);
    teffi.SetMarkerSize(2)
    teffi.SetLineWidth(3)
    teffi.SetMarkerColor(col)
    teffi.SetLineColor(col)
    
    if histo3 and f3:
        # numerator: selected events for the other WP
        tf3 = ROOT.TFile.Open(f3)
        hpass2 = tf3.Get(histo3)
        hpass2.GetYaxis().SetTitle('efficiency')
        hpass2.GetXaxis().SetRangeUser(0,xmax)
        
        ## default is 68% CL
        col2 = ROOT.kAzure-6
        teffi2 = ROOT.TEfficiency(hpass2,htotal)
        teffi2.SetStatisticOption(ROOT.TEfficiency.kFCP);
        teffi2.SetFillColorAlpha(col2,0.1);
        teffi2.SetMarkerStyle(ROOT.kFullCircle);
        teffi2.SetMarkerSize(2)
        teffi2.SetLineWidth(1)
        teffi2.SetMarkerColor(col2);
        teffi2.SetLineColor(col2);

    
    c = getCanvas()
    if histo1.startswith('fe'):
        c.SetLogy()
    teffi.Draw("A3")
    teffi.Draw("same pe1")
    teffi2.Draw("same pe1")
    teffi2.Draw("same3")

    # to change the x-axis range
    ROOT.gPad.Update() 
    teffi.GetPaintedGraph().GetXaxis().SetTitle("E (keV)")
    teffi.GetPaintedGraph().GetXaxis().SetRangeUser(0,xmax)
    teffi.GetPaintedGraph().GetXaxis().SetTitleOffset(1.3)
    teffi.GetPaintedGraph().GetXaxis().SetTitleFont(42)
    teffi.GetPaintedGraph().GetXaxis().SetTitleSize(0.05)

    teffi.GetPaintedGraph().GetYaxis().SetTitle('signal efficiency ( #varepsilon_{S}^{total} )')
    teffi.GetPaintedGraph().GetYaxis().SetRangeUser(0,1)
    teffi.GetPaintedGraph().GetYaxis().SetTitleFont(42)
    teffi.GetPaintedGraph().GetYaxis().SetTitleSize(0.05)
    ROOT.gPad.Update()    

    
    ## add legend
    histos = [teffi,teffi2]
    labels = ['#varepsilon_{B}^{total}=1%','#varepsilon_{B}^{total}=0.1%']
    styles = ['F','F']
    legend = doLegend(histos,labels,styles,corner='BR')
    #legend.Draw('same')
    
    for ext in ['png','pdf']:
        c.SaveAs('{odir}/{var}_effi.{ext}'.format(odir=plotdir,var=histo1,ext=ext))
    

def compareROCs(f1,g1,f2,g2,plotdir):
    tf1 = ROOT.TFile.Open(f1)
    graph1 = tf1.Get(g1)
    graph1.SetName("roc1")

    tf2 = ROOT.TFile.Open(f2)
    graph2 = tf2.Get(g2)
    graph2.SetName("roc2")

    c = getCanvas()
    graph1.SetMarkerStyle(ROOT.kOpenSquare)
    graph1.SetMarkerSize(0.8)
    graph1.SetMarkerColor(ROOT.kGray+2)
    graph1.SetLineColor(ROOT.kGray+2)

    graph2.SetMarkerColor(ROOT.kRed+1)
    graph2.SetLineColor(ROOT.kRed+1)

    graph1.Draw('APC')
    #graph2.Draw('PC')
    graph1.GetXaxis().SetTitle("Signal efficiency ( #varepsilon_{S}^{#delta} )")
    graph1.GetXaxis().SetTitleFont(42)
    graph1.GetXaxis().SetTitleSize(0.05)
    graph1.GetXaxis().SetDecimals()
    graph1.GetYaxis().SetTitle("Background rejection ( 1 - #varepsilon_{B}^{#delta} )")
    graph1.GetYaxis().SetTitleFont(42)
    graph1.GetYaxis().SetTitleSize(0.05)
    graph1.GetYaxis().SetDecimals()
    
    
    lat50 = ROOT.TLatex()
    lat50.SetNDC(); lat50.SetTextFont(42); lat50.SetTextSize(0.03); lat50.SetTextColor(ROOT.kRed+1)
    lat50.DrawLatex(0.6, 0.79, "#delta > 10 photons/pixel") 
    ar50 = ROOT.TArrow(0.52,0.95,0.6,0.95,0.02,"|>");
    ar50.SetLineWidth(2);
    ar50.SetLineColor(ROOT.kRed+1)
    ar50.SetFillColor(ROOT.kRed+1)
    ar50.Draw();

    lat40 = ROOT.TLatex()
    lat40.SetNDC(); lat40.SetTextFont(42); lat40.SetTextSize(0.03); lat40.SetTextColor(ROOT.kAzure+7)
    lat40.DrawLatex(0.3, 0.5, "#delta > 11 photons/pixel") 
    ar40 = ROOT.TArrow(0.39,0.95,0.39,0.6,0.02,"|>");
    ar40.SetLineWidth(2);
    ar40.SetLineColor(ROOT.kAzure+7)
    ar40.SetFillColor(ROOT.kAzure+7)
    ar40.Draw();

    
    for ext in ['png','pdf']:
        c.SaveAs("{plotdir}/comp_roc.{ext}".format(plotdir=plotdir,ext=ext))
        
def fitFeVsPosition(filename,outfile):

    tf_out = ROOT.TFile(outfile,'recreate')
    
    tf_in = ROOT.TFile(filename,'read')
    tree = tf_in.Get('Events')
    
    presel = 'sc_length*0.152<10 && sc_width/sc_length>0.9'

    xmax = 2304
    binsize = 256
    nbinsx = int(xmax/binsize)
    cenb = int(nbinsx/2)+1
    scaleMap = ROOT.TH2F('scale','',nbinsx,0,xmax,nbinsx,0,xmax)
    scaleMap.GetXaxis().SetTitle('x')
    scaleMap.GetYaxis().SetTitle('y')
    resMap = scaleMap.Clone('res')

    counts = ROOT.TH1F('counts','',24,0,5000)
    counts.SetFillStyle(3005)
    counts.SetMarkerStyle(ROOT.kFullCircle)
    counts.SetMarkerSize(2)
    counts.SetLineColor(ROOT.kBlack)
    counts.SetMarkerColor(ROOT.kBlack)

    for ix in range(nbinsx):
        for iy in range(nbinsx):
            hname = 'bin_ix{ix}_iy{iy}'.format(ix=ix,iy=iy)
            print ("Fitting ",hname)
            xlo = scaleMap.GetXaxis().GetBinLowEdge(ix+1)
            xhi = xlo+binsize
            ylo = scaleMap.GetYaxis().GetBinLowEdge(iy+1)
            yhi = ylo+binsize
            fullsel = '{pres} && sc_xmean>{xlo} && sc_xmean<{xhi} && sc_ymean>{ylo} && sc_ymean<{yhi}'.format(pres=presel,xlo=xlo,xhi=xhi,ylo=ylo,yhi=yhi)
            #print ("Fullsel = ",fullsel)
            hist = counts.Clone(hname)
            tree.Draw("sc_integral>>{hname}".format(hname=hname),fullsel)
            #print ("number of entries = ",hist.GetEntries())
            nentries = hist.GetEntries()
            if nentries>20:
                if nentries<50:
                    hist.Rebin(2)
                mean = counts.GetMean()
                rms  = counts.GetRMS()
                f = ROOT.TF1('f','gaus',mean-rms,mean+rms) # there is the double-spot peak, skip it
                f.SetParameter(1,mean);
                f.SetParLimits(1,mean-rms,mean+rms);
                f.SetParameter(2,rms/2.); # there is a tail
                fitr_xmin = mean-rms/2.
                fitr_xmax = mean+rms/2.
                fitRe = hist.Fit(f,'SQ','',fitr_xmin,fitr_xmax)
                rMean  = f.GetParameter(1); rMeanErr = f.GetParError(1)
                rSigma = f.GetParameter(2); rSigmaErr = f.GetParError(1)
                scaleMap.SetBinContent(iy+1,ix+1,rMean); scaleMap.SetBinError(iy+1,ix+1,rMeanErr)
                resMap.SetBinContent(iy+1,ix+1,rSigma/rMean if rMean>0 else 0); resMap.SetBinError(iy+1,ix+1,rSigmaErr/rMean if rMean>0 else 0)
                tf_out.cd()
                hist.Write()
            else:
                scaleMap.SetBinContent(iy+1,ix+1,0)
                resMap.SetBinContent(iy+1,ix+1,0)

    print ("central bin = ",cenb)
    normval = scaleMap.GetBinContent(cenb,cenb)
    for ix in range(nbinsx):
        for iy in range(nbinsx):
            val = scaleMap.GetBinContent(iy+1,ix+1)
            err = scaleMap.GetBinError(iy+1,ix+1)
            scaleMap.SetBinContent(iy+1,ix+1,val/normval)
            scaleMap.SetBinError(iy+1,ix+1,err/normval)

    tf_out.cd()
    scaleMap.Write()
    resMap.Write()

    ROOT.gStyle.SetPaintTextFormat("2.2f");

    cScale = getCanvas('cscale')
    scaleMap.SetMinimum(0.3)
    scaleMap.SetMaximum(1.5)
    scaleMap.Draw('colz')
    scaleMap.Draw('texte same')
    cScale.SaveAs(outfile.replace('.root','_scale.pdf'))
    cScale.Write()
    
    cRes = getCanvas('cres')
    resMap.SetMinimum(0.1)
    resMap.SetMaximum(0.4)
    resMap.Draw('colz')
    resMap.Draw('texte same')
    cRes.SaveAs(outfile.replace('.root','_resolution.pdf'))
    cRes.Write()
    
    tf_out.Close()
    tf_in.Close()
            
def fitIntegral(rfile,xmin,xmax,outfilename,hname='integral_diff'):
    tf = ROOT.TFile.Open(rfile)
    histo_sig = tf.Get(hname)
    histo_sig.SetMarkerColor(ROOT.kBlack)
    histo_sig.SetLineColor(ROOT.kBlack)
    histo_sig.GetXaxis().SetTitle('cluster photons')
    histo_sig.GetXaxis().SetTitleSize(0.05)
    histo_sig.GetYaxis().SetTitle('superclusters (bkg subtracted)')
    
    c = getCanvas()
    histo_sig.Draw("pe 1")

    par = array( 'd', 3*[0.] )

    g1 = ROOT.TF1("g1","gaus",xmin,xmax);
    g1.SetLineColor(ROOT.kBlue+1)
    histo_sig.Fit('g1','RS')
    mean  = g1.GetParameter(1); mErr = g1.GetParError(1)
    sigma = g1.GetParameter(2); sErr = g1.GetParError(2)
    lat = ROOT.TLatex()
    unit = 'counts'
    ndigits = 0
    lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.03)
    lat.DrawLatex(0.55, 0.70, "m_{{1}} = {m:.{nd}f} #pm {em:.{nd}f} {unit}".format(m=mean,em=mErr,unit=unit,nd=ndigits))
    lat.DrawLatex(0.55, 0.65, "#sigma_{{1}} = {s:.{nd}f} #pm {es:.{nd}f} {unit}".format(s=sigma,es=sErr,unit=unit,nd=ndigits))

    c.SaveAs(outfilename)
    return [mean,mErr,sigma,sErr]


def plotVignette1D(histfile="../data/vignette_runs03930to03932.root"):
    f = ROOT.TFile.Open(histfile)
    h_data = f.Get('vign1d')

    c=getCanvas()
    ROOT.gStyle.SetOptStat(0)
    h_data.SetMarkerStyle(ROOT.kFullCircle)
    h_data.SetMarkerSize(1.5)
    h_data.GetXaxis().SetTitle('R (pixels)')
    h_data.GetYaxis().SetTitle('LY / LY_{R=0}')
    formatHisto(h_data)
    h_data.Draw('pe')

    doTinyCmsPrelim("#bf{CYGNO}","(LIME)",lumi=0,textSize=0.05,xoffs=-0.04)
    c.SaveAs('vignette1d.pdf')

def plotPedRms1D(histfile="../pedestals/pedmap_run5857_rebin1.root"):
    f = ROOT.TFile.Open(histfile)
    h_data = f.Get('pedrms')

    c=getCanvas()
    c.SetLeftMargin(0.2)
    ROOT.gStyle.SetOptStat(0)
    h_data.SetFillStyle(1001)
    h_data.SetFillColor(ROOT.kAzure-9)
    h_data.GetXaxis().SetTitle('Pedestal #sigma')
    h_data.GetYaxis().SetTitle('Pixels')
    formatHisto(h_data)
    h_data.GetYaxis().SetTitleOffset(2.1)
    h_data.Draw('hist')

    doTinyCmsPrelim("#bf{CYGNO}","(LIME)",lumi=0,textSize=0.05,xoffs=0.0)
    c.SaveAs('pedrms1d.pdf')

def plotPedHistory(var,tmin=0,tmax=100):
    f = ROOT.TFile.Open('ped%s.root' % var)
    g_data = f.Get('history')
    g_data.SetMarkerStyle(ROOT.kOpenCircle)
    g_data.SetMarkerSize(1.5)
    
    start_time = ROOT.TDatime(2022,6,23,12,0,0).Convert()
    stop_time  = ROOT.TDatime(2022,7,6,12,0,0).Convert()

    g_data.GetXaxis().SetRangeUser(start_time,stop_time)
    g_data.GetXaxis().SetTimeFormat("%d-%m");
    g_data.GetXaxis().SetNdivisions(-502)

    if var=='mean':
        g_data.GetYaxis().SetRangeUser(100.3,100.7)
    else:
        g_data.GetYaxis().SetRangeUser(3.2,3.5)        

    g_data.GetYaxis().SetDecimals()
            
    c=getCanvas()
    c.SetLeftMargin(0.2)
    ROOT.gStyle.SetOptStat(0)
    g_data.GetYaxis().SetTitle('pedestal mean' if var=='mean' else 'pedestal #sigma')
    formatHisto(g_data)
    g_data.GetYaxis().SetTitleOffset(2.1)
    g_data.Draw('AP')

    doTinyCmsPrelim("#bf{CYGNO}","(LIME)",lumi=0,textSize=0.05,xoffs=0)
    c.SaveAs('ped%s.pdf' % var)

    
def compareDeDx(histfile="plots_lnf/2022-12-15-cosmics/plots-bkgonly.root"):
    f = ROOT.TFile.Open(histfile)
    h_data = f.Get('dedx_background')
    h_data.SetTitle('')
    h_toymc_arr = [0.0, 0.0, 0.0, 0.0, 3.0, 6.0, 4.0, 9.0, 18.0, 26.0, 74.0, 122.0, 196.0, 309.0, 476.0, 619.0, 838.0, 1170.0, 1494.0, 1884.0, 2224.0, 2735.0, 3054.0, 3359.0, 3378.0, 3605.0, 3395.0, 3042.0, 2932.0, 2673.0, 2362.0, 2037.0, 1656.0, 1458.0, 1232.0, 1012.0, 784.0, 608.0, 586.0, 440.0, 341.0, 342.0, 245.0, 244.0, 216.0, 158.0, 160.0, 144.0, 143.0, 108.0, 86.0, 79.0, 80.0, 56.0, 79.0, 56.0, 42.0, 34.0, 36.0, 36.0, 34.0, 28.0, 18.0, 27.0, 22.0, 15.0, 25.0, 16.0, 23.0, 16.0, 10.0, 11.0, 15.0, 11.0, 6.0, 13.0, 7.0, 6.0, 7.0, 11.0, 3.0, 7.0, 12.0, 7.0, 7.0, 6.0, 6.0, 5.0, 2.0, 5.0, 3.0, 3.0, 3.0, 2.0, 6.0, 3.0, 2.0, 3.0, 2.0, 3.0]
    binsize = 12./100. # keV/cm

    h_toymc = h_data.Clone('dedx_toymc')
    h_toymc.SetFillStyle(1001)
    for ibin,val in enumerate(h_toymc_arr):
        h_toymc.SetBinContent(ibin+1,val)

    h_toymc.Sumw2()
    h_toymc.Scale(h_data.Integral()/h_toymc.Integral())

    c=getCanvas()
    ROOT.gStyle.SetOptStat(0)
    c.SetLogy()
    h_data.SetMarkerStyle(ROOT.kFullCircle)
    h_data.SetMarkerSize(1.5)
    h_toymc.SetFillColor(ROOT.kAzure+6)
    h_toymc.GetYaxis().SetTitle('Superclusters')
    formatHisto(h_toymc)
    h_toymc.Draw('hist')
    h_data.Draw('pe same')

    ## add legend
    histos = [h_data,h_toymc]
    labels = ['data','simulation']
    styles = ['pe','F']
    legend = doLegend(histos,labels,styles,corner='TR')
    legend.Draw()

    doTinyCmsPrelim("#bf{CYGNO}","(LIME)",lumi=0,textSize=0.05,xoffs=-0.04)

    c.SaveAs('dedx_cosmics_datamc.pdf')
    
### this is meant to be run on top of ROOT files produced by simple_plots, not on the trees
if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('', '--make'   , type='string'       , default='fitfe' , help='run ambe_miscellanea.py (options = fitfe, efficiency, tworocs)')
    parser.add_option('', '--outdir' , type='string'       , default='./'    , help='output directory with directory structure and plots')
    parser.add_option('', '--source' , type='string'       , default='ambe'  , help='in case of efficiency plotting, make it for fe/ambe')
    (options, args) = parser.parse_args()

    ROOT.gROOT.ProcessLine(".x tdrstyle.cc")

    if  options.make == 'fitfe':
        fitFe('plots/ambe/clusters_3sourcesNloCalNeutronsFex1_2020_05_05/energy.root')
    elif options.make == 'fitfeuncalib':
        fitFe('plots/ambe/clusters_3sourcesNloCalNeutronsFex1_2020_05_05/integral.root',calib=False)
    elif options.make == 'fitdensity':
        fitDensity('plots/ambeV3_cosmsel_compareCosm_01-12-2020/density.root',options.outdir)

    ## usages:
    # AmBe efficiency:> python ambe_miscellanea.py --make efficiency --source ambe --outdir './'
    # Fe55 efficiency:> python ambe_miscellanea.py --make efficiency --source fe --outdir './'
    elif options.make == 'efficiency':
        var = 'energy' if options.source=='fe' else 'energyFull'
        prefix = '' if options.source=='ambe' else ('fe_' if options.source=='fe' else 'cosm_')
        ## OLD binning
        # makeEff('plots/ambe/clusters_3sourcesNloCalNeutronsDensityGt11_2020_05_06/'+var+'.root',prefix+var,
        #         'plots/ambe/clusters_3sourcesNloCalNeutrons_2020_05_05/'+var+'.root',           prefix+var,
        #         options.outdir)
        ## PAPER binning
        # makeEff('plots/ambe/clusters_3sources_WP50Paper_2020_05_29/'+var+'.root',prefix+var,
        #         'plots/ambe/clusters_3sources_FullSelPaper_2020_05_29/'+var+'.root',prefix+var,
        #         'plots/ambe/clusters_3sources_WP40Paper_2020_05_29/'+var+'.root',prefix+var,
        #         options.outdir)
        ## LIME
        makeEff('ambeplots/2021-03-29-ambeV6_cosmveto_WPbkg10em2/'+var+'.root',prefix+var,
                'ambeplots/2021-03-29-ambeV6_cosmveto/'+var+'.root',prefix+var,
                'ambeplots/2021-03-31-ambeV6_cosmveto_WPDNNbkg10em2/'+var+'.root',prefix+var,
                options.outdir)
                
    elif options.make == 'tworocs':
        compareROCs('plots/ambe/clusters_3sources_FullSelPaper_2020_05_29/density_roc.root','Graph',
                    'plots/ambe/clusters_3sources_FullSelAndPMTCutPaper_2020_05_29/density_roc.root','Graph',
                    options.outdir)
    elif options.make == 'scalevspos':
        fitFeVsPosition('~/Work/data/cygnus/RECO/lime2020/v2/fe_lime.root','fe_novign.root')
        fitFeVsPosition('~/Work/data/cygnus/RECO/lime2020/v4/fe55_runs3686to3691_v4.root','fe_v4.root')
        fitFeVsPosition('~/Work/data/cygnus/RECO/lime2020/v4/fe55_runs3686to3691_v4_OptVignetting.root','fe_v4_OptVignetting.root')

    elif options.make == 'linearity':
        #filepatt = 'ambeplots/2021-07-21-xrays-%s/integral.root'
        filepatt = 'ambeplots/2021-10-12-xrays-%s/integral.root'
        sources = ['Ti','Fe','Cu','Rb','Mo','Ag','Ba','Tb']
        energies = [4.51,5.9,8,15,18,23.5,34.4,47]
        files =  [filepatt%s for s in sources]
        # old data
        # ranges = [(5e3,1e4),(1e4,1.5e4),(1.3e4,2e4),(1.5e4,3e4),(2.5e4,4e4)]
        # July 2021 data
        ranges = [(6.5e3,8.5e3),(8e3,1.2e4),(1e4,1.4e4),(1.7e4,2.5e4),(1.8e4,2.9e4),(2.2e4,3.8e4),(3.5e4,6.2e4),(5.2e4,8.0e4)]
        pars = []
        for i,f in enumerate(files):
            pars.append(fitIntegral(f,ranges[i][0],ranges[i][1],sources[i]+'.pdf',hname='integral_diff'))

        fout = ROOT.TFile.Open("linearity.root","recreate")
        resp = ROOT.TGraphErrors(len(sources))
        reso = ROOT.TGraphErrors(len(sources))
        for i in range(len(sources)):
            resp.SetPoint(i,energies[i],pars[i][0])
            resp.SetPointError(i,0,pars[i][1])
            reso.SetPoint(i,energies[i],100*pars[i][2]/pars[i][0])
            reso.SetPointError(i,0,100*pars[i][3]/pars[i][0])


        ## response linearity
        c = getCanvas()
        resp.SetTitle('')
        resp.Draw("AP")
        resp.SetMarkerStyle(ROOT.kFullCircle)
        resp.SetMarkerSize(2)
        resp.GetXaxis().SetLimits(0,50)
        resp.GetXaxis().SetTitle("Energy (keV)")
        resp.GetYaxis().SetLimits(0,8e4)
        resp.GetYaxis().SetRangeUser(0,8e4)
        resp.GetYaxis().SetTitle("Cluster integral (counts)")

        resp.Fit("pol1")
        
        line = ROOT.TLine(0,0,50,50*30105/23.5) # Calibration from Ag
        line.SetLineColor(ROOT.kBlue)
        line.SetLineStyle(ROOT.kDashed)
        line.Draw()

        fout.cd()
        resp.Write()
        c.Write()
        c.SaveAs("linearity.pdf")
        fout.Close()
        
        # energy resolution
        c = getCanvas()
        reso.SetTitle('')
        reso.Draw("AP")
        reso.SetMarkerStyle(ROOT.kFullCircle)
        reso.SetMarkerSize(2)
        reso.GetXaxis().SetLimits(0,50)
        reso.GetXaxis().SetTitle("Energy (keV)")
        reso.GetYaxis().SetLimits(0,30)
        reso.GetYaxis().SetRangeUser(0,30)
        reso.GetYaxis().SetTitle("#delta E/E (%)")        
        c.SaveAs("resolution.pdf")

    elif options.make == 'dedx-comp':
        compareDeDx()

    elif options.make == 'plotvignette':
        plotVignette1D()

    elif options.make == 'pedrms1d':
        plotPedRms1D()

    elif options.make == 'pedhistory':
        plotPedHistory('mean')
        plotPedHistory('rms')
        
    else:
        print ("make ",options.make," not implemented.")

    
