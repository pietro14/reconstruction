import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
import math, os, joblib

from gbr_trainer import GBRLikelihoodTrainer,fill_hist

Zsteps = [5,15,25,36,48]

def doSpam(text,x1,y1,x2,y2,align=12,fill=False,textSize=0.033,_noDelete={}):
    cmsprel = ROOT.TPaveText(x1,y1,x2,y2,"NDC");
    cmsprel.SetTextSize(textSize);
    cmsprel.SetFillColor(0);
    cmsprel.SetFillStyle(1001 if fill else 0);
    cmsprel.SetLineStyle(2);
    cmsprel.SetLineColor(0);
    cmsprel.SetTextAlign(align);
    cmsprel.SetTextFont(42);
    cmsprel.AddText(text);
    cmsprel.Draw("same");
    _noDelete[text] = cmsprel; ## so it doesn't get deleted by PyROOT
    return cmsprel

def doTinyCmsPrelim(textLeft="#bf{CYGNO}",textRight="(LIME)",hasExpo=False,textSize=0.033,lumi=None, xoffs=0, options=None, doWide=False):
    if textLeft  == "_default_": textLeft  = "#bf{CYGNO}"
    if textRight == "_default_": textRight = "(LIME)"
    if lumi      == None       : lumi      = 1
    if   lumi > 3.54e+1: lumitext = "%.0f fb^{-1}" % lumi
    elif lumi > 3.54e+0: lumitext = "%.1f fb^{-1}" % lumi
    elif lumi > 3.54e-1: lumitext = "%.2f fb^{-1}" % lumi
    elif lumi > 3.54e-2: lumitext = "%.0f pb^{-1}" % (lumi*1000)
    elif lumi > 3.54e-3: lumitext = "%.1f pb^{-1}" % (lumi*1000)
    else               : lumitext = "%.2f pb^{-1}" % (lumi*1000)
    lumitext = "%.1f fb^{-1}" % lumi
    textLeft = textLeft.replace("%(lumi)",lumitext)
    textRight = textRight.replace("%(lumi)",lumitext)
    if textLeft not in ['', None]:
        doSpam(textLeft, (.28 if hasExpo else 0.2 if doWide else .12)+xoffs, .94, .60+xoffs, .94, align=12, textSize=textSize)
    if textRight not in ['', None]:
        doSpam(textRight,(0.5 if doWide else .55)+xoffs, .94, .82+xoffs if doWide else .91+xoffs, .94, align=32, textSize=textSize)

def getCanvas(name='c'):

    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ROOT.gStyle.SetNumberContours(51)
    ROOT.gErrorIgnoreLevel = 100
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas(name,'',1200,1200)
    lMargin = 0.14
    rMargin = 0.10
    bMargin = 0.15
    tMargin = 0.10
    c.SetLeftMargin(lMargin)
    c.SetRightMargin(rMargin)
    c.SetTopMargin(tMargin)
    c.SetBottomMargin(bMargin)
    c.SetFrameBorderMode(0);
    c.SetBorderMode(0);
    c.SetBorderSize(0);
    return c

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


def getFriend(inputfile):
    friendfile = "{dirn}/{base}_Friend.root".format(dirn=os.path.dirname(inputfile),
                                                    base=os.path.basename(inputfile).split('.root')[0])
    return friendfile

def response_vsrun(inputfile,paramsfile,Z,panda=None):
    GBR = GBRLikelihoodTrainer(paramsfile)
    eps=0.1
    zSel = {'sc_truez': (Z-eps,Z+eps)}
    X,y = GBR.get_dataset(inputfile,getFriend(inputfile),loadPanda=panda,addCuts=zSel)
    filename = "gbrLikelihood_mse.sav"
    model = joblib.load(filename)
    y_pred = model.predict(X)

    h = ROOT.TH1F('h','',800,0,2)
    histos = {}

    histos['uncorr'] = h.Clone('h_uncorr')
    histos['regr'] = h.Clone('h_regr')

    fill_hist(histos['uncorr'],y)
    fill_hist(histos['regr'],y_pred)

    return (histos['uncorr'],histos['regr'])
    
    
def response2D(inputfile,paramsfile,panda=None):

    friendfile = "{dirn}/{base}_Friend.root".format(dirn=os.path.dirname(inputfile),
                                                    base=os.path.basename(inputfile).split('.root')[0])
    GBR = GBRLikelihoodTrainer(paramsfile)
    X,y = GBR.get_dataset(inputfile,getFriend(inputfile),loadPanda=panda)

    xy_indices = []
    vars = GBR.variables()
    for i,v in enumerate(vars):
        if v=='sc_xmean': xy_indices.append(i)
        if v=='sc_ymean': xy_indices.append(i)

    xy = X[:,xy_indices[0]:xy_indices[1]+1]
    print("xy = ",xy)
    print("y = ",y)

    energy_2D = ROOT.TH3D('energy_2D','',144,0,2304,144,0,2304,200,0,2)
    energy_1D = ROOT.TH2D('energy_1D','',20,0,2304/math.sqrt(2.),70,0,2)
    
    filename = "gbrLikelihood_mse.sav"
    model = joblib.load(filename)
    y_pred = model.predict(X)

    energy_2Ds = {}
    energy_1Ds = {}

    energy_2Ds['uncorr'] = energy_2D.Clone('energy_2D_uncorr')
    energy_2Ds['regr'] = energy_2D.Clone('energy_2D_regr')

    energy_1Ds['uncorr'] = energy_1D.Clone('energy_1D_uncorr')
    energy_1Ds['regr'] = energy_1D.Clone('energy_1D_regr')

    hresp = ROOT.TH1F('hresp','',50,0.,2.0)
    center = 2304/2.
    
    for i in range(len(y)):
        #print("ev {iev} has x,y=({x},{y}) and y = {z}".format(iev=i,x=xy[i][0],y=xy[i][1],z=y[i]))
        energy_2Ds['uncorr'].Fill(xy[i][1],xy[i][0],y[i])
        energy_2Ds['regr'].Fill(xy[i][1],xy[i][0],y_pred[i])
        energy_1Ds['uncorr'].Fill(math.hypot(xy[i][1]-center,xy[i][0]-center),y[i])
        energy_1Ds['regr'].Fill(math.hypot(xy[i][1]-center,xy[i][0]-center),y_pred[i])

    energy_2D_mode = ROOT.TH2D('energy_2D_mode','',144,0,2304,144,0,2304)
    energy_1D_mode = ROOT.TH1D('energy_1D_mode','',20,200,2100/math.sqrt(2.))
    energy_2D_mode.GetXaxis().SetTitle("ix")
    energy_2D_mode.GetYaxis().SetTitle("iy")
    energy_2D_mode.GetZaxis().SetTitle("mode")
    energy_2D_mode.GetZaxis().SetRangeUser(0.2,1.4)
    
    energy_2D_modes = {}
    energy_2D_modes['uncorr'] = energy_2D_mode.Clone('energy_2D_mode_uncorr')
    energy_2D_modes['regr'] = energy_2D_mode.Clone('energy_2D_mode_regr')

    ROOT.gStyle.SetPaintTextFormat("1.2f");

    c = getCanvas('c')
    for k,h3D in energy_2Ds.items():
        for ix in range(1,h3D.GetNbinsX()+1):
            for iy in range(1,h3D.GetNbinsY()+1):
                # hz = ROOT.TH1F('hz','',200,0,2)
                # for iz in range(1,h3D.GetNbinsZ()+1):
                #     hz.SetBinContent(iz,h3D.GetBinContent(ix,iy,iz))
                # zmean = hz.GetMean()
                # print ("zmean = ",zmean)
                # energy_2D_modes[k].SetBinContent(ix,iy,zmean)
                                                 
                # h_integral = sum([h3D.GetBinContent(ix,iy,iz) for iz in range(1,h3D.GetNbinsZ()+1)])
                # h_runint = 0; zbinmedian = -1
                # for iz in range(1,h3D.GetNbinsZ()+1):
                #     h_runint = h_runint + h3D.GetBinContent(ix,iy,iz)
                #     if h_runint > 0.5 * h_integral:
                #         zbinmedian = iz
                #         break

                maxZ = -1; zbinmax = -1
                for iz in range(1,h3D.GetNbinsZ()+1):
                    if h3D.GetBinContent(ix,iy,iz) > maxZ:
                        maxZ = h3D.GetBinContent(ix,iy,iz)
                        zbinmax = iz
                energy_2D_modes[k].SetBinContent(ix,iy,h3D.GetZaxis().GetBinCenter(zbinmax))
        energy_2D_modes[k].Draw("colz") # text45")
        c.SaveAs("energy_2D_{name}.png".format(name=k))

    energy_2D_modes['ratio'] = energy_2D_modes['regr'].Clone('energy_2D_mode_ratio')
    energy_2D_modes['ratio'].Divide(energy_2D_modes['uncorr'])
    energy_2D_modes['ratio'].Draw("colz") # text45")
    c.SaveAs("energy_2D_ratio.png")
    
    
    # energy_1Ds['uncorr'].SetMarkerColor(ROOT.kBlack)
    # energy_1Ds['regr'].SetMarkerColor(ROOT.kRed)
    # energy_1Ds['uncorr'].Draw("pe1")
    # energy_1Ds['regr'].Draw("pe1 same")
    # c.SaveAs("energy_1D.png")

def makeResponseHistos(treeFile,params_txt,outfileHistos="response_histos.root",panda=None):
    outfile = ROOT.TFile.Open(outfileHistos,"recreate")
    response2D(treeFile,params_txt,panda=panda)
    outfile.cd()
    for i,dist in enumerate(Zsteps):
        print("filling graph with point ",i," and Z = ",dist)
        histos = response_vsrun(args[0],params_txt,dist,panda=panda)
        for ih,h in enumerate(histos):
            h.Draw("hist")
            suff = 'uncorr' if ih==0 else 'regr'
            name = "resp_{suf}_{dist}".format(suf=suff,dist=dist)
            h.SetName(name)
            outfile.cd()
            h.Write()                        
    outfile.Close()

def fitResponseHisto(histo,xmin=0.3,xmax=1.3,rebin=4,marker=ROOT.kFullCircle,color=ROOT.kRed+1):

    histo.Rebin(rebin)
    work = ROOT.RooWorkspace()
    work.factory('CBShape::cb(x[{xmin},{xmax}],mean[0.5,1.4],sigma[0.05,0.01,0.50],alpha[1,0.05,10],n[5,1,10])'.format(xmin=xmin,xmax=xmax))
    work.Print()
    
    x = work.var('x')

    histoname = histo.GetName()
    rooData = ROOT.RooDataHist(histoname,histoname,ROOT.RooArgList(work.var("x")),histo)
    getattr(work,'import')(rooData)

    frame = x.frame()
    frame.SetTitle('')

    # fit histogram
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(marker),ROOT.RooFit.MarkerSize(2))
    pdf = work.pdf('cb')
    pdf.fitTo(rooData,ROOT.RooFit.Save())#,ROOT.RooFit.PrintLevel(-1))
    pdf.plotOn(frame,ROOT.RooFit.LineColor(color))
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(marker),ROOT.RooFit.MarkerSize(2))

    frame.GetXaxis().SetTitleFont(42)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetLabelFont(42)
    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetLabelOffset(0.007)
    frame.GetYaxis().SetTitleFont(42)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.GetYaxis().SetLabelFont(42)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetLabelOffset(0.007)
    frame.GetYaxis().SetTitle("Events")
    frame.GetXaxis().SetTitle("E/E^{mpv}")
    frame.GetXaxis().SetNdivisions(510)

    m = work.var('mean').getVal()
    s = work.var('sigma').getVal()

    em = work.var('mean').getError()
    es = work.var('sigma').getError()

    return (frame,m,s,em,es)

def fitAllResponses(inputfile="response_histos.root"):

    limits = {5 : (0.3,1.1),
              15  : (0.4,1.4),
              25  : (0.5,1.6),
              36  : (0.6,1.6),
              48  : (0.6,1.7),
              }

    dummy_uncorr = ROOT.TH1F("dummy_uncorr","",1,0,1)
    dummy_regr = ROOT.TH1F("dummy_regr","",1,0,1)

    dummy_uncorr.SetMarkerStyle(ROOT.kOpenSquare)
    dummy_uncorr.SetLineColor(ROOT.kRed+2)
    dummy_uncorr.SetMarkerSize(2)
    
    dummy_regr.SetMarkerStyle(ROOT.kFullCircle)
    dummy_regr.SetLineColor(ROOT.kBlue)
    dummy_regr.SetMarkerSize(2)

    results = {}
    tf = ROOT.TFile(inputfile,"read")
    for i,dist in enumerate(Zsteps):
        c = getCanvas()
        for corr in ["regr","uncorr"]:
            hname = "resp_{corr}_{dist}".format(corr=corr,dist=dist)
            histo = tf.Get(hname)
            xmin = limits[dist][0]; xmax=limits[dist][1]
            if corr=="uncorr":
                marker = ROOT.kOpenSquare
                linecol = ROOT.kRed+2
            else:
                marker = ROOT.kFullCircle
                linecol = ROOT.kBlue
            (frame,mean,sigma,meanerr,sigmaerr)=fitResponseHisto(histo,xmin,xmax,8,marker,linecol)
            frame.Draw("same")
            results[(corr,dist)] = (mean,sigma,meanerr,sigmaerr)
            
        responses = [dummy_uncorr,dummy_regr]
        titles = ['E_{rec}','E']
        styles = ['pl','pl']    
        legend = doLegend(responses,titles,styles,corner="TL")    

        doTinyCmsPrelim(hasExpo = False,textSize=0.04, options=None,doWide=False)
        for ext in ['pdf','png','root']:
            c.SaveAs('respcomp_{dist}_fit.{ext}'.format(dist=dist,ext=ext))
    return results

def graphVsR(typeInput='histo',histoFile="response_histos.root"):

    resp_vs_z_uncorr = ROOT.TGraphErrors(len(Zsteps))
    reso_vs_z_uncorr = ROOT.TGraph(len(Zsteps))
    resp_vs_z_regr = ROOT.TGraphErrors(len(Zsteps))
    reso_vs_z_regr = ROOT.TGraph(len(Zsteps))

    reso_vs_z_uncorr_fullrms = ROOT.TGraph(len(Zsteps))
    reso_vs_z_regr_fullrms = ROOT.TGraph(len(Zsteps))

    resp_vs_z_uncorr.GetXaxis().SetTitle("z (cm)")
    resp_vs_z_uncorr.GetYaxis().SetTitle("E/E^{peak}")
    resp_vs_z_uncorr.SetTitle("")

    reso_vs_z_uncorr.GetXaxis().SetTitle("z (cm)")
    reso_vs_z_uncorr.GetYaxis().SetTitle("E/E^{peak} resolution")
    reso_vs_z_uncorr.SetTitle("")

    if typeInput=="histo":
        tf = ROOT.TFile(histoFile,"read")
        for i,dist in enumerate(Zsteps):
            huncorr = tf.Get("resp_uncorr_{dist}".format(dist=dist))
            hregr   = tf.Get("resp_regr_{dist}".format(dist=dist))
            resp_vs_z_uncorr.SetPoint(i,dist,huncorr.GetMean())
            resp_vs_z_uncorr.SetPointError(i,0,huncorr.GetMeanError())

            resp_vs_z_regr.SetPoint(i,dist,hregr.GetMean())
            resp_vs_z_regr.SetPointError(i,0,hregr.GetMeanError())
        
            reso_vs_z_uncorr.SetPoint(i,dist,huncorr.GetRMS()/huncorr.GetMean())
            reso_vs_z_regr.SetPoint(i,dist,hregr.GetRMS()/hregr.GetMean())
        tf.Close()

    elif typeInput=="fit":
        results = fitAllResponses()
        tf = ROOT.TFile(histoFile,"read")
        for i,dist in enumerate(Zsteps):

            mean,sigma,meanerr,sigmaerr = results[('uncorr',dist)]
            resp_vs_z_uncorr.SetPoint(i,dist,mean)
            resp_vs_z_uncorr.SetPointError(i,0,meanerr)
            reso_vs_z_uncorr.SetPoint(i,dist,sigma/mean)

            mean,sigma,meanerr,sigmaerr = results[('regr',dist)]
            resp_vs_z_regr.SetPoint(i,dist,mean)
            resp_vs_z_regr.SetPointError(i,0,meanerr)
            reso_vs_z_regr.SetPoint(i,dist,sigma/mean)

            huncorr = tf.Get("resp_uncorr_{dist}".format(dist=dist))
            hregr   = tf.Get("resp_regr_{dist}".format(dist=dist))
            reso_vs_z_uncorr_fullrms.SetPoint(i,dist,huncorr.GetRMS()/huncorr.GetMean())
            reso_vs_z_regr_fullrms.SetPoint(i,dist,hregr.GetRMS()/hregr.GetMean())
        tf.Close()

    c = getCanvas('c')

    # response
    resp_vs_z_uncorr.SetMarkerStyle(ROOT.kFullCircle)
    resp_vs_z_regr.SetMarkerStyle(ROOT.kFullCircle)
    resp_vs_z_uncorr.SetMarkerSize(2)
    resp_vs_z_regr.SetMarkerSize(2)
    resp_vs_z_uncorr.SetMarkerColor(ROOT.kRed)
    resp_vs_z_uncorr.SetLineColor(ROOT.kRed)
    resp_vs_z_regr.SetMarkerColor(ROOT.kBlack)
    resp_vs_z_uncorr.Draw("Ape")
    resp_vs_z_regr.Draw("pe")

    responses = [resp_vs_z_uncorr,resp_vs_z_regr]
    titles = ['uncorrected','regression']
    styles = ['pl','pl']
    
    resp_vs_z_uncorr.GetYaxis().SetRangeUser(0.6,1.4)
    legend = doLegend(responses,titles,styles,corner="TL")    
    for ext in ['pdf','png']:
        c.SaveAs("response.%s"%ext)

    # resolution
    reso_vs_z_uncorr.SetMarkerStyle(ROOT.kFullCircle)
    reso_vs_z_regr.SetMarkerStyle(ROOT.kFullCircle)
    reso_vs_z_uncorr.SetMarkerColor(ROOT.kRed)
    reso_vs_z_regr.SetMarkerColor(ROOT.kBlack)
    reso_vs_z_uncorr.SetMarkerSize(2)
    reso_vs_z_regr.SetMarkerSize(2)
    reso_vs_z_uncorr.Draw("Ape")
    reso_vs_z_regr.Draw("pe")
    legCols = 1
    if typeInput=='fit':
        reso_vs_z_uncorr_fullrms.SetMarkerStyle(ROOT.kOpenSquare)
        reso_vs_z_regr_fullrms.SetMarkerStyle(ROOT.kOpenSquare)
        reso_vs_z_uncorr_fullrms.SetMarkerColor(ROOT.kRed)
        reso_vs_z_regr_fullrms.SetMarkerColor(ROOT.kBlack)
        reso_vs_z_uncorr_fullrms.SetMarkerSize(2)
        reso_vs_z_regr_fullrms.SetMarkerSize(2)
        responses = [resp_vs_z_uncorr,resp_vs_z_regr,reso_vs_z_uncorr_fullrms,reso_vs_z_regr_fullrms]
        titles = ['uncorrected (#sigma_{G})','regression (#sigma_{G})','uncorrected (rms)','regression (rms)']
        styles = ['pl','pl','pl','pl']
        reso_vs_z_uncorr_fullrms.Draw("pe")
        reso_vs_z_regr_fullrms.Draw("pe")
        legCols=2
        
    reso_vs_z_uncorr.GetYaxis().SetRangeUser(0.0,0.4)
    legend = doLegend(responses,titles,styles,corner="TL",textSize=0.03,legWidth=0.8,nColumns=legCols)
    for ext in ['pdf','png']:
        c.SaveAs("resolution.%s"%ext)

            
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog input.root params_gbrtrain.txt [opts] ')
    parser.add_option(        '--loadPanda', dest='loadPanda', default=None, type='string', help='file with regression data as panda dataframe to be loaded instead of the full ROOT files')
#    parser.add_option(        '--truthmap', dest='truthmap', default="../postprocessing/data/fitm-zhvscans-ext.pkl", type='string', help='pickle file with the results of the fit: table with HV,z,LY,sigma(LY)')
    (options, args) = parser.parse_args()

    makeResponseHistos(args[0],args[1],panda=options.loadPanda)
    graphVsR("fit")
