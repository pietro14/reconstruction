import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
from root_numpy import tree2array,fill_hist,fill_profile
from gbr_trainer import GBRLikelihoodTrainer, getCanvas, doLegend
import math, joblib

# map run - z distance
#          no source        8.7       8.4        8.1k       5k
#runs = {4433: 50, 4441: 46, 4448: 36, 4455: 26, 4463: 6}
#runs = {4441: 46, 4448: 36, 4455: 26}
runs = {4120: 46,
        4128: 41,
        4136: 36,
        4144: 31,
        4152: 26,
        4160: 21,
        4168: 16,
        4176: 11,
        4184: 6.5}


def response_vsrun(inputfile,params,runN):
    params['selection'] = "({sel}) && run=={run}".format(sel=params['selection'],run=runN)
    print ("New selection SEL = {sel}".format(sel=params['selection']))
    GBR = GBRLikelihoodTrainer(params)
    X,y = GBR.get_dataset(inputfile)
    filename = "gbrLikelihood_mse.sav"
    model = joblib.load(filename)
    y_pred = model.predict(X)

    h = ROOT.TH1F('h','',400,0,2)
    histos = {}

    histos['uncorr'] = h.Clone('h_uncorr')
    histos['regr'] = h.Clone('h_regr')

    fill_hist(histos['uncorr'],y)
    fill_hist(histos['regr'],y_pred)

    return (histos['uncorr'],histos['regr'])
    
    
def response2D(inputfile,params):

    GBR = GBRLikelihoodTrainer(params)
    X,y = GBR.get_dataset(inputfile)

    xy_indices = []
    vars = params["inputs"].split("|")
    for i,v in enumerate(vars):
        if v=='sc_xmean': xy_indices.append(i)
        if v=='sc_ymean': xy_indices.append(i)

    xy = X[:,xy_indices[0]:xy_indices[1]+1]
    print("xy = ",xy)
    print("y = ",y)

    energy_2D = ROOT.TH3D('energy_2D','',10,0,2304,10,0,2304,200,0,2)
    energy_1D = ROOT.TH2D('energy_1D','',20,0,2304/math.sqrt(2.),70,0,2)
    
    filename = "gbrLikelihood_q0.50.sav"
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

    energy_2D_mode = ROOT.TH2D('energy_2D_mode','',10,0,2304,10,0,2304)
    energy_1D_mode = ROOT.TH1D('energy_1D_mode','',20,200,2100/math.sqrt(2.))
    energy_2D_mode.GetXaxis().SetTitle("ix")
    energy_2D_mode.GetYaxis().SetTitle("iy")
    energy_2D_mode.GetZaxis().SetTitle("mode")
    energy_2D_mode.GetZaxis().SetRangeUser(0.8,1.3)
    
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
        energy_2D_modes[k].Draw("colz text45")
        c.SaveAs("energy_2D_{name}.png".format(name=k))

    energy_2D_modes['ratio'] = energy_2D_modes['regr'].Clone('energy_2D_mode_ratio')
    energy_2D_modes['ratio'].Divide(energy_2D_modes['uncorr'])
    energy_2D_modes['ratio'].Draw("colz text45")
    c.SaveAs("energy_2D_ratio.png")
    
    
    # energy_1Ds['uncorr'].SetMarkerColor(ROOT.kBlack)
    # energy_1Ds['regr'].SetMarkerColor(ROOT.kRed)
    # energy_1Ds['uncorr'].Draw("pe1")
    # energy_1Ds['regr'].Draw("pe1 same")
    # c.SaveAs("energy_1D.png")

def makeResponseHistos(treeFile,params_txt,outfileHistos="response_histos.root"):
    outfile = ROOT.TFile.Open(outfileHistos,"recreate")
    config = open(params_txt, "r")
    params = eval(config.read())
    config.close()
    response2D(treeFile,params)
    outfile.cd()
    for i,v in enumerate(runs.items()):
        config = open(args[1], "r")
        params = eval(config.read()) # re-read to reload run
        config.close()
        run,dist = v
        print("filling graph with point ",i," and run = ",run," corresponding to distance ",dist)
        histos = response_vsrun(args[0],params,run)
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
    work.factory('CBShape::cb(x[{xmin},{xmax}],mean[0.7,1.4],sigma[0.05,0.1,0.30],alpha[1,0.1,10],n[5,1,10])'.format(xmin=xmin,xmax=xmax))
    work.Print()
    
    x = work.var('x')

    histoname = histo.GetName()
    rooData = ROOT.RooDataHist(histoname,histoname,ROOT.RooArgList(work.var("x")),histo)
    getattr(work,'import')(rooData)

    frame = x.frame()
    frame.SetTitle('')

    # fit histogram
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(marker))
    pdf = work.pdf('cb')
    pdf.fitTo(rooData,ROOT.RooFit.Save())#,ROOT.RooFit.PrintLevel(-1))
    pdf.plotOn(frame,ROOT.RooFit.LineColor(color))
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(marker))

    frame.GetYaxis().SetTitle("superclusters")
    frame.GetXaxis().SetTitle("E/E^{raw}_{peak}")
    frame.GetXaxis().SetTitleOffset(1.2)
    
    m = work.var('mean').getVal()
    s = work.var('sigma').getVal()

    em = work.var('mean').getError()
    es = work.var('sigma').getError()

    return (frame,m,s,em,es)

def fitAllResponses(inputfile="response_histos.root"):

    limits = {6.5 : (0.3,1.3),
              11  : (0.3,1.5),
              16  : (0.3,1.5),
              21  : (0.3,1.6),
              26  : (0.3,1.6),
              31  : (0.3,1.5),
              36  : (0.3,1.5),
              41  : (0.3,1.5),
              46  : (0.3,1.4),
              }

    results = {}
    tf = ROOT.TFile(inputfile,"read")
    for i,v in enumerate(runs.items()):
        run,dist = v
        c = getCanvas()
        for corr in ["regr","uncorr"]:
            hname = "resp_{corr}_{dist}".format(corr=corr,dist=dist)
            histo = tf.Get(hname)
            if corr=="uncorr":
                marker = ROOT.kOpenSquare
                linecol = ROOT.kRed+2
            else:
                marker = ROOT.kFullCircle
                linecol = ROOT.kBlue
            (frame,mean,sigma,meanerr,sigmaerr)=fitResponseHisto(histo,limits[dist][0],limits[dist][1],8,marker,linecol)
            frame.Draw("same")
            results[(corr,dist)] = (mean,sigma,meanerr,sigmaerr)
        for ext in ['pdf','png']:
            c.SaveAs('respcomp_{dist}_fit.{ext}'.format(dist=dist,ext=ext))
    return results

def graphVsR(typeInput='histo',histoFile="response_histos.root"):


    resp_vs_z_uncorr = ROOT.TGraphErrors(len(runs))
    reso_vs_z_uncorr = ROOT.TGraph(len(runs))
    resp_vs_z_regr = ROOT.TGraphErrors(len(runs))
    reso_vs_z_regr = ROOT.TGraph(len(runs))

    reso_vs_z_uncorr_fullrms = ROOT.TGraph(len(runs))
    reso_vs_z_regr_fullrms = ROOT.TGraph(len(runs))

    resp_vs_z_uncorr.GetXaxis().SetTitle("z (cm)")
    resp_vs_z_uncorr.GetYaxis().SetTitle("E/E^{peak}")
    resp_vs_z_uncorr.SetTitle("")

    reso_vs_z_uncorr.GetXaxis().SetTitle("z (cm)")
    reso_vs_z_uncorr.GetYaxis().SetTitle("E/E^{peak} resolution")
    reso_vs_z_uncorr.SetTitle("")

    if typeInput=="histo":
        tf = ROOT.TFile(histoFile,"read")
        for i,v in enumerate(runs.items()):
            run,dist = v
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
        for i,v in enumerate(runs.items()):
            run,dist = v

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
    
    resp_vs_z_uncorr.GetYaxis().SetRangeUser(0.85,1.3)
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
        
    reso_vs_z_uncorr.GetYaxis().SetRangeUser(0.05,0.4)
    legend = doLegend(responses,titles,styles,corner="TL",textSize=0.03,legWidth=0.8,nColumns=legCols)
    for ext in ['pdf','png']:
        c.SaveAs("resolution.%s"%ext)

            
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog input.root params_gbrtrain.txt [opts] ')
    (options, args) = parser.parse_args()

    makeResponseHistos(args[0],args[1])
    graphVsR("fit")
