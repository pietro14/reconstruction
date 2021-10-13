import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
from root_numpy import tree2array,fill_hist,fill_profile
from gbr_trainer import GBRLikelihoodTrainer, getCanvas, doLegend
import math, joblib

# map run - z distance
#          no source        8.7       8.4        8.1k       5k
#runs = {4433: 50, 4441: 46, 4448: 36, 4455: 26, 4463: 6}
runs = {4441: 46, 4448: 36, 4455: 26}

def response_vsrun(inputfile,params,runN):
    params['selection'] = "({sel}) && run=={run}".format(sel=params['selection'],run=runN)
    print ("New selection SEL = {sel}".format(sel=params['selection']))
    GBR = GBRLikelihoodTrainer(params)
    X,y = GBR.get_dataset(inputfile)
    filename = "gbrLikelihood_q0.50.sav"
    model = joblib.load(filename)
    y_pred = model.predict(X)

    h = ROOT.TH1F('h','',200,0,2)
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

    energy_2D = ROOT.TH3D('energy_2D','',16,0,2304,16,0,2304,40,0,2)
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
        print("ev {iev} has x,y=({x},{y}) and y = {z}".format(iev=i,x=xy[i][0],y=xy[i][1],z=y[i]))
        energy_2Ds['uncorr'].Fill(xy[i][1],xy[i][0],y[i])
        energy_2Ds['regr'].Fill(xy[i][1],xy[i][0],y_pred[i])
        energy_1Ds['uncorr'].Fill(math.hypot(xy[i][1]-center,xy[i][0]-center),y[i])
        energy_1Ds['regr'].Fill(math.hypot(xy[i][1]-center,xy[i][0]-center),y_pred[i])

    energy_2D_mode = ROOT.TH2D('energy_2D_mode','',16,0,2304,16,0,2304)
    energy_1D_mode = ROOT.TH1D('energy_1D_mode','',20,200,2100/math.sqrt(2.))
    energy_2D_mode.GetXaxis().SetTitle("ix")
    energy_2D_mode.GetYaxis().SetTitle("iy")
    energy_2D_mode.GetZaxis().SetTitle("mode")
    energy_2D_mode.GetZaxis().SetRangeUser(0.8,1.3)
    
    energy_2D_modes = {}
    energy_2D_modes['uncorr'] = energy_2D_mode.Clone('energy_2D_mode_uncorr')
    energy_2D_modes['regr'] = energy_2D_mode.Clone('energy_2D_mode_regr')
    
    c = getCanvas('c')
    for k,h3D in energy_2Ds.items():
        for ix in range(1,h3D.GetNbinsX()+1):
            for iy in range(1,h3D.GetNbinsY()+1):
                maxZ = -1; zbinmax = -1
                for iz in range(1,h3D.GetNbinsZ()+1):
                    if h3D.GetBinContent(ix,iy,iz) > maxZ:
                        maxZ = h3D.GetBinContent(ix,iy,iz)
                        zbinmax = iz
                energy_2D_modes[k].SetBinContent(ix,iy,h3D.GetZaxis().GetBinCenter(zbinmax))
        energy_2D_modes[k].Draw("colz")
        c.SaveAs("energy_2D_{name}.png".format(name=k))

    # energy_1Ds['uncorr'].SetMarkerColor(ROOT.kBlack)
    # energy_1Ds['regr'].SetMarkerColor(ROOT.kRed)
    # energy_1Ds['uncorr'].Draw("pe1")
    # energy_1Ds['regr'].Draw("pe1 same")
    # c.SaveAs("energy_1D.png")
        
    
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog input.root params_gbrtrain.txt [opts] ')
    (options, args) = parser.parse_args()

    config = open(args[1], "r")
    params = eval(config.read())
    config.close()
    response2D(args[0],params)

    c = getCanvas('c')

    resp_vs_z_uncorr = ROOT.TGraphErrors(len(runs))
    reso_vs_z_uncorr = ROOT.TGraph(len(runs))
    resp_vs_z_regr = ROOT.TGraphErrors(len(runs))
    reso_vs_z_regr = ROOT.TGraph(len(runs))

    resp_vs_z_uncorr.GetXaxis().SetTitle("z (cm)")
    resp_vs_z_uncorr.GetYaxis().SetTitle("E/E^{peak}")
    resp_vs_z_uncorr.SetTitle("")

    reso_vs_z_uncorr.GetXaxis().SetTitle("z (cm)")
    reso_vs_z_uncorr.GetYaxis().SetTitle("E/E^{peak} resolution")
    reso_vs_z_uncorr.SetTitle("")

    for i,v in enumerate(runs.items()):
        config = open(args[1], "r")
        params = eval(config.read()) # re-read to reload run
        config.close()
        run,dist = v
        print("filling graph with point ",i," and run = ",run," corresponding to distance ",dist)
        histos = response_vsrun(args[0],params,run)

        resp_vs_z_uncorr.SetPoint(i,dist,histos[0].GetXaxis().GetBinCenter(histos[0].GetMaximumBin()))
        resp_vs_z_uncorr.SetPointError(i,0,histos[0].GetXaxis().GetBinWidth(histos[0].GetMaximumBin()))

        resp_vs_z_regr.SetPoint(i,dist,histos[1].GetXaxis().GetBinCenter(histos[1].GetMaximumBin()))
        resp_vs_z_regr.SetPointError(i,0,histos[0].GetXaxis().GetBinWidth(histos[1].GetMaximumBin()))
        
        reso_vs_z_uncorr.SetPoint(i,dist,histos[0].GetRMS())
        reso_vs_z_regr.SetPoint(i,dist,histos[1].GetRMS())

    # response
    resp_vs_z_uncorr.SetMarkerStyle(ROOT.kFullCircle)
    resp_vs_z_regr.SetMarkerStyle(ROOT.kFullCircle)
    resp_vs_z_uncorr.SetMarkerSize(2)
    resp_vs_z_regr.SetMarkerSize(2)
    resp_vs_z_uncorr.SetMarkerColor(ROOT.kRed)
    resp_vs_z_uncorr.SetLineColor(ROOT.kRed)
    resp_vs_z_regr.SetMarkerColor(ROOT.kBlack)
    resp_vs_z_uncorr.Draw("Ape")
    resp_vs_z_uncorr.Fit("pol1")
    resp_vs_z_regr.Draw("pe")
    resp_vs_z_regr.Fit("pol1")

    responses = [resp_vs_z_uncorr,resp_vs_z_regr]
    titles = ['uncorrected','regression']
    styles = ['pl','pl']
    
    f = resp_vs_z_regr.GetListOfFunctions().FindObject("pol1")
    if f!=None:
        f.SetLineColor(ROOT.kBlack)
        f.SetLineWidth(3)
        f.SetLineStyle(2)
    resp_vs_z_uncorr.GetYaxis().SetRangeUser(0.8,1.2)
    legend = doLegend(responses,titles,styles,corner="TL")    
    c.SaveAs("response.png")

    # resolution
    reso_vs_z_uncorr.SetMarkerStyle(ROOT.kFullCircle)
    reso_vs_z_regr.SetMarkerStyle(ROOT.kFullCircle)
    reso_vs_z_uncorr.SetMarkerColor(ROOT.kRed)
    reso_vs_z_regr.SetMarkerColor(ROOT.kBlack)
    reso_vs_z_uncorr.SetMarkerSize(2)
    reso_vs_z_regr.SetMarkerSize(2)
    reso_vs_z_uncorr.Draw("Ape")
    reso_vs_z_regr.Draw("pe")
    reso_vs_z_uncorr.GetYaxis().SetRangeUser(0.10,0.25)
    legend = doLegend(responses,titles,styles,corner="TL")

    c.SaveAs("resolution.png")
