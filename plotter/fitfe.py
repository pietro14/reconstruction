import os, math, optparse, ROOT
from array import array
from simple_plot import getCanvas

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

def fitFe(rfile):
    tf = ROOT.TFile.Open(rfile)
    histo_sig = tf.Get("energyfe_diff")
    histo_sig.SetMarkerColor(ROOT.kBlack)
    histo_sig.SetLineColor(ROOT.kBlack)
    histo_sig.GetXaxis().SetTitle('energy (keV)')
    histo_sig.GetXaxis().SetTitleSize(0.05)
    histo_sig.GetYaxis().SetTitle('superclusters (bkg subtracted)')
    
    c = getCanvas()
    histo_sig.Draw("pe 1")

    par = array( 'd', 6*[0.] )

    g1 = ROOT.TF1("g1","gaus",3,9);
    g2 = ROOT.TF1("g2","gaus",11,16);
    total = ROOT.TF1('total','gaus(0)+gaus(3)',3,16)
    total.SetLineColor(ROOT.kBlue+1)
    histo_sig.Fit('g1','RS')
    histo_sig.Fit('g2','R+S')
    mean  = g1.GetParameter(1); mErr = g1.GetParError(1)
    sigma = g1.GetParameter(2); mSigma = g1.GetParError(2)
    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.03)
    lat.DrawLatex(0.55, 0.70, "m_{{1}} = {m:.2f} #pm {em:.2f} keV".format(m=mean,em=mErr))
    lat.DrawLatex(0.55, 0.65, "#sigma_{{1}} = {s:.2f} #pm {es:.2f} keV".format(s=sigma,es=mSigma))

    mean  = g2.GetParameter(1); mErr = g2.GetParError(1)
    sigma = g2.GetParameter(2); mSigma = g2.GetParError(2)
    lat.DrawLatex(0.55, 0.50, "m_{{2}} = {m:.2f} #pm {em:.2f} keV".format(m=mean,em=mErr))
    lat.DrawLatex(0.55, 0.45, "#sigma_{{2}} = {s:.2f} #pm {es:.2f} keV".format(s=sigma,es=mSigma))



    par1 = g1.GetParameters()
    par2 = g2.GetParameters()

    par[0], par[1], par[2] = par1[0], par1[1], par1[2]
    par[3], par[4], par[5] = par2[0], par2[1], par2[2]
    total.SetParameters( par )
    histo_sig.Fit( total, 'R+' )
    
    c.SaveAs('fe_diff_simplefit.pdf')
    

if __name__ == "__main__":

    fitFe('plots/ambe/clusters_3sourcesNloCalNeutronsFex1_2020_05_05/energy.root')
