import os, math, optparse, ROOT
from array import array
from simple_plot import getCanvas, doLegend

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
    

def makeEff(f1,histo1,f2,histo2,plotdir):
    # numerator: selected events
    tf1 = ROOT.TFile.Open(f1)
    hpass = tf1.Get(histo1)
    hpass.GetYaxis().SetTitle('efficiency')
    hpass.GetXaxis().SetRangeUser(0,140)
    
    # denominator: all events
    tf2 = ROOT.TFile.Open(f2)
    htotal = tf2.Get(histo2)
    htotal.GetYaxis().SetTitle('efficiency')
    htotal.GetXaxis().SetRangeUser(0,140)

    ## default is 68% CL
    teffi68 = ROOT.TEfficiency(hpass,htotal)
    teffi68.SetStatisticOption(ROOT.TEfficiency.kFCP);
    teffi68.SetFillStyle(3004);
    teffi68.SetFillColor(ROOT.kRed);
    teffi68.SetMarkerStyle(ROOT.kFullSquare);
    teffi68.SetMarkerSize(2)
    teffi68.SetLineWidth(2)
    teffi68.SetMarkerColor(ROOT.kBlack);

    ## copy current TEfficiency object and set new confidence level
    teffi90 = ROOT.TEfficiency(teffi68);
    teffi90.SetStatisticOption(ROOT.TEfficiency.kFCP);
    teffi90.SetConfidenceLevel(0.90);
    teffi90.SetFillStyle(3005);
    teffi90.SetFillColor(ROOT.kBlue);
    
    c = getCanvas()
    teffi68.Draw("A4")
    teffi68.Draw("same pe1")
    teffi90.Draw("same4")

    
    ## add legend
    histos = [teffi68,teffi90]
    labels = ['95%','68%']
    styles = ['F','F']
    legend = doLegend(histos,labels,styles,corner='BL')
    legend.Draw('same')
    
    for ext in ['png','pdf']:
        c.SaveAs('{odir}/{var}_effi.{ext}'.format(odir=plotdir,var=histo1,ext=ext))
    

### this is meant to be run on top of ROOT files produced by simple_plots, not on the trees
if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('', '--make'   , type='string'       , default='fitfe' , help='run ambe_miscellanea.py (options = fitfe, efficiency)')
    parser.add_option('', '--outdir' , type='string'       , default='./'      , help='output directory with directory structure and plots')
    (options, args) = parser.parse_args()
                 
    if  options.make == 'fitfe':
        fitFe('plots/ambe/clusters_3sourcesNloCalNeutronsFex1_2020_05_05/energy.root')

    if options.make == 'efficiency':
        var = 'energyExt'
        makeEff('plots/ambe/clusters_3sourcesNloCalNeutronsDensityGt11_2020_05_06/'+var+'.root',var,
                'plots/ambe/clusters_3sourcesNloCalNeutrons_2020_05_05/'+var+'.root',           var,
                options.outdir)
                
                
