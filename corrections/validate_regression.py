import ROOT
ROOT.gROOT.SetBatch(True)

from test_regression import getFriend
from gbr_trainer import getCanvas,doLegend

def onePlot(source, params):
    tf = ROOT.TFile.Open(params[0])
    tree = tf.Get("Events")
    frf = getFriend(params[0])
    tree.AddFriend("Friends",frf)

    selection = "run>={rmin} && run<={rmax}".format(rmin=params[1],rmax=params[2])
    c = getCanvas('c')
    histo1 = ROOT.TH1F("histo1","",200,0,params[3])
    histo2 = histo1.Clone("histo2")
    tree.Draw("sc_integral >> histo1",selection)
    tree.Draw("sc_regr_integral >> histo2",selection)

    histo1.GetXaxis().SetLabelSize(0.03)
    histo1.GetXaxis().SetLabelFont(42)
    histo1.GetXaxis().SetTitleSize(0.04)
    histo1.GetYaxis().SetLabelSize(0.03)
    histo1.GetYaxis().SetLabelFont(42)
    histo1.GetYaxis().SetTitleSize(0.04)
    histo1.GetYaxis().SetTitleFont(42)
    histo1.GetYaxis().SetTitleOffset(1.5)
    histo1.GetXaxis().SetTitle("integral")
    histo1.GetYaxis().SetTitle("superclusters")

    
    histo1.SetFillColorAlpha(ROOT.kViolet-4,0.5)
    histo2.SetMarkerStyle(ROOT.kFullDotLarge)
    histo2.SetLineColor(ROOT.kBlack)

    histo1.Draw("hist")
    histo2.Draw("pe same")
    ycut = 0.2 if source=='fe55' else 0.05
    if source=="bkg": c.SetLogy() # to make the distribution agreement well visible in the full spectrum
    else: histo1.GetYaxis().SetRangeUser(0,ycut * histo1.GetMaximum()) # truncate the peak of fakes and make the tail visible, with also the peaks (invisible in log scale)

    histos = [histo1,histo2]
    labels = ["raw","mean regression"]
    styles = ['f','pe']

    legend = doLegend(histos,labels,styles,corner="TC")
    legend.Draw()
    
    for ext in ['png','pdf']:
        c.SaveAs('integral_regression_{source}.{ext}'.format(source=source,ext=ext))
    
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog ')
    (options, args) = parser.parse_args()

    #        source   [file, runmin, runmax, max_integral]
    sources = {'fe55' : ('trees/fe55.root',        5867, 5911, 2e4),
               'CuMo' : ('trees/multiSource.root', 5811, 5820, 5e4),
               'CuAg' : ('trees/multiSource.root', 5821, 5830, 5e4),
               'bkg'  : ('trees/noSource.root',    5861, 5866, 2e4) }
    
    for s,params in sources.items():
        onePlot(s,params)
        

    
             
