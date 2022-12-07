import os, math, optparse, ROOT
ROOT.gROOT.SetBatch(1)
from simple_plot import getCanvas, doLegend
from mcPlots import doTinyCmsPrelim

ROOT.gROOT.LoadMacro('./fitter/libCpp/RooDoubleCBFast.cc+')

def fitOne(filein,selection,erange,suffix,energies,nbins):

    var = 'sc_integral*6/8000'

    xmin,xmax=erange
    mark = ROOT.kFullSquare
    line = ROOT.kAzure-6

    rf = ROOT.TFile.Open(filein)
    t = rf.Get('Events')
    spectrum = ROOT.TH1F('spectrum','',nbins,xmin,xmax)
    spectrum.SetMarkerStyle(mark)
    spectrum.SetLineColor(ROOT.kBlack)
    t.Draw('{var}>>spectrum'.format(var=var),'({var}>{xmin}) * ({var}<{xmax}) * ({sel})'.format(var=var,xmin=xmin,xmax=xmax,sel=selection))
    deltae = energies[1]-energies[0]
    av=energies[0]

    work = ROOT.RooWorkspace()
    work.factory("mean1[{av},{mmin},{mmax}]".format(av=av,mmin=0.7*av,mmax=1.3*av))
    #work.factory('Gaussian::cb1(x[{xmin},{xmax}],mean1,sigma[1,0.5,1])'.format(xmin=xmin,xmax=xmax))
    if suffix in ['Cu','Rb','Ag']:
        work.factory('CBShape::cb1(x[{xmin},{xmax}],mean1,sigma[1,0.5,3],alpha[1,0.1,10],n[9,1,20])'.format(xmin=xmin,xmax=xmax))
    else:
        work.factory('DoubleCBFast::cb1(x[{xmin},{xmax}],mean1,sigma[1,0.5,1],alpha1[1,0.1,10],n1[5,1,10],alpha2[1,0.05,10],n2[5,1,10])'.format(xmin=xmin,xmax=xmax))
    work.factory("expr::mean2('mean1+delta',mean1,delta[{delta}])".format(delta=deltae))
    #work.factory('Gaussian::cb2(x[{xmin},{xmax}],mean2,sigma[1,0.5,1])'.format(xmin=xmin,xmax=xmax))
    if suffix in ['Cu','Rb','Ag']:
        work.factory('CBShape::cb2(x[{xmin},{xmax}],mean2,sigma,alpha,n)'.format(xmin=xmin,xmax=xmax))
    else:
        work.factory('DoubleCBFast::cb2(x[{xmin},{xmax}],mean2,sigma[1,0.5,1],alpha1,n1,alpha2,n2)'.format(xmin=xmin,xmax=xmax))

    work.factory("a[0,-100,1000]")
    work.factory("b[0,-100,100]")
    work.factory("c[0,-100,100]")
    work.factory("d[0,-100,100]")
    work.factory("e[0,-100,100]")
    #work.factory('Chebychev::bkg(x,{a,b,c,d,e})')
    work.factory('Bernstein::bkg(x,{a,b,c,d})')
    
    x = work.var('x')

    print ("ssss = ",spectrum.Integral())
    work.factory('nSig1[{ns}]'.format(ns=0.2*spectrum.Integral()))
    work.factory("expr::nSig2('nSig1*frac',nSig1,frac[0.2,0,0.2])")
    work.factory('nBkg[{ns}]'.format(ns=0.5*spectrum.Integral()))
    work.factory("SUM::pdfTot(nSig1*cb1,nSig2*cb2,nBkg*bkg)");

    work.Print()
    
    rooData = ROOT.RooDataHist("spectrum","spectrum",ROOT.RooArgList(work.var("x")),spectrum)
    getattr(work,'import')(rooData)

    frame = x.frame()
    frame.SetTitle('')
    
    # fit spectrum
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(mark),ROOT.RooFit.MarkerSize(2),ROOT.RooFit.Name('data'))
    pdf = work.pdf('pdfTot')
    pdf.fitTo(rooData,ROOT.RooFit.Save(),ROOT.RooFit.PrintLevel(10))
    pdf.plotOn(frame,ROOT.RooFit.LineColor(line), ROOT.RooFit.Name('ptot'))
    pdf.plotOn(frame,ROOT.RooFit.Components("bkg"),ROOT.RooFit.LineStyle(ROOT.kDotted), ROOT.RooFit.LineColor(ROOT.kCyan), ROOT.RooFit.Name('pbkg'))
    pdf.plotOn(frame,ROOT.RooFit.Components("cb1"), ROOT.RooFit.LineStyle(ROOT.kDotted), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name('pcb1'))
    pdf.plotOn(frame,ROOT.RooFit.Components("cb2"), ROOT.RooFit.LineStyle(ROOT.kDotted), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name('pcb2'))
    rooData.plotOn(frame,ROOT.RooFit.MarkerStyle(mark))

    frame.GetYaxis().SetTitle("superclusters")
    frame.GetXaxis().SetTitle("Energy (keV)")
    frame.GetXaxis().SetTitleOffset(1.2)

    c = getCanvas()
    frame.Draw()

    # Add legend
    leg = ROOT.TLegend(0.6,0.6,0.85,0.88)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetTextSize(0.02)
    leg.AddEntry("pdata","data with {source}".format(source=suffix),"ep")
    leg.AddEntry("tot","Parametric Total Model","L")
    leg.AddEntry("pcb1","Signal Model 1","L")
    leg.AddEntry("pcb2","Signal Model 2","L")
    leg.AddEntry("pbkg","Background Model","L")
    leg.Draw("Same")

    for ext in ['pdf','png','root']:
        c.SaveAs('energy_{suff}.{ext}'.format(ext=ext,suff=suffix))

    m    = work.var('mean1').getVal()
    merr = work.var('mean1').getError()
    s    = work.var('sigma').getVal()
    serr = work.var('sigma').getError()

    return [m,merr,s,serr]

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    (options, args) = parser.parse_args()

    base_sel = "(TMath::Hypot(sc_xmean-2304/2,sc_ymean-2304/2)<800)"

    energies = {"Fe" : (5.9,6),
                "Cu" : (8.04,  8.91),
                "Rb" : (13.37, 14.97),
                "Mo" : (17.44, 19.63),
                "Ag" : (22.10, 24.99),
                "Ba" : (32.06, 36.55)
                }
        
    
    extras = {"Cu" : "(0.152*sc_length < 80)  * (sc_rms>7.5) * (run>=5801 && run<=5805)",
              "Rb" : "(0.152*sc_length < 80)  * (sc_rms>6.5) * (run>=5806 && run<=5810)",
              "Mo" : "(0.152*sc_length < 80)  * (sc_rms>6.5) * (run>=5811 && run<=5820)",
              "Ag" : "(0.152*sc_length < 80)  * (sc_rms>6.5) * (run>=5821 && run<=5830)",
              "Ba" : "(0.152*sc_length < 120) * (sc_rms>6.5) * (run>=5831 && run<=5845)"}


    eranges = {"Cu" : (2,16),
               "Rb" : (1,25),
               "Mo" : (8,32),
               "Ag" : (9,40),
               "Ba" : (10,60)}

    nbins = {"Cu" : 45,
             "Rb" : 45,
             "Mo" : 60,
             "Ag" : 60,
             "Ba" : 60}
        
    respdic = {}
    respdic["Fe"] = [energies['Fe'],0,0.9,0]

    mat = 'Ba'
    fitOne("trees-lnf/reco_multisource.root ",'({base})*({extra})'.format(base=base_sel,extra=extras[mat]),eranges[mat],mat,energies[mat],nbins[mat])
"""    
    for mat,extracut in extras.items():
        respdic[mat] = fitOne("/Users/emanuele/data/cygnus/RECO/lnf/reco_multisource.root ",'({base})*({extra})'.format(base=base_sel,extra=extras[mat]),eranges[mat],mat)


    fout = ROOT.TFile.Open("linearity.root","recreate")
    resp = ROOT.TGraphErrors(len(energies))
    reso = ROOT.TGraphErrors(len(energies))
    i=0
    for mat,etrue in energies.items():
        resp.SetPoint(i,energies[mat],respdic[mat][0])
        resp.SetPointError(i,0,respdic[mat][1])
        reso.SetPoint(i,energies[mat],100*respdic[mat][2]/respdic[mat][0])
        reso.SetPointError(i,0,100*respdic[mat][3]/respdic[mat][0])
        i+=1

    
    ## response linearity
    c = getCanvas()
    resp.SetTitle('')
    resp.Draw("AP")
    resp.SetMarkerStyle(ROOT.kFullCircle)
    resp.SetMarkerSize(2)
    resp.GetXaxis().SetLimits(0,50)
    resp.GetXaxis().SetTitle("Expected energy (keV)")
    resp.GetYaxis().SetLimits(0,8e4)
    resp.GetYaxis().SetRangeUser(0,50)
    resp.GetYaxis().SetTitle("Measured energy (keV)")

    resp.Fit("pol1")

    line = ROOT.TLine(0,0,50,50)
    line.SetLineColor(ROOT.kBlue)
    line.SetLineStyle(ROOT.kDashed)
    line.Draw()
    
    fout.cd()
    resp.Write()
    c.Write()
    for ext in ['pdf','png','root']:
        c.SaveAs('linearity.{ext}'.format(ext=ext))
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
    reso.GetYaxis().SetTitle("#sigma_{E}/E (%)")        
    for ext in ['pdf','png','root']:
        c.SaveAs('resolution.{ext}'.format(ext=ext))
"""
