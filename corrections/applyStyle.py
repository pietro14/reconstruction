import ROOT
ROOT.gROOT.SetBatch(True)

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
        doSpam(textRight,(0.5 if doWide else .55)+xoffs, .94, .82+xoffs if doWide else .93+xoffs, .94, align=32, textSize=textSize)

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


def doResponse():
    tf = ROOT.TFile.Open("response_graphs.root")

    c = getCanvas('c')
    
    resp_vs_z_uncorr = tf.Get("resp_vs_z_uncorr")
    resp_vs_z_regr = tf.Get("resp_vs_z_regr")
    
    resp_vs_z_uncorr.GetXaxis().SetTitleFont(42)
    resp_vs_z_uncorr.GetXaxis().SetTitleSize(0.05)
    resp_vs_z_uncorr.GetXaxis().SetTitleOffset(1.1)
    resp_vs_z_uncorr.GetXaxis().SetLabelFont(42)
    resp_vs_z_uncorr.GetXaxis().SetLabelSize(0.045)
    resp_vs_z_uncorr.GetXaxis().SetLabelOffset(0.007)
    resp_vs_z_uncorr.GetYaxis().SetTitleFont(42)
    resp_vs_z_uncorr.GetYaxis().SetTitleSize(0.05)
    resp_vs_z_uncorr.GetYaxis().SetTitleOffset(1.3)
    resp_vs_z_uncorr.GetYaxis().SetLabelFont(42)
    resp_vs_z_uncorr.GetYaxis().SetLabelSize(0.045)
    resp_vs_z_uncorr.GetYaxis().SetLabelOffset(0.007)
    resp_vs_z_uncorr.GetYaxis().SetTitle("E/E^{mpv}")
    resp_vs_z_uncorr.GetXaxis().SetTitle("z (cm)")
    resp_vs_z_uncorr.GetXaxis().SetNdivisions(510)


    resp_vs_z_uncorr.Draw("Ape")
    resp_vs_z_regr.Draw("pe")
    
    responses = [resp_vs_z_uncorr,resp_vs_z_regr]
    titles = ['E_{rec}','E']
    styles = ['pl','pl']

    resp_vs_z_uncorr.GetYaxis().SetRangeUser(0.5,1.3)
    legend = doLegend(responses,titles,styles,corner="TL")

    doTinyCmsPrelim(hasExpo = False,textSize=0.04, options=None,doWide=False)

    for ext in ['pdf','png','root']:
        c.SaveAs("response.%s"%ext)




def doResolution():
    tf = ROOT.TFile.Open("response_graphs.root")

    c = getCanvas('c')


    reso_vs_z_uncorr = tf.Get("reso_vs_z_uncorr")
    reso_vs_z_regr = tf.Get("reso_vs_z_regr")
    reso_vs_z_uncorr_fullrms = tf.Get("reso_vs_z_uncorr_fullrms")
    reso_vs_z_regr_fullrms = tf.Get("reso_vs_z_regr_fullrms")
    
    reso_vs_z_uncorr.GetXaxis().SetTitleFont(42)
    reso_vs_z_uncorr.GetXaxis().SetTitleSize(0.05)
    reso_vs_z_uncorr.GetXaxis().SetTitleOffset(1.1)
    reso_vs_z_uncorr.GetXaxis().SetLabelFont(42)
    reso_vs_z_uncorr.GetXaxis().SetLabelSize(0.045)
    reso_vs_z_uncorr.GetXaxis().SetLabelOffset(0.007)
    reso_vs_z_uncorr.GetYaxis().SetTitleFont(42)
    reso_vs_z_uncorr.GetYaxis().SetTitleSize(0.05)
    reso_vs_z_uncorr.GetYaxis().SetTitleOffset(1.3)
    reso_vs_z_uncorr.GetYaxis().SetLabelFont(42)
    reso_vs_z_uncorr.GetYaxis().SetLabelSize(0.045)
    reso_vs_z_uncorr.GetYaxis().SetLabelOffset(0.007)
    reso_vs_z_uncorr.GetYaxis().SetTitle("#sigma_{E} / E^{mpv}")
    reso_vs_z_uncorr.GetXaxis().SetTitle("z (cm)")
    reso_vs_z_uncorr.GetXaxis().SetNdivisions(510)

    
    responses = [reso_vs_z_uncorr,reso_vs_z_regr,reso_vs_z_uncorr_fullrms,reso_vs_z_regr_fullrms]
    titles = ['E_{rec} (#sigma_{G})','E (#sigma_{G})','E_{rec} (rms)','E (rms)']
    styles = ['pl','pl','pl','pl']
    
    reso_vs_z_uncorr.Draw("Ape")
    reso_vs_z_regr.Draw("pe")
    reso_vs_z_uncorr_fullrms.Draw("pe")
    reso_vs_z_regr_fullrms.Draw("pe")
    legCols=2

    reso_vs_z_uncorr.GetYaxis().SetRangeUser(0.0,0.4)
    reso_vs_z_uncorr.GetXaxis().SetRangeUser(0.0,50)
    legend = doLegend(responses,titles,styles,corner="TL",textSize=0.03,legWidth=0.8,nColumns=legCols)

    doTinyCmsPrelim(hasExpo = False,textSize=0.04, options=None,doWide=False)

    for ext in ['pdf','png','root']:
        c.SaveAs("resolution.%s"%ext)


def doHistory(tfile='pedmean.root'):
    tf = ROOT.TFile.Open(tfile)
    history = tf.Get("history")

    c = getCanvas('c')
    
    history.GetXaxis().SetTitleFont(42)
    history.GetXaxis().SetTitleSize(0.05)
    history.GetXaxis().SetTitleOffset(1.1)
    history.GetXaxis().SetLabelFont(42)
    history.GetXaxis().SetLabelSize(0.045)
    history.GetXaxis().SetLabelOffset(0.007)

    history.GetXaxis().SetTimeDisplay(1)
    da = ROOT.TDatime(1970,6,9,12,00,00);
    history.GetXaxis().SetTimeOffset(da.Convert())
    history.GetXaxis().SetTimeFormat("#splitline{%Y}{%d\/%m}")
    
    history.GetYaxis().SetTitleFont(42)
    history.GetYaxis().SetTitleSize(0.05)
    history.GetYaxis().SetTitleOffset(1.3)
    history.GetYaxis().SetLabelFont(42)
    history.GetYaxis().SetLabelSize(0.045)
    history.GetYaxis().SetLabelOffset(0.007)
    #history.GetYaxis().SetTitle("E/E^{mpv}")
    #history.GetXaxis().SetTitle("z (cm)")
    #history.GetXaxis().SetNdivisions(510)

    history.Draw("ape")

    doTinyCmsPrelim(hasExpo = False,textSize=0.04, options=None,doWide=False)

    for ext in ['pdf','png','root']:
        c.SaveAs("%s.%s"%(tfile.split('.')[0]+"_history",ext))
    
    
if __name__ == '__main__':
    doHistory()

    






