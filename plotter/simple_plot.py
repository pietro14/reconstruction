import ROOT

ROOT.gStyle.SetOptStat(0)

tf = ROOT.TFile('reco_run815.root')
tree = tf.Get('Events')

histos = []
colors = [ROOT.kRed,ROOT.kBlue,ROOT.kOrange]
#histo = ROOT.TH1F('density','',70,0,2500)
histo = ROOT.TH1F('density','',70,0,20)
for it in xrange(1,4):
    h = histo.Clone('h_iter{it}'.format(it=it))
    h.Sumw2()
    #tree.Draw('track_integral/track_nhits>>h_iter{it}'.format(it=it),'track_iteration=={it}'.format(it=it))
    #tree.Draw('track_integral/track_length>>h_iter{it}'.format(it=it),'track_iteration=={it}'.format(it=it))
    tree.Draw('track_length>>h_iter{it}'.format(it=it),'track_iteration=={it}'.format(it=it))
    h.Scale(1./h.Integral())
    h.SetFillColor(colors[it-1])
    h.SetLineColor(colors[it-1])
    h.SetFillStyle(3005)
    histos.append(h)


# legend
(x1,y1,x2,y2) = (0.7, .70, .9, .87)
leg = ROOT.TLegend(x1,y1,x2,y2)
leg.SetFillColor(0)
leg.SetFillColorAlpha(0,0.6)
leg.SetShadowColor(0)
leg.SetLineColor(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)

c = ROOT.TCanvas('c','',600,600)
for ih,h in enumerate(histos):
    h.Draw('hist' if ih==0 else 'hist same')
    h.GetYaxis().SetRangeUser(0,0.25)
    h.GetXaxis().SetTitle('length (mm)')
    leg.AddEntry(h,'iteration {it}'.format(it=ih+1),'f')

leg.Draw()

c.SaveAs('density.pdf')


    
