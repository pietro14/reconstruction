from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText
from ROOT import gROOT
import ROOT

#ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptStat(1111111)

cons = 0.06 # for 440V
#cons = 0.10 # for 450V
#cons = 0.12 # for 460V

lenFac = 117E-3

ifac = [(5.9/800)*(cons),(5.9/800)*(cons)/(10),1/10,1/10,1,lenFac**2]

name   = ["energy","energyPlength","length","width","ratioLW","Linearity"]
var    = ["track_integral","track_integral/track_length","track_length","track_width","track_length/track_width","track_nhits/track_length/track_width"]
x_axis = ["Energy (keV)","Energy (keV/cm)","Length (cm)","Width (cm)","Length/Width","Linearity"]
x_lim  = [ 10,0.25,   2,  1,   7,  1]
y_lim  = [0.3,0.35,0.10,0.1,0.08,0.05]

#vr2 = "track_nhits";
#vr = vr1+":"+vr2

run = 815

filename2 = 'reco_%d_2d.root' % (run)
filename3 = 'reco_%d_3d.root' % (run)

tf2 = ROOT.TFile(filename2)
tree2 = tf2.Get('Events')

tf3 = ROOT.TFile(filename3)
tree3 = tf3.Get('Events')

hist1 = '2D method'
hist2 = '3D method'

colors = [ROOT.kRed,ROOT.kBlue]

for i in range(0,6):
    #i = 5
    vr = var[i]
    fac = ifac[i]

    for it in range(1,4):

        c = ROOT.TCanvas('','',800,600)

        #-------2D Method------#
        h = ROOT.TH1F('{hist1}'.format(hist1=hist1),'Histograms for iteration {it} - Run {run}'.format(it=it,run=run),200,0,x_lim[i])
        h1 = h.Clone('{hist2}'.format(hist2=hist2))
        h.Sumw2()
        tree2.Draw("{fac}*{vr}>>{hist1}".format(fac=fac,vr=vr,hist1=hist1),"track_iteration=={it}".format(it=it))
        h.Scale(1./h.Integral())
        h.SetFillColor(colors[0])
        h.SetLineColor(colors[0])
        h.SetFillStyle(3005)

        #-------3D Method------#

        tree3.Draw("{fac}*{vr}>>{hist2}".format(fac=fac,vr=vr,hist2=hist2),"track_iteration=={it}".format(it=it))
        h1.Sumw2()
        h1.Scale(1./h1.Integral())
        h1.SetFillColor(colors[1])
        h1.SetLineColor(colors[1])
        h1.SetFillStyle(3005)

        #-------Legend------#
        difx = 0.245
        xpos = 0.5
        (x1,y1,x2,y2) = (xpos, .75, xpos+difx, .89)
        leg = ROOT.TLegend(x1,y1,x2,y2)

        h.Draw('hist')
        h1.Draw('hist sames')
        ROOT.gPad.Update()

        h.GetYaxis().SetRangeUser(0,y_lim[i])
        h.GetXaxis().SetRangeUser(0,x_lim[i])
        h.GetXaxis().SetTitle(x_axis[i])
        h.GetYaxis().SetTitle("Normalized histogram")

        #-------Stat Box --------#
        dif = 0.28
        ypos = 0.35
        st = h1.FindObject("stats")
        #st.SetX1NDC(0.2); #new x start position
        #st.SetX2NDC(0.4); #new x end position
        st.SetY1NDC(ypos); #new y start position
        st.SetY2NDC(ypos+dif); #new y end position

        leg.AddEntry(h,'{hist1}'.format(hist1=hist1),'f')
        leg.AddEntry(h1,'{hist2}'.format(hist2=hist2),'f')
        leg.Draw()
        c.SetGrid()

        #c.Draw()
        c.SaveAs('plots/Plot_Comp_Run%d_i%d_%s.pdf' % (run,it,name[i]))
