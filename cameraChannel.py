#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)

class cameraGeometry:
    def __init__(self):
        self.pixelwidth = 100E-3 # mm
        
class cameraTools:
    def __init__(self):
        pass
        
    def isGoodChannelSafe(self,pedmap,ix,iy):
        pedvals = []; pedrmss = []
        for x in xrange(ix-1,ix+2):
            for y in xrange(iy-1,iy+2):
                pedvals.append(pedmap.GetBinContent(x,y))
                pedrmss.append(pedmap.GetBinError(x,y))
        if any([p>110 for p in pedvals]): return False
        if any([r<0.2 for r in pedvals]): return False
        #if pedrms > 5: return False
        return True

    def isGoodChannelFast(self,pedval,pedrms):
        if pedval > 110: return False
        if pedrms < 0.2: return False
        if pedrms > 5: return False
        return True
    
    def zs(self,th2,pedmap,nsigma=2,plot=False):
        nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
        xmin,xmax=(th2.GetXaxis().GetXmin(),th2.GetXaxis().GetXmax())
        ymin,ymax=(th2.GetYaxis().GetXmin(),th2.GetYaxis().GetXmax())
        th2_zs = ROOT.TH2D(th2.GetName()+'_zs',th2.GetName()+'_zs',nx,xmin,xmax,ny,ymin,ymax)
        th2_zs.SetDirectory(None)
        for ix in xrange(1,nx+1):
            for iy in xrange(1,ny+1):
                x = th2.GetXaxis().GetBinCenter(ix)
                y = th2.GetYaxis().GetBinCenter(iy)
                ped_ixb = pedmap.GetXaxis().FindBin(x)
                ped_iyb = pedmap.GetYaxis().FindBin(y)
                ped = pedmap.GetBinContent(ped_ixb,ped_iyb)
                noise = pedmap.GetBinError(ped_ixb,ped_iyb)
                if not self.isGoodChannelFast(ped,noise): continue                
                z = max(th2.GetBinContent(ix,iy)-ped,0)
                if z>nsigma*noise:
                    th2_zs.SetBinContent(ix,iy,z)
                    #print "x,y,z=",ix," ",iy," ",z,"   noise = ",noise
        #th2_zs.GetZaxis().SetRangeUser(0,1)
        if plot:
            canv = ROOT.TCanvas('zs','',600,600)
            th2_zs.Draw('colz')
            for ext in ['png','pdf']:
                canv.SaveAs('{name}.{ext}'.format(name=th2.GetName()+'_zs',ext=ext))
        return th2_zs
