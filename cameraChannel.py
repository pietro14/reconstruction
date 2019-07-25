#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
import numpy as np
from root_numpy import hist2array
import math
import debug_code.tools_lib as tl

class cameraGeometry:
    def __init__(self):
        # 240 mm = major axis of the ellipse = 2048 pixels
        # ORANGE
        # 1 px - 55*10-6 m
        # 1 ph/ev - 440 V

        # LEMON
        # 1 px - 125*10-6 m
        # 0,12 ph/ev - 460V
        # 0,10 ph/ev - 450V
        # 0,06 ph/ev - 440V
        
        self.pixelwidth = 125E-3 # mm
        
class cameraTools:
    def __init__(self):
        pass
        
    def isGoodChannelSafe(self,pedmap,ix,iy):
        pedvals = []; pedrmss = []
        for x in range(ix-1,ix+2):
            for y in range(iy-1,iy+2):
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
    
    def zs(self,th2,pedmap,nsigma=3.2,plot=False):
        cimax = 200       # valori del cut sull'imagine
        #print "zero suppressing..."
        nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
        xmin,xmax=(th2.GetXaxis().GetXmin(),th2.GetXaxis().GetXmax())
        ymin,ymax=(th2.GetYaxis().GetXmin(),th2.GetYaxis().GetXmax())
        th2_zs = ROOT.TH2D(th2.GetName()+'_zs',th2.GetName()+'_zs',nx,xmin,xmax,ny,ymin,ymax)
        th2_unzs = th2_zs.Clone(th2.GetName()+'_unzs')
        th2_zs.SetDirectory(None); th2_unzs.SetDirectory(None);
        for ix in range(1,nx+1):
            for iy in range(1,ny+1):
                ped   = pedmap.GetBinContent(ix,iy)
                noise = pedmap.GetBinError(ix,iy)
                if not self.isGoodChannelFast(ped,noise): continue
                z = th2.GetBinContent(ix,iy)-ped
                z1 = th2.GetBinContent(ix,iy)
                th2_unzs.SetBinContent(ix,iy,z)
                #print("x,y,z=",ix," ",iy," ",z,"   noise = ",noise , "ped-ib = ", ped_ixb," - ", ped_ixb)
                if (z>nsigma*noise) & (z1<cimax):
                    th2_zs.SetBinContent(ix,iy,z)
                    #print "x,y,z=",ix," ",iy," ",z,"   noise = ",noise
        #th2_zs.GetZaxis().SetRangeUser(0,1)
        if plot:
            canv = ROOT.TCanvas('zs','',600,600)
            th2_zs.Draw('colz')
            for ext in ['pdf','root']:
                canv.SaveAs('{name}.{ext}'.format(name=th2.GetName()+'_zs',ext=ext))
        return (th2_zs,th2_unzs)

    def zsfullres(self,img,pedarr,noisearr,nsigma=1):
        img_sub   = img - (pedarr + nsigma * noisearr)
        img_zs = np.where(img_sub > 0, img_sub, 0)
        return img_zs
        
    def arrrebin(self,img,rebin):
        newshape = 2048/rebin
        img_rebin = tl.rebin(img,(newshape,newshape))
        return img_rebin
        
    def zshits(self,hits,pedmap,nsigma=3,plot=False):
        print("zero suppressing a subset of ",len(hits)," hits...")
        ret=[]
        for h in hits:
            x,y,z=h[0],h[1],h[2]
            # warning: this works only for FR coordinates, otherwise need the findbin etc...
            ped = pedmap.GetBinContent(x+1,y+1)
            noise = pedmap.GetBinError(x+1,y+1)
            z_ps = z-ped
            if self.isGoodChannelFast(ped,noise) and z_ps>nsigma*noise and z_ps<30:
                ret.append(x,y,z_ps)
            else:
                ret.append(x,y,0)
        return ret
    
    def getData(self,th2):
        if not th2.InheritsFrom("TH2"):
            print("ERROR! The input object should be a TH2")
            return []
        x_bins = th2.GetNbinsX()
        y_bins = th2.GetNbinsY()
        bins = np.zeros((x_bins,y_bins))
        for y_bin in range(y_bins): 
            for x_bin in range(x_bins): 
                z = th2.GetBinContent(x_bin + 1,y_bin + 1)
                bins[x_bin,y_bin] = z
        return bins

    # this returns only x,y,z after zero suppression
    def getActiveCoords(self,th2):
        ret = []
        if not th2.InheritsFrom("TH2"):
            print("ERROR! The input object should be a TH2")
            return ret
        x_bins = th2.GetNbinsX()
        y_bins = th2.GetNbinsY()
        for y_bin in range(y_bins): 
            for x_bin in range(x_bins): 
                x = th2.GetXaxis().GetBinCenter(x_bin+1)
                y = th2.GetYaxis().GetBinCenter(y_bin+1)
                z = th2.GetBinContent(x_bin + 1,y_bin + 1)
                if z>0:
                    #ret.append((x,y,math.log(z)))
                    ret.append((x,y,z))
        return np.array(ret)

    # this returns x,y,z before the zero suppression
    def getImage(self,th2):
        return hist2array(th2)

    def noisearray(self,th2):
        noisearr = np.zeros( (th2.GetNbinsX(),th2.GetNbinsY()) )
        for ix in range(th2.GetNbinsX()):
            for iy in range(th2.GetNbinsY()):
                noisearr[ix][iy] = th2.GetBinError(ix+1,iy+1)
        return noisearr

    def getRestrictedImage(self,th2,xmin,xmax,ymin,ymax):
        nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
        nxp = xmax-xmin
        nyp = ymax-ymin
        th2_rs = ROOT.TH2D(th2.GetName()+'_rs',th2.GetName()+'_rs',nxp,xmin,xmax,nyp,ymin,ymax)
        th2_rs.SetDirectory(None)
        for ix,x in enumerate(range(xmin,xmax)):
            for iy,y in enumerate(range(ymin,ymax)):
                orig_ixb = th2.GetXaxis().FindBin(x)
                orig_iyb = th2.GetYaxis().FindBin(y)
                z = th2.GetBinContent(orig_ixb,orig_iyb)
                th2_rs.SetBinContent(ix+1,iy+1,z)
        # canv = ROOT.TCanvas('zs','',600,600)
        # th2_rs.Draw('colz')
        # for ext in ['png','pdf','root']:
        #     canv.SaveAs('{name}_rs.{ext}'.format(name=th2.GetName(),ext=ext))
        return th2_rs

