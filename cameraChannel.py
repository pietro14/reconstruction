#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
import uproot
import numpy as np
import math
import debug_code.tools_lib as tl
import sys

class cameraGeometry:
    def __init__(self,params):
        self.pixelwidth = params['pixelwidth']
        self.cameratype = params['cameratype']
        self.npixx = 0
        self.npixy = 0
        if self.cameratype == 'Flash':
            self.npixx = 2048
            self.npixy = 2048
        if self.cameratype == 'Fusion':
            self.npixx = 2304
            self.npixy = 2304
        if self.cameratype == 'Quest':
            self.npixx = 4096
            self.npixy = 2304
        if self.npixx==0:
            print('\nIt seems you are trying to use a camera which is not supported. Correct this in modules_config/geometry_<detector>.txt.\n ANALYSIS FAILED')
            sys.exit()
        self.name = params['name']
        self.vignette = params['vignette']
        self.xmin = params['xmin'] if 'xmin' in params else 0
        self.xmax = params['xmax'] if 'xmax' in params else self.npixx
        self.ymin = params['ymin'] if 'ymin' in params else 0
        self.ymax = params['ymax'] if 'ymax' in params else self.npixy
        
        
class cameraTools:
    def __init__(self,geometry):
        self.geometry = geometry
        # attach to a dict to make it persistent
        # the matrix is the max size possible, still ok if rebinned (because it is redone from the TH2D when it is readout)
        self.vignetteMap = { self.geometry.name : np.ones((int(self.geometry.npixy),int(self.geometry.npixx))) }

    def pedsub(self,img,pedarr):
        return img - pedarr

    def satur_corr(self,img):
        e = 1.60217662e-7 # electric charge in C
        d2 = 0.015625 # mm^2
        omega = 0.00018 #solid  angle  covered  by  thephotocamera
        alpha = 0.08 # coefficient for varied gas mixture
        sigma0 = 2.5 # small variable of density of charge, we set it as limit for different fit parameters, unit of pC/mm^2
        a0 = 0.1855 # 0.2368 
        #fit function is y=ax^2+bx+c 
        a = a0*e/(d2*alpha*omega) # from fit of charge densities functions, converted to photon units
        b = (1.- 2*a0*sigma0)
        c= a0*sigma0*sigma0*(d2*alpha*omega)/e
        img0 = sigma0*(d2*alpha*omega)/e  # minimal desity of photons, like sigma0 but in photon unit 
        img_sat_corr = np.where(img > img0, (a*img*img + b*img + c), img)
        return img_sat_corr
    
    def zsfullres(self,img_sub,noisearr,nsigma=1):
        img_zs = np.where(img_sub > nsigma * noisearr, img_sub, 0)
        return img_zs
        
    def arrrebin(self,img,rebin):
        newshapex = int(self.geometry.npixx/rebin)
        newshapey = int(self.geometry.npixy/rebin)
        img_rebin = tl.rebin(img,(newshapey,newshapex))
        return img_rebin

    def acceptance(self,img,rowmin,rowmax,colmin,colmax):
        img[:rowmin,:]=0
        img[rowmax:,:]=0
        img[:,:colmin]=0
        img[:,colmax:]=0
        return img
        
    # this returns x,y,z before the zero suppression
    def getImage(self,th2):
        return np.array(th2)

    def noisearray(self,th2):
        noisearr = np.zeros( (th2.GetNbinsX(),th2.GetNbinsY()) )
        for ix in range(th2.GetNbinsX()):
            for iy in range(th2.GetNbinsY()):
                noisearr[ix][iy] = th2.GetBinError(ix+1,iy+1)
        return noisearr

    def loadVignettingMap(self):
        print ("Loading vignette map from: {vf}...".format(vf=self.geometry.vignette))
        det = self.geometry.name
        if self.geometry.cameratype == 'Fusion' or self.geometry.cameratype == 'Quest':
                maptocamera = open('data/Vign_cam_matcher.txt')
       	        maptocameraParams = eval(maptocamera.read())

                if self.geometry.cameratype != maptocameraParams[(self.geometry.vignette).split('/')[1]]: #check that the vignetting map is made with the camera you set for the analysis
                       print ('ERROR! The camera and the vignetting map do not match. Check what you wrote in the configFile or in modules_config/geometry_xxx.txt.\nYou can aso check the vignetting readme in data folder.\nAnalysis FAILED')
                       sys.exit()
                tf = uproot.open(self.geometry.vignette)
                namehmap = 'normmap'
                
                vignetteMapRebinned = tf[namehmap].values().T
                print("py ndim = ",vignetteMapRebinned.ndim)
                rebiny = int(self.geometry.npixy/vignetteMapRebinned.shape[0])
                rebinx = int(self.geometry.npixx/vignetteMapRebinned.shape[1])
                macroPixel = np.zeros((rebiny,rebinx))
                print ("Macro-pixel of the vignetting map has size = ",macroPixel.shape)
                
                for iby in range(vignetteMapRebinned.shape[0]):
                    for ibx in range(vignetteMapRebinned.shape[1]):
                        macroPixel[:,:] = 1./vignetteMapRebinned[iby,ibx]
                        (self.vignetteMap[det])[int(iby*rebiny):int((iby+1)*rebiny),int(ibx*rebinx):int((ibx+1)*rebinx)] = macroPixel
                return self.vignetteMap[det]
                
        else:
            print ('\nWARNING! There is no ',self.geometry.cameratype,' vignetting map. Flat mask of 1 will be used\n')
            return self.vignetteMap[det]

    def vignette_corr(self,img,vignette):
        return np.multiply(img,vignette)