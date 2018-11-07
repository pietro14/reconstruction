#!/usr/bin/env python

import os,math,sys,random
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.SetBatch(True)
from cameraChannel import cameraTools

from snakes import SnakesProducer
from output import OutputTree
from treeVars import AutoFillTreeProducer

import utilities
utilities = utilities.utils()

class analysis:

    def __init__(self,rfile,options):
        self.xmax = 2048
        self.rebin = options.rebin        
        self.rfile = rfile
        self.options = options
        self.pedfile_name = '{base}_ped_rebin{rb}.root'.format(base=os.path.splitext(self.rfile)[0],rb=self.rebin)
        if not os.path.exists(self.pedfile_name):
            print "WARNING: pedestal file ",self.pedfile_name, " not existing. First calculate them..."
            self.calcPedestal(options)
        self.pedfile_fullres_name = '{base}_ped_rebin1.root'.format(base=os.path.splitext(self.rfile)[0])
        if not os.path.exists(self.pedfile_fullres_name):
            print "WARNING: pedestal file with full resolution ",self.pedfile_fullres_name, " not existing. First calculate them..."
            self.calcPedestal(options,1)
        if not options.justPedestal:
           print "Pulling pedestals..."
           # first the one for clustering with rebin
           pedrf = ROOT.TFile.Open(self.pedfile_name)
           self.pedmap = pedrf.Get('pedmap').Clone()
           self.pedmap.SetDirectory(None)
           self.pedmean = pedrf.Get('pedmean').GetMean()
           self.pedrms = pedrf.Get('pedrms').GetMean()
           pedrf.Close()
           # then the full resolution one
           pedrf_fr = ROOT.TFile.Open(self.pedfile_fullres_name)
           self.pedmap_fr = pedrf_fr.Get('pedmap').Clone()
           self.pedmap_fr.SetDirectory(None)
           pedrf_fr.Close()

    # the following is needed for multithreading
    def __call__(self,evrange=(-1,-1,-1)):
        if evrange[0]==-1:
            outfname = self.options.outFile
        else:
            outfname = '{base}_chunk{ij}.root'.format(base=self.options.outFile.split('.')[0],ij=evrange[0])
        self.beginJob(outfname)
        self.reconstruct(evrange)
        self.endJob()
        
    def beginJob(self,outfname):
        # prepare output file
        self.outputFile = ROOT.TFile.Open(outfname, "RECREATE")
        # prepare output tree
        self.outputTree = ROOT.TTree("Events","Tree containing reconstructed quantities")
        self.outTree = OutputTree(self.outputFile,self.outputTree)
        self.autotree = AutoFillTreeProducer(self.outTree)

        self.outTree.branch("run", "I")
        self.outTree.branch("event", "I")
        if self.options.daq == 'midas': self.autotree.createPMTVariables()
        self.autotree.createCameraVariables()

    def endJob(self):
        self.outTree.write()
        self.outputFile.Close()
        
    def getNEvents(self):
        tf = ROOT.TFile.Open(self.rfile)
        ret = len(tf.GetListOfKeys()) if self.options.daq!='midas' else int(len(tf.GetListOfKeys())/2) 
        tf.Close()
        return ret

    def calcPedestal(self,options,alternativeRebin=-1):
        maxImages=options.numPedEvents
        nx=ny=self.xmax
        rebin = self.rebin if alternativeRebin<0 else alternativeRebin
        nx=int(nx/rebin); ny=int(ny/rebin); 
        pedfile = ROOT.TFile.Open(self.pedfile_name,'recreate')
        pedmap = ROOT.TProfile2D('pedmap','pedmap',nx,0,self.xmax,ny,0,self.xmax,'s')
        tf = ROOT.TFile.Open(self.rfile)
        if options.pedExclRegion:
            print "Excluding the following region from the calculation of pedestals:"
            xmin,xmax = options.pedExclRegion.split(',')[0].split(':')
            ymin,ymax = options.pedExclRegion.split(',')[1].split(':')
            print "xrange excl = ({xmin},{xmax})".format(xmin=xmin,xmax=xmax)
            print "yrange excl = ({ymin},{ymax})".format(ymin=ymin,ymax=ymax)
            
        for i,e in enumerate(tf.GetListOfKeys()):
            if maxImages>-1 and i<len(tf.GetListOfKeys())-maxImages: continue
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print "Calc pedestal with event: ",name
            if rebin>1:
                obj.RebinX(rebin);
                obj.RebinY(rebin); 
            for ix in xrange(nx+1):
                for iy in xrange(ny+1):
                    x = obj.GetXaxis().GetBinCenter(ix+1)
                    y = obj.GetYaxis().GetBinCenter(iy+1)
                    if options.pedExclRegion and x>int(xmin) and x<int(xmax) and y>int(ymin) and y<int(ymax):
                        ix_rnd = iy_rnd = 0
                        # need to take only hits within the ellipse not to bias pedestals
                        while math.hypot(ix_rnd-self.xmax/2, iy_rnd-self.xmax/2)>600:
                            ix_rnd,iy_rnd = utilities.gen_rand_limit(int(xmin)/rebin,int(xmax)/rebin,
                                                                     int(ymin)/rebin,int(ymax)/rebin,
                                                                     maxx=self.xmax/rebin,maxy=self.xmax/rebin)
                        pedmap.Fill(x,y,obj.GetBinContent(ix_rnd+1,iy_rnd+1)/float(math.pow(self.rebin,2)))                        
                    else:
                        pedmap.Fill(x,y,obj.GetBinContent(ix+1,iy+1)/float(math.pow(self.rebin,2)))
        tf.Close()
        pedfile.cd()
        pedmap.Write()
        pedmean = ROOT.TH1D('pedmean','pedestal mean',500,97,103)
        pedrms = ROOT.TH1D('pedrms','pedestal RMS',500,0,5)
        for ix in xrange(nx):
            for iy in xrange(ny):
               pedmean.Fill(pedmap.GetBinContent(ix,iy)) 
               pedrms.Fill(pedmap.GetBinError(ix,iy)) 
        pedmean.Write()
        pedrms.Write()
        pedfile.Close()
        print "Pedestal calculated and saved into ",self.pedfile_name

    def reconstruct(self,evrange=(-1,-1,-1)):

        ROOT.gROOT.Macro('rootlogon.C')
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kWarning

        tf = ROOT.TFile.Open(self.rfile)
        #c1 = ROOT.TCanvas('c1','',600,600)
        ctools = cameraTools()
        print "Reconstructing event range: ",evrange[1],"-",evrange[2]
        # loop over events (pictures)
        for iobj,key in enumerate(tf.GetListOfKeys()) :
            iev = iobj if self.options.daq == 'btf' else iobj/2 # when PMT is present
            if self.options.maxEntries>0 and iev==max(evrange[0],0)+self.options.maxEntries: break
            if sum(evrange[1:])>-2:
                if iev<evrange[1] or iev>evrange[2]: continue
                
            name=key.GetName()
            obj=key.ReadObj()

            ###### DEBUG #########
            # if iev!=9 and iev!=4 and iev!=162: continue
            if iev<2: continue
            ######################
            
            if obj.InheritsFrom('TH2'):
                # DAQ convention
                # BTF convention
                if self.options.daq == 'btf':
                    run,event=(int(name.split('_')[0].split('run')[-1].lstrip("0")),int(name.split('_')[-1].lstrip("0")))
                else:
                    run,event=(int(name.split('_')[1].split('run')[-1].lstrip("0")),int(name.split('_')[-1].split('ev')[-1]))
                print "Processing run: ",run," event ",event,"..."
                self.outTree.fillBranch("run",run)
                self.outTree.fillBranch("event",event)

                pic_fullres = obj.Clone(obj.GetName()+'_fr')
                # rebin to make the cluster finding fast
                obj.RebinX(self.rebin); obj.RebinY(self.rebin)
                obj.Scale(1./float(math.pow(self.rebin,2)))

                # restrict to a subrange
                if self.options.pedExclRegion:
                    xmin,xmax = self.options.pedExclRegion.split(',')[0].split(':')
                    ymin,ymax = self.options.pedExclRegion.split(',')[1].split(':')
                    h2rs           = ctools.getRestrictedImage(obj,int(xmin),int(xmax),int(ymin),int(ymax))
                    pic_fullres_rs = ctools.getRestrictedImage(pic_fullres,int(xmin),int(xmax),int(ymin),int(ymax))
                    pedmap_fr_rs   = ctools.getRestrictedImage(self.pedmap_fr,int(xmin),int(xmax),int(ymin),int(ymax))
                else:
                    h2rs = obj
                    pic_fullres_rs = pic_fullres
                    pedmap_fr_rs = self.pedmap_fr
                # applying zero-suppression
                h2zs,h2unzs = ctools.zs(h2rs,self.pedmap,plot=False)
                #print "Zero-suppression done. Now clustering..."
                
                # Cluster reconstruction on 2D picture
                algo = 'DBSCAN'
                if self.options.type in ['beam','cosmics']: algo = 'HOUGH'
                snprod_inputs = {'picture': h2zs, 'pictureHD': pic_fullres_rs, 'pedmapHD': pedmap_fr_rs, 'name': name, 'algo': algo}
                snprod_params = {'snake_qual': 3, 'plot2D': True, 'plotpy': False, 'plotprofiles': True}
                snprod = SnakesProducer(snprod_inputs,snprod_params,self.options)
                snakes = snprod.run()                
                self.autotree.fillCameraVariables(h2zs,snakes)
                
                if self.options.daq != 'btf':
                   # PMT waveform reconstruction
                   from waveform import PeakFinder,PeaksProducer
                   wform = tf.Get('wfm_'+'_'.join(name.split('_')[1:]))
                   # sampling was 5 GHz (5/ns). Rebin by 5 (1/ns)
                   pkprod_inputs = {'waveform': wform}
                   pkprod_params = {'threshold': 0, # min threshold for a signal
                                    'minPeakDistance': 1, # number of samples (1 sample = 1ns )
                                    'prominence': 0.5, # noise after resampling very small
                                    'width': 1, # minimal width of the signal
                                    'resample': 5,  # to sample waveform at 1 GHz only
                                    'rangex': (6160,6500),
                                    'plotpy': False
                   }
                   pkprod = PeaksProducer(pkprod_inputs,pkprod_params,self.options)
                   peaksfinder = pkprod.run()
                   self.autotree.fillPMTVariables(peaksfinder,0.2*pkprod_params['resample'])
                
                # fill reco tree
                self.outTree.fill()

        ROOT.gErrorIgnoreLevel = savErrorLevel

                
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog h5file1,...,h5fileN [opts] ')
    parser.add_option("-o","--out", dest="outFile", type="string", default="reco.root", help="name of the output root file");
    parser.add_option('-j', '--jobs', dest='jobs', default=1, type='int', help='Jobs to be run in parallel')
    parser.add_option('-r', '--rebin', dest='rebin', default=4, type='int', help='Rebin factor (same in x and y)')
    parser.add_option(      '--numPedEvents', dest='numPedEvents', default=-1, type='float', help='Use the last n events to calculate the pedestal. Default is all events')
    parser.add_option(      '--max-entries', dest='maxEntries', default=-1, type='float', help='Process only the first n entries')
    parser.add_option(      '--pdir', dest='plotDir', default='./', type='string', help='Directory where to put the plots')
    parser.add_option('-p', '--pedestal', dest='justPedestal', default=False, action='store_true', help='Just compute the pedestals, do not run the analysis')
    parser.add_option(      '--exclude-region', dest='pedExclRegion', default=None, type='string', help='Exclude a rectangular region for pedestals. In the form "xmin:xmax,ymin:ymax"')
    parser.add_option(       '--daq', dest="daq", type="string", default="midas", help="DAQ type (btf/midas)");
    parser.add_option(       '--type', dest="type", type="string", default="neutrons", help="events type (beam/cosmics/neutrons)");
    
    (options, args) = parser.parse_args()

    inputf = args[0]

    if options.justPedestal:
        ana = analysis(inputf,options)
        print "Pedestals with rebin factor = ",options.rebin, "done. Exiting."
        sys.exit(0)
        
    ana = analysis(inputf,options)
    nev = ana.getNEvents() if options.maxEntries == -1 else int(options.maxEntries)
    print "This run has ",nev," events."
    print "Will save plots to ",options.plotDir
    
    if options.jobs>1:
        nj = int(nev/options.jobs)
        chunks = [(ichunk,i,min(i+nj-1,nev)) for ichunk,i in enumerate(xrange(0,nev,nj))]
        print chunks
        from multiprocessing import Pool
        pool = Pool(options.jobs)
        ret = pool.map(ana, chunks)
    else:
        ana.beginJob(options.outFile)
        ana.reconstruct()
        ana.endJob()
