#!/usr/bin/env python3.8
from multiprocessing import Pool,set_start_method

import os,math,sys,random
import numpy as np

import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import hist2array
from cameraChannel import cameraTools, cameraGeometry


from snakes import SnakesProducer
from output import OutputTree
from treeVars import AutoFillTreeProducer
import onlines3lib as s3on

# this kills also the still running subprocesses.
# use with a safe MAX TIMEOUT duration, since it will kill everything
def terminate_pool_2(pool):
    print ("Some subprocess timed out. Killing it brutally.")
    os.system('killall -9 python3.8')

# this still stucks
def terminate_pool(pool):
    print ("Some subprocess timed out. Killing it brutally.")
    for p in pool._pool:
        print ("KILLING PID ",p.pid)
        os.kill(p.pid, 9)
    pool.close()
    pool.terminate()
    pool.join()

import utilities
utilities = utilities.utils()

class analysis:

    def __init__(self,options):
        self.rebin = options.rebin        
        self.options = options
        self.pedfile_fullres_name = options.pedfile_fullres_name
        self.tmpname = options.tmpname
        geometryPSet   = open('modules_config/geometry_{det}.txt'.format(det=options.geometry),'r')
        geometryParams = eval(geometryPSet.read())
        self.cg = cameraGeometry(geometryParams)
        self.xmax = self.cg.npixx

        if not os.path.exists(self.pedfile_fullres_name):
            print("WARNING: pedestal file with full resolution ",self.pedfile_fullres_name, " not existing. First calculate them...")
            self.calcPedestal(options,1)
        if not options.justPedestal:
           print("Pulling pedestals...")
           # first the one for clustering with rebin
           ctools = cameraTools(self.cg)
           # then the full resolution one
           pedrf_fr = ROOT.TFile.Open(self.pedfile_fullres_name)
           self.pedmap_fr = pedrf_fr.Get('pedmap').Clone()
           self.pedmap_fr.SetDirectory(0)
           self.pedarr_fr = hist2array(self.pedmap_fr).T
           self.noisearr_fr = ctools.noisearray(self.pedmap_fr).T
           pedrf_fr.Close()
           if options.vignetteCorr:
               self.vignmap = ctools.loadVignettingMap()
           else:
               self.vignmap = np.ones((self.xmax, self.xmax))
            

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
        self.autotree = AutoFillTreeProducer(self.outTree,self.options.scfullinfo)

        self.outTree.branch("run", "I")
        self.outTree.branch("event", "I")
        self.outTree.branch("pedestal_run", "I")
        if self.options.Save_MC_data:
#            self.outTree.branch("MC_track_len","F")
            self.outTree.branch("eventnumber","I")
            self.outTree.branch("particle_type","I")
            self.outTree.branch("energy","F")
            self.outTree.branch("ioniz_energy","F")
            self.outTree.branch("drift","F")
            self.outTree.branch("phi_initial","F")
            self.outTree.branch("theta_initial","F")
            self.outTree.branch("MC_x_vertex","F")
            self.outTree.branch("MC_y_vertex","F")
            self.outTree.branch("MC_z_vertex","F")
            self.outTree.branch("MC_x_vertex_end","F")
            self.outTree.branch("MC_y_vertex_end","F")
            self.outTree.branch("MC_z_vertex_end","F")
            self.outTree.branch("MC_3D_pathlength","F")
            self.outTree.branch("MC_2D_pathlength","F")

        if self.options.camera_mode:
            self.autotree.createCameraVariables()
        #    self.autotree.createClusterVariables('cl')
            self.autotree.createClusterVariables('sc')
        #if self.options.pmt_mode:
        #    self.autotree.createPMTVariables()

    def endJob(self):
        self.outTree.write()
        self.outputFile.Close()
        
    def getNEvents(self):
        tf = ROOT.TFile.Open(self.tmpname)
        ret = int(len(tf.GetListOfKeys())/(options.waveform_number+1)) if (self.options.daq=='midas' and self.options.pmt_mode) else len(tf.GetListOfKeys())
        tf.Close()
        return ret

    def calcPedestal(self,options,alternativeRebin=-1):
        maxImages=options.maxEntries
        nx=ny=self.xmax
        rebin = self.rebin if alternativeRebin<0 else alternativeRebin
        nx=int(nx/rebin); ny=int(ny/rebin); 
        #pedfilename = 'pedestals/pedmap_ex%d_rebin%d.root' % (options.pedexposure,rebin)
        pedfilename = 'pedestals/pedmap_run%s_rebin%d.root' % (options.run,rebin)
        
        pedfile = ROOT.TFile.Open(pedfilename,'recreate')
        pedmap = ROOT.TH2D('pedmap','pedmap',nx,0,self.xmax,ny,0,self.xmax)
        pedmapS = ROOT.TH2D('pedmapsigma','pedmapsigma',nx,0,self.xmax,ny,0,self.xmax)

        pedsum = np.zeros((nx,ny))
        
        tf = ROOT.TFile.Open(self.tmpname)

        # first calculate the mean 
        numev = 0
        for i,e in enumerate(tf.GetListOfKeys()):
            iev = i if self.options.daq != 'midas' and self.options.pmt_mode else i/(options.waveform_number+1) # when PMT is present
            if iev in self.options.excImages: continue
            if maxImages>-1 and i>min(len(tf.GetListOfKeys()),maxImages): break
            
            name=e.GetName()
            obj=e.ReadObj()

            if not obj.InheritsFrom('TH2'): continue
            print("Calc pedestal mean with event: ",name)
            if rebin>1:
                obj.RebinX(rebin);
                obj.RebinY(rebin); 
            arr = hist2array(obj)
            pedsum = np.add(pedsum,arr)
            numev += 1
        pedmean = pedsum / float(numev)

        # now compute the rms (two separate loops is faster than one, yes)
        numev=0
        pedsqdiff = np.zeros((nx,ny))
        for i,e in enumerate(tf.GetListOfKeys()):
            iev = i if self.options.daq != 'midas' and self.options.pmt_mode else i/(options.waveform_number+1) # when PMT is present
            if iev in self.options.excImages: continue
            if maxImages>-1 and i>min(len(tf.GetListOfKeys()),maxImages): break
            
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print("Calc pedestal rms with event: ",name)
            if rebin>1:
                obj.RebinX(rebin);
                obj.RebinY(rebin); 
            arr = hist2array(obj)
            pedsqdiff = np.add(pedsqdiff, np.square(np.add(arr,-1*pedmean)))
            numev += 1
        pedrms = np.sqrt(pedsqdiff/float(numev-1))

        # now save in a persistent ROOT object
        for ix in range(nx):
            for iy in range(ny):
                pedmap.SetBinContent(ix+1,iy+1,pedmean[ix,iy]);
                pedmap.SetBinError(ix+1,iy+1,pedrms[ix,iy]);
                pedmapS.SetBinContent(ix+1,iy+1,pedrms[ix,iy]);
        tf.Close()

        pedfile.cd()
        pedmap.Write()
        pedmapS.Write()
        pedmean1D = ROOT.TH1D('pedmean','pedestal mean',500,97,103)
        pedrms1D = ROOT.TH1D('pedrms','pedestal RMS',500,0,5)
        for ix in range(nx):
            for iy in range(ny):
               pedmean1D.Fill(pedmap.GetBinContent(ix,iy)) 
               pedrms1D.Fill(pedmap.GetBinError(ix,iy)) 
        pedmean1D.Write()
        pedrms1D.Write()
        pedfile.Close()
        print("Pedestal calculated and saved into ",pedfilename)


    def reconstruct(self,evrange=(-1,-1,-1)):

        #ROOT.gROOT.Macro('rootlogon.C')
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kWarning

        tf = ROOT.TFile.Open(self.tmpname)
        #tf = ROOT.TFile.Open("https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/cygnus/Data/MAN/histograms_Run00677.root") #it would be nice if ot worked
        
        #c1 = ROOT.TCanvas('c1','',600,600)
        ctools = cameraTools(self.cg)
        print("Reconstructing event range: ",evrange[1],"-",evrange[2])
        # loop over events (pictures)
        for iobj,key in enumerate(tf.GetListOfKeys()) :
            iev = int(iobj/(options.waveform_number+1)) if self.options.daq == 'midas' and self.options.pmt_mode else iobj
            #print("max entries = ",self.options.maxEntries)
            if self.options.maxEntries>0 and iev==max(evrange[0],0)+self.options.maxEntries: break
            if sum(evrange[1:])>-2:
                if iev<evrange[1] or iev>evrange[2]: continue

            name=key.GetName()
            obj=key.ReadObj()

            # Routine to skip some images if needed
            if iev in self.options.excImages: continue

            if self.options.debug_mode == 1 and iev != self.options.ev: continue

            if obj.InheritsFrom('TH2'):
                if self.options.daq == 'btf':
                    run,event=(int(name.split('_')[0].split('run')[-1].lstrip("0")),int(name.split('_')[-1].lstrip("0")))
                elif self.options.daq == 'h5':
                    run,event=(int(name.split('_')[0].split('run')[-1]),int(name.split('_')[-1]))
                else:
                    run,event=(int(name.split('_')[1].split('run')[-1].lstrip("0")),int(name.split('_')[-1].split('ev')[-1]))
                print("Processing Run: ",run,"- Event ",event,"...")
                
                testspark=100*self.cg.npixx*self.cg.npixx+9000000
                if obj.Integral()>testspark:
                          print("Run ",run,"- Event ",event," has spark, will not be analyzed!")
                          continue
                            
                self.outTree.fillBranch("run",run)
                self.outTree.fillBranch("event",event)
                self.outTree.fillBranch("pedestal_run", int(self.options.pedrun))
                if self.options.Save_MC_data:
                    mc_tree = tf.Get('event_info/info_tree')
                    mc_tree.GetEntry(event)
#                    self.outTree.fillBranch("MC_track_len",mc_tree.MC_track_len)
                    self.outTree.fillBranch("eventnumber",mc_tree.eventnumber)
                    self.outTree.fillBranch("particle_type",mc_tree.particle_type)
                    self.outTree.fillBranch("energy",mc_tree.energy_ini)
                    self.outTree.fillBranch("ioniz_energy",mc_tree.ioniz_energy)
                    self.outTree.fillBranch("drift",mc_tree.drift)
                    self.outTree.fillBranch("phi_initial",mc_tree.phi_ini)
                    self.outTree.fillBranch("theta_initial",mc_tree.theta_ini)
                    self.outTree.fillBranch("MC_x_vertex",mc_tree.x_vertex)
                    self.outTree.fillBranch("MC_y_vertex",mc_tree.y_vertex)
                    self.outTree.fillBranch("MC_z_vertex",mc_tree.z_vertex)
                    self.outTree.fillBranch("MC_x_vertex_end",mc_tree.x_vertex_end)
                    self.outTree.fillBranch("MC_y_vertex_end",mc_tree.y_vertex_end)
                    self.outTree.fillBranch("MC_z_vertex_end",mc_tree.z_vertex_end)
                    self.outTree.fillBranch("MC_2D_pathlength",mc_tree.proj_track_2D)
                    self.outTree.fillBranch("MC_3D_pathlength",mc_tree.track_length_3D)

            if self.options.camera_mode:
                if obj.InheritsFrom('TH2'):
     
                    pic_fullres = obj.Clone(obj.GetName()+'_fr')
                    img_fr = hist2array(pic_fullres).T

                    # Upper Threshold full image
                    img_cimax = np.where(img_fr < self.options.cimax, img_fr, 0)
                    
                    # zs on full image + saturation correction on full image
                   # if self.options.saturation_corr:
                   # 	#print("you are in saturation correction mode")
                   # 	img_fr_sub = ctools.pedsub(img_cimax,self.pedarr_fr)
                   # 	img_fr_satcor = ctools.satur_corr(img_fr_sub) 
                   # 	img_fr_zs  = ctools.zsfullres(img_fr_satcor,self.noisearr_fr,nsigma=self.options.nsigma)
                   # 	img_rb_zs  = ctools.arrrebin(img_fr_zs,self.rebin)
                        
                   # # skip saturation and set satcor =img_fr_sub 
                   # else:
                   #     #print("you are in poor mode")
                    img_fr_sub = ctools.pedsub(img_cimax,self.pedarr_fr)
                   # img_fr_satcor = img_fr_sub  
                    img_fr_zs  = ctools.zsfullres(img_fr_sub,self.noisearr_fr,nsigma=self.options.nsigma)
                    img_rb_zs  = ctools.arrrebin(img_fr_zs,self.rebin)
                    
                    
                    # Cluster reconstruction on 2D picture
                    algo = 'DBSCAN'
                    #if self.options.type in ['beam','cosmics']: algo = 'HOUGH'
                    snprod_inputs = {'picture': img_rb_zs, 'pictureHD': img_fr_sub, 'picturezsHD': img_fr_zs, 'pictureOri': img_fr, 'vignette': self.vignmap, 'name': name, 'algo': algo}
                    plotpy = options.jobs < 2 # for some reason on macOS this crashes in multicore
                    snprod_params = {'snake_qual': 3, 'plot2D': False, 'plotpy': False, 'plotprofiles': False}
                    snprod = SnakesProducer(snprod_inputs,snprod_params,self.options,self.cg)
                    clusters,snakes = snprod.run()
                    self.autotree.fillCameraVariables(img_fr_zs)
                    self.autotree.fillClusterVariables(snakes,'sc')
                    #self.autotree.fillClusterVariables(clusters,'cl')
                    
            #if self.options.pmt_mode:
            #    if obj.InheritsFrom('TGraph'):
            #        # PMT waveform reconstruction
            #        from waveform import PeakFinder,PeaksProducer
            #        wform = tf.Get('wfm_'+'_'.join(name.split('_')[1:]))
            #        # sampling was 5 GHz (5/ns). Rebin by 5 (1/ns)
            #        pkprod_inputs = {'waveform': wform}
            #        pkprod_params = {'threshold': options.threshold, # min threshold for a signal (baseline is -20 mV)
            #                         'minPeakDistance': options.minPeakDistance, # number of samples (1 sample = 1ns )
            #                         'prominence': options.prominence, # noise after resampling very small
            #                         'width': options.width, # minimal width of the signal
            #                         'resample': options.resample,  # to sample waveform at 1 GHz only
            #                         'rangex': (self.options.time_range[0],self.options.time_range[1]),
            #                         'plotpy': options.pmt_plotpy
            #        }
            #        pkprod = PeaksProducer(pkprod_inputs,pkprod_params,self.options)
            #        
            #        peaksfinder = pkprod.run()
            #        self.autotree.fillPMTVariables(peaksfinder,0.2*pkprod_params['resample'])
            #        
            # fill reco tree (just once/event, and the TGraph is analyses as last)
            if (self.options.daq == 'midas' and self.options.pmt_mode):
                #if obj.InheritsFrom('TGraph'):			#works only if you have one waveform
                nome = obj.GetName()
                if "ch8" in nome:			#if you have mesh as last waveform use "mesh" or "wfm".. look first at the data
                    self.outTree.fill()
            else:
                self.outTree.fill()

        ROOT.gErrorIgnoreLevel = savErrorLevel

                
if __name__ == '__main__':
    from optparse import OptionParser
    
    parser = OptionParser(usage='%prog h5file1,...,h5fileN [opts] ')
    parser.add_option('-r', '--run', dest='run', default='00000', type='string', help='run number with 5 characteres')
    parser.add_option('-j', '--jobs', dest='jobs', default=1, type='int', help='Jobs to be run in parallel (-1 uses all the cores available)')
    parser.add_option(      '--max-entries', dest='maxEntries', default=-1, type='float', help='Process only the first n entries')
    parser.add_option(      '--pdir', dest='plotDir', default='./', type='string', help='Directory where to put the plots')
    parser.add_option(      '--tmp',  dest='tmpdir', default=None, type='string', help='Directory where to put the input file. If none is given, /tmp/<user> is used')
    parser.add_option(      '--max-hours', dest='maxHours', default=-1, type='float', help='Kill a subprocess if hanging for more than given number of hours.')
    
    (options, args) = parser.parse_args()
    
    f = open(args[0], "r")
    params = eval(f.read())
    
    for k,v in params.items():
        setattr(options,k,v)

    run = int(options.run)
    
    if options.debug_mode == 1:
        setattr(options,'outFile','reco_run%d_%s_debug.root' % (run, options.tip))
        if options.ev: options.maxEntries = options.ev + 1
        #if options.daq == 'midas': options.ev +=0.5 
    else:
        setattr(options,'outFile','reco_run%05d_%s.root' % (run, options.tip))
        
    if not hasattr(options,"pedrun"):
        pf = open("pedestals/pedruns.txt","r")
        peddic = eval(pf.read())
        options.pedrun = -1
        for runrange,ped in peddic.items():
            if int(runrange[0])<=run<=int(runrange[1]):
                options.pedrun = int(ped)
                print("Will use pedestal run %05d, valid for run range [%05d - %05d]" % (int(ped), int(runrange[0]), (runrange[1])))
                break
        assert options.pedrun>0, ("Didn't find the pedestal corresponding to run ",run," in the pedestals/pedruns.txt. Check the dictionary inside it!")
            
    setattr(options,'pedfile_fullres_name', 'pedestals/pedmap_run%s_rebin1.root' % (options.pedrun))

    USER = os.environ['USER']
    tmpdir = '/mnt/ssdcache/' if os.path.exists('/mnt/ssdcache/') else '/tmp/'
    os.system('mkdir -p {tmpdir}/{user}'.format(tmpdir=tmpdir,user=USER))
    if os.path.isfile("%s/%s/histograms_Run%05d.root" % (tmpdir,USER,int(options.run))):
        options.tmpname = "%s/%s/histograms_Run%05d.root" % (tmpdir,USER,int(options.run))
    else:
        print ('Downloading file: ' + s3on.s3_root_file(options.tag, int(options.run)))
        options.tmpname = s3on.s3_download_root_file(s3on.s3_root_file(options.tag, int(options.run)),int(options.run))
    
    if options.justPedestal:
        ana = analysis(options)
        print("Pedestals done. Exiting.")
        if options.donotremove == False:
            os.remove(options.tmpname)
            print("tmp file removed")
        sys.exit(0)     
    
    ana = analysis(options)
    nev = ana.getNEvents() if options.maxEntries == -1 else int(options.maxEntries)
    print("This run has ",nev," events.")
    print("Will save plots to ",options.plotDir)
    os.system('cp utils/index.php {od}'.format(od=options.plotDir))
    
    nThreads = 1
    if options.jobs==-1:
        import multiprocessing
        nThreads = multiprocessing.cpu_count()
    else:
        nThreads = options.jobs

    if nThreads>1:
        print ("RUNNING USING ",nThreads," THREADS.")
        nj = int(nev/nThreads)
        chunks = [(ichunk,i,min(i+nj-1,nev)) for ichunk,i in enumerate(range(0,nev,nj))]
        print(chunks)
        pool = Pool(nThreads)
        ret = list(pool.apply_async(ana,args=(c, )) for c in chunks)
        try:
            if options.maxHours>0:
                maxTime = options.maxHours * 3600
            else:
                # 64-bit integer, converted from nanoseconds to seconds, and subtracting 0.1 just to be in bounds.
                maxTime = 2 ** 63 / 1e9 - 0.1
            print([r.get(timeout=maxTime) for r in ret])
            pool.close()
            pool.terminate()
            pool.join()
        except TimeoutError:
            print("except")
            terminate_pool_2(pool)
        print("Now hadding the chunks...")
        base = options.outFile.split('.')[0]
        os.system('{rootsys}/bin/hadd -k -f {base}.root {base}_chunk*.root'.format(rootsys=os.environ['ROOTSYS'],base=base))
        os.system('rm {base}_chunk*.root'.format(base=base))
    else:
        ana.beginJob(options.outFile)
        ana.reconstruct()
        ana.endJob()

    #### FOR SOME REASON THIS DOESN'T WORK IN BATCH.
    # now add the git commit hash to track the version in the ROOT file
    # tf = ROOT.TFile.Open(options.outFile,'update')
    # githash = ROOT.TNamed("gitHash",str(utilities.get_git_revision_hash()).replace('\n',''))
    # githash.Write()
    # tf.Close()
    
    #if options.donotremove == False:
    #     os.remove(options.tmpname)
    #     print("tmp file removed")
