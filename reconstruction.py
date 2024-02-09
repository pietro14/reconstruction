from concurrent import futures
from subprocess import Popen, PIPE
import signal,time

import os,math,sys,random,re,gc
import numpy as np

import ROOT
ROOT.gROOT.SetBatch(True)
import uproot
from cameraChannel import cameraTools, cameraGeometry
import midas.file_reader
import h5py


from snakes import SnakesProducer
from output import OutputTree
from treeVars import AutoFillTreeProducer
import swiftlib as sw
import cygno as cy

import pandas as pd

import utilities
utilities = utilities.utils()

from waveform import PMTreco

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

        eventContentPSet = open('modules_config/reco_eventcontent.txt')
        self.eventContentParams = eval(eventContentPSet.read())
        for k,v in self.eventContentParams.items():
            setattr(self.options,k,v)
        
        if not os.path.exists(self.pedfile_fullres_name):
            print("WARNING: pedestal file with full resolution ",self.pedfile_fullres_name, " not existing. First calculate them...")
            self.calcPedestal(options,1)
        if not options.justPedestal:
           print("Pulling pedestals...")
           # first the one for clustering with rebin
           ctools = cameraTools(self.cg)
           # then the full resolution one
           pedrf_fr = uproot.open(self.pedfile_fullres_name)
           self.pedarr_fr   = pedrf_fr['pedmap'].values().T
           self.noisearr_fr = pedrf_fr['pedmap'].errors().T
           if options.vignetteCorr:
               self.vignmap = ctools.loadVignettingMap()
           else:
               self.vignmap = np.ones((self.xmax, self.xmax))

        ## Dictionary with the PMT parameters found in config_file
        self.pmt_params = {
        'threshold': options.threshold,
        'height_RMS': options.height_RMS, 
        'minPeakDistance': options.minPeakDistance, 
        'prominence': options.prominence,
        'fixed_prom': options.fixed_prom, 
        'width': options.width,
        'resample': options.resample,
        'plotpy': options.pmt_plotpy,
        'wf_in_tree': options.pmt_wf_in_tree,
        'digit_tag': options.digitizer_tag,
        'pmt_verb':  options.pmt_verbose
        }

    # the following is needed for multithreading
    def __call__(self,evrange=(-1,-1,-1)):
        if evrange[0]==-1:
            outfname = '{outdir}/{base}'.format(base=self.options.outFile,outdir=options.outdir)
        else:
            outfname = '{outdir}/{base}_chunk{ij}.root'.format(base=self.options.outFile.split('.')[0],ij=evrange[0],outdir=self.options.outdir)
        self.beginJob(outfname)
        self.reconstruct(evrange)
        self.endJob()
        
    def beginJob(self,outfname):
        # prepare output file
        ROOT.EnableThreadSafety()
        self.outputFile = ROOT.TFile.Open(outfname, "RECREATE")
        print("Opening out file: ",outfname," self.outputFile = ",self.outputFile)
        ROOT.gDirectory.cd()

        # prepare output tree
        self.outputTree = ROOT.TTree("Events","Tree containing reconstructed quantities")
        self.outTree = OutputTree(self.outputFile,self.outputTree)
        self.autotree = AutoFillTreeProducer(self.outTree,self.eventContentParams)

        ## Prepare PMT waveform Tree (1 event = 1 waveform)
        self.outputTree_pmt = ROOT.TTree("PMT_Events","Tree containing reconstructed PMT quantities")
        self.outTree_pmt = OutputTree(self.outputFile,self.outputTree_pmt)
        self.autotree_pmt = AutoFillTreeProducer(self.outTree_pmt,self.eventContentParams)

        ## Prepare PMT average waveform Tree (1 event = 1 averaged waveform using 4 PMTs)
        self.outputTree_pmt_avg = ROOT.TTree("PMT_Avg_Events","Tree containing the average PMT waveforms of 4 channels")
        self.outTree_pmt_avg = OutputTree(self.outputFile,self.outputTree_pmt_avg)
        self.autotree_pmt_avg = AutoFillTreeProducer(self.outTree_pmt_avg,self.eventContentParams)

        self.outTree.branch("run", "I", title="run number")
        self.outTree.branch("event", "I", title="event number")
        self.outTree.branch("pedestal_run", "I", title="run number used for pedestal subtraction")

        if self.options.save_MC_data:
            #self.outTree.branch("MC_track_len","F")
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
            self.autotree.createTimeCameraVariables()
            self.autotree.createClusterVariables('sc')
            if self.options.cosmic_killer:
                self.autotree.addCosmicKillerVariables('sc')
        
        ## Create PMT branchs
        if self.options.pmt_mode:
            self.autotree_pmt.createPMTVariables(self.pmt_params)                   ## Individual waveform
            self.autotree_pmt_avg.createPMTVariables_average(self.pmt_params)       ## Average Waveform
            
            # time variables -- see if it can be merged with above functions
            self.autotree_pmt.createTimePMTVariables()                          
            self.autotree_pmt_avg.createTimePMTVariables_average()

        if options.environment_variables: 
            self.autotree.createEnvVariables()
        

    def endJob(self):
        self.outTree.write()
        self.outTree_pmt.write()
        self.outTree_pmt_avg.write()
        self.outputFile.Close()
        
    def getNEvents(self,options):
        if options.rawdata_tier == 'root':
            tf = sw.swift_read_root_file(self.tmpname)
            pics = [k for k in tf.keys() if 'pic' in k]
            return len(pics)
        elif options.rawdata_tier == 'h5':
            tf = sw.swift_read_h5_file(self.tmpname)
            pics = [k for k in tf.keys() if 'pic' in k]
            print("n events:", len(pics))
            return len(pics)
        else:
            run,tmpdir,tag = self.tmpname
            mf = sw.swift_download_midas_file(run,tmpdir,tag)
            evs =0
            for mevent in mf:
                if mevent.header.is_midas_internal_event():
                    continue
                else:
                    keys = mevent.banks.keys()
                for iobj,key in enumerate(keys):
                    name=key
                    if name.startswith('CAM'):
                        evs += 1
            return evs

    def calcPedestal(self,options,alternativeRebin=-1):
        maxImages=options.maxEntries
        nx=ny=self.xmax
        rebin = self.rebin if alternativeRebin<0 else alternativeRebin
        nx=int(nx/rebin); ny=int(ny/rebin); 
        pedfilename = 'pedestals/pedmap_run%s_rebin%d.root' % (options.pedrun,rebin)
        
        pedfile = ROOT.TFile.Open(pedfilename,'recreate')
        pedmap = ROOT.TH2D('pedmap','pedmap',nx,0,self.xmax,ny,0,self.xmax)
        pedmapS = ROOT.TH2D('pedmapsigma','pedmapsigma',nx,0,self.xmax,ny,0,self.xmax)

        pedsum = np.zeros((nx,ny))

        if options.rawdata_tier == 'root':
            tmpdir = '{tmpdir}'.format(tmpdir=options.tmpdir if options.tmpdir else "/tmp/")
            if not sw.checkfiletmp(int(options.pedrun),'root',tmpdir):
                print ('Downloading file: ' + sw.swift_root_file(options.tag, int(options.pedrun)))
                pedfilename = sw.swift_download_root_file(sw.swift_root_file(options.tag, int(options.pedrun)),int(options.pedrun),tmpdir)
            else:
                pedfilename = sw.swift_download_root_file(sw.swift_root_file(options.tag, int(options.pedrun)),int(options.pedrun),tmp=tmpdir,justName=True)                
            tf = sw.swift_read_root_file(pedfilename)
            keys = tf.keys()
            mf = [0] # dummy array to make a common loop with MIDAS case
        else:
            sigrun,tmpdir,tag = self.tmpname
            mf = sw.swift_download_midas_file(options.pedrun,tmpdir,tag)
            #mf = self.tmpname

        # first calculate the mean 
        numev = 0
        if  options.rawdata_tier == 'midas':
            mf.jump_to_start()
            for mevent in mf:
                if mevent.header.is_midas_internal_event():
                    continue
                else:
                    keys = mevent.banks.keys()
                for iobj,key in enumerate(keys):
                    name=key
                    if name.startswith('CAM'):
                        if options.rawdata_tier == 'root':
                            arr = tf[key].values()
                        else:
                            arr,_,_ = cy.daq_cam2array(mevent.banks[key])
                            arr = np.rot90(arr)
                        justSkip=False
                        if numev in self.options.excImages: justSkip=True
                        if maxImages>-1 and numev>min(len(keys),maxImages): break
                        if numev>100: break # no need to compute pedestals with >100 evts (avoid large RAM usage)
                            
                        if numev%20 == 0:
                            print("Calc pedestal mean with event: ",numev)
                        if justSkip:
                             continue
                        if rebin>1:
                            ctools.arrrebin(arr,rebin)
                        pedsum = np.add(pedsum,arr)
                        numev += 1
        else:
            print ("keys = ",keys)
            for i,name in enumerate(keys):
                if 'pic' in name:
                    patt = re.compile('\S+run(\d+)_ev(\d+)')
                    m = patt.match(name)
                    run = int(m.group(1))
                    event = int(m.group(2))
                justSkip=False
                if event in self.options.excImages: justSkip=True
                if maxImages>-1 and event>min(len(keys),maxImages): break
                if numev>100: break # no need to compute pedestals with >100 evts (avoid large RAM usage)
                if 'pic' not in name: justSkip=True
                if justSkip:
                    continue
                if event%20 == 0:
                    print("Calc pedestal mean with event: ",event)
                arr = tf[name].values()
                pedsum = np.add(pedsum,arr)
                numev += 1
        pedmean = pedsum / float(numev)

        # now compute the rms (two separate loops is faster than one, yes)
        numev=0
        pedsqdiff = np.zeros((nx,ny))
        numev = 0
        if  options.rawdata_tier == 'midas':
            mf.jump_to_start()
            for mevent in mf:
                if  options.rawdata_tier == 'midas':
                    if mevent.header.is_midas_internal_event():
                        continue
                    else:
                        keys = mevent.banks.keys()
                for iobj,key in enumerate(keys):
                    name=key
                    if name.startswith('CAM'):
                        if options.rawdata_tier == 'root':
                            arr = tf[key].values()
                        else:
                            arr,_,_ = cy.daq_cam2array(mevent.banks[key])
                            arr = np.rot90(arr)
                        justSkip=False
                        if numev in self.options.excImages: justSkip=True
                        if maxImages>-1 and numev>min(len(keys),maxImages): break
                        if numev>100: break # no need to compute pedestals with >100 evts (avoid large RAM usage)
             
                        if numev%20 == 0:
                            print("Calc pedestal rms with event: ",numev)
                        if justSkip:
                             continue
                        if rebin>1:
                            ctools.arrrebin(arr,rebin)
                        pedsqdiff = np.add(pedsqdiff, np.square(np.add(arr,-1*pedmean)))
                        numev += 1
        else:
            for i,name in enumerate(keys):
                if 'pic' in name:
                    patt = re.compile('\S+run(\d+)_ev(\d+)')
                    m = patt.match(name)
                    run = int(m.group(1))
                    event = int(m.group(2))
                justSkip=False
                if event in self.options.excImages: justSkip=True
                if maxImages>-1 and event>min(len(keys),maxImages): break
                if numev>100: break # no need to compute pedestals with >100 evts (avoid large RAM usage)
     
                if 'pic' not in name: justSkip=True
                if justSkip:
                     continue
                if event%20 == 0:
                    print("Calc pedestal rms with event: ",event)
                arr = tf[name].values()
                pedsqdiff = np.add(pedsqdiff, np.square(np.add(arr,-1*pedmean)))
                numev += 1
        pedrms = np.sqrt(pedsqdiff/float(numev-1))

        # now save in a persistent ROOT object
        for ix in range(nx):
            for iy in range(ny):
                pedmap.SetBinContent(ix+1,iy+1,pedmean[ix,iy]);
                pedmap.SetBinError(ix+1,iy+1,pedrms[ix,iy]);
                pedmapS.SetBinContent(ix+1,iy+1,pedrms[ix,iy]);

        pedfile.cd()
        pedmap.Write()
        pedmapS.Write()
        pedmean1D = ROOT.TH1D('pedmean','pedestal mean',500,97,103)
        pedrms1D = ROOT.TH1D('pedrms','pedestal RMS',1000,0,10)
        for ix in range(nx):
            for iy in range(ny):
               pedmean1D.Fill(pedmap.GetBinContent(ix,iy)) 
               pedrms1D.Fill(pedmap.GetBinError(ix,iy)) 
        pedmean1D.Write()
        pedrms1D.Write()
        pedfile.Close()
        print("Pedestal calculated and saved into ",pedfilename)


    def reconstruct(self,evrange=(-1,-1,-1)):

        ROOT.gROOT.Macro('rootlogon.C')
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kWarning
        
        ctools = cameraTools(self.cg)
        print("Reconstructing event range: ",evrange[1],"-",evrange[2])
        self.outputFile.cd()
        
        if self.options.rawdata_tier == 'root':
            tf = sw.swift_read_root_file(self.tmpname)
            keys = tf.keys()
            mf = [0] # dummy array to make a common loop with MIDAS case
        elif self.options.rawdata_tier == 'h5':
            tf = sw.swift_read_h5_file(self.tmpname)
            keys = tf.keys()
            mf = [0] # dummy array to make a common loop with MIDAS case

        elif options.rawdata_tier == 'midas':
            run,tmpdir,tag = self.tmpname
            mf = sw.swift_download_midas_file(run,tmpdir,tag)
            
            ## Necessary to read the ODB to retrieve some info necessary for the waveform analysis
            ## Seems to repeat the opening process but *doesn't* slow down the code.
            if self.options.pmt_mode == 1:        
                
                if run > 7790:   ## There is an issue with Run1 because it doesn't contain all the info present in Run2
                        #TO FIX -- ADD TAG CHECK since this only applies for LIME @ LNGS
                    odb=cy.get_bor_odb(mf)
                    corrected  = odb.data['Configurations']['DRS4Correction']
                    channels_offsets  = odb.data['Configurations']['DigitizerOffset']
                    camera_exposure   = odb.data['Configurations']['Exposure']
                    mf.jump_to_start()
                else:
                    corrected = True
                    camera_exposure = 300
                    channels_offsets = 0
                    mf.jump_to_start()


            # mf.jump_to_start()
            dslow = pd.DataFrame()
            if options.environment_variables:
        
                odb = cy.get_bor_odb(mf)
                header_environment = odb.data['Equipment']['Environment']['Settings']['Names Input']
                value_variables = odb.data['Equipment']['Environment']['Variables']
                dslow = pd.DataFrame(columns = header_environment)
                dslow.loc[len(dslow)] = value_variables['Input']
                for i in dslow.keys():
                    dslow = utilities.conversion_env_variables(dslow, odb, i, j = 0)
                self.autotree.fillEnvVariables(dslow.take([0]))
                j = 1

        numev = 0
        event=0

        for mevent in mf:
            if self.options.rawdata_tier == 'midas':
                if mevent.header.is_midas_internal_event():
                    continue
                else:
                    keys = mevent.banks.keys()
                    
            for bank_name, bank in mevent.banks.items():
                name=bank_name
                camera = False
                pmt = False

                if self.options.rawdata_tier == 'root':
                    if 'pic' in name:
                        patt = re.compile('\S+run(\d+)_ev(\d+)')
                        m = patt.match(name)
                        run = int(m.group(1))
                        event = int(m.group(2))
                        obj = tf[key].values()
                        #obj = np.rot90(obj)
                        camera=True
                    else:
                        camera=False
                elif self.options.rawdata_tier == 'h5':
                    if 'pic' in name:
                        patt = re.compile('\S+run(\d+)_ev(\d+)')
                        m = patt.match(name)
                        run = int(m.group(1))
                        event = int(m.group(2))
                        obj = np.array(tf[key])
                        #obj = np.rot90(obj)
                        camera=True
                    else:
                        camera=False

                elif self.options.rawdata_tier == 'midas':
                    
                    run = int(self.options.run)

                    if bank_name=='CAM0' and options.camera_mode:

                        obj,_,_ = cy.daq_cam2array(bank, dslow)
                        obj = np.rot90(obj)
                    
                        camera=True
                    
                    elif name.startswith('INPT') and options.environment_variables: # SLOW channels array
                        #try:
                        dslow = utilities.read_env_variables(mevent.banks[key], dslow, odb, j=j)
                           #print(dslow)
                        self.autotree.fillEnvVariables(dslow.take([j]))
                        j = j+1
                           #print(dslow)
                        #except:
                        #   print("WARNING: INPT bank is not as expected.")
                    
                    elif bank_name=='DGH0' and options.pmt_mode:
                        header=cy.daq_dgz_full2header(bank, verbose=False)
                        SIC = header.SIC
                        sample_rate = header.sampling_rate

                        nChannels_f  = header.nchannels[0]
                        nTriggers_f = len(header.TTT[0])
                        TTTs_f = header.TTT[0]

                        nChannels_s  = header.nchannels[1]
                        nTriggers_s = len(header.TTT[1])
                        TTTs_s = header.TTT[1]

                        ## TO FIX : catch the dgitizer general 'tag' in the config file. Should work
                        waveform_f, waveform_s = cy.daq_dgz_full2array(mevent.banks['DIG0'], header, verbose=False, corrected=corrected, ch_offset=channels_offsets,tag=self.pmt_params['digit_tag'])

                        pmt = True
                    else:
                        camera = False
                        pmt = False

                    event=numev

                justSkip = False
                if event<evrange[1]: justSkip=True
                if event>evrange[2]: return # avoids seeking up to EOF which with MIDAS is slow
                if event in self.options.excImages: justSkip=True
                if self.options.debug_mode == 1 and event != self.options.ev: justSkip=True
                if justSkip:
                    continue
                
                # Event is identified by EITHER camera or pmt
                # FIX: Perhaps we should change this with the same scheme as the filling at the end
                # FIX: I'm filling the Branchs twice, but that should not be a problem (I think it overwrites)
                if camera == True or pmt == True:
                    print("Processing Run: ",run,"- Event ",event,"...")

                    self.outTree.fillBranch("run",run)
                    self.outTree.fillBranch("event",event)
                    self.outTree.fillBranch("pedestal_run", int(self.options.pedrun))

                if self.options.camera_mode:
                    if camera==True:

                        testspark=2*100*self.cg.npixx*self.cg.npixx+9000000		#for ORCA QUEST data multiply also by 2: 2*100*....
                        if np.sum(obj)>testspark:
                            print("Run ",run,"- Event ",event," has spark, will not be analyzed!")
                            continue

                        if self.options.save_MC_data:
                            mc_tree = tf.Get('event_info/info_tree')
                            mc_tree.GetEntry(event)
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

             
                        img_fr = obj.T
         
                        # Upper Threshold full image
                        img_cimax = np.where(img_fr < self.options.cimax, img_fr, 0)
                        
                        # zs on full image + saturation correction on full image
                        if self.options.saturation_corr:
                            #print("you are in saturation correction mode")
                            t_pre0 = time.perf_counter()
                            img_fr_sub = ctools.pedsub(img_cimax,self.pedarr_fr)
                            t_pre1 = time.perf_counter()
                            img_fr_satcor = ctools.satur_corr(img_fr_sub) 
                            t_pre2 = time.perf_counter()
                            img_fr_zs  = ctools.zsfullres(img_fr_satcor,self.noisearr_fr,nsigma=self.options.nsigma)
                            t_pre3 = time.perf_counter()
                            img_fr_zs_acc = ctools.acceptance(img_fr_zs,self.cg.ymin,self.cg.ymax,self.cg.xmin,self.cg.xmax)
                            t_pre4 = time.perf_counter()
                            img_rb_zs  = ctools.arrrebin(img_fr_zs_acc,self.rebin)
                            t_pre5 = time.perf_counter()
                            
                        # skip saturation and set satcor =img_fr_sub 
                        else:
                            #print("you are in poor mode")
                            t_pre0 = time.perf_counter()
                            img_fr_sub = ctools.pedsub(img_cimax,self.pedarr_fr)
                            t_pre1 = time.perf_counter()
                            img_fr_satcor = img_fr_sub  
                            t_pre2 = time.perf_counter()
                            img_fr_zs  = ctools.zsfullres(img_fr_satcor,self.noisearr_fr,nsigma=self.options.nsigma)
                            t_pre3 = time.perf_counter()
                            img_fr_zs_acc = ctools.acceptance(img_fr_zs,self.cg.ymin,self.cg.ymax,self.cg.xmin,self.cg.xmax)
                            t_pre4 = time.perf_counter()
                            img_rb_zs  = ctools.arrrebin(img_fr_zs_acc,self.rebin)
                            t_pre5 = time.perf_counter()

                        t_pedsub = t_pre1 - t_pre0
                        t_saturation = t_pre2 - t_pre1
                        t_zerosup = t_pre3 - t_pre2
                        t_xycut = t_pre4 - t_pre3
                        t_rebin = t_pre5 - t_pre4
                            
                        # Cluster reconstruction on 2D picture
                        algo = 'DBSCAN'
                        if self.options.type in ['beam','cosmics']: algo = 'HOUGH'
                        snprod_inputs = {'picture': img_rb_zs, 'pictureHD': img_fr_satcor, 'picturezsHD': img_fr_zs, 'pictureOri': img_fr, 'vignette': self.vignmap, 'name': name, 'algo': algo}
                        plotpy = self.options.jobs < 2 # for some reason on macOS this crashes in multicore
                        snprod_params = {'snake_qual': 3, 'plot2D': False, 'plotpy': False, 'plotprofiles': False}
                        t_DBSCAN_0 = time.perf_counter()
                        snprod = SnakesProducer(snprod_inputs,snprod_params,self.options,self.cg)
                        t_DBSCAN_1 = time.perf_counter()
                        snakes, t_DBSCAN, t_variables, lp_len, t_medianfilter, t_noisered = snprod.run()
                        t_DBSCAN_2 = time.perf_counter()
                        if options.debug_mode == 1:
                            print(f"1. DBSCAN run + variables calculation in {t_DBSCAN_2 - t_DBSCAN_1:0.4f} seconds")
                        self.autotree.fillCameraVariables(img_fr_zs)
                        t_DBSCAN_3 = time.perf_counter()
                        if options.debug_mode == 1:
                            print(f"fillCameraVariables in {t_DBSCAN_3 - t_DBSCAN_2:0.4f} seconds")
                        self.autotree.fillClusterVariables(snakes,'sc')
                        t_DBSCAN_4 = time.perf_counter()
                        self.autotree.fillTimeCameraVariables(t_variables, t_DBSCAN, lp_len, t_pedsub, t_saturation, t_zerosup, t_xycut, t_rebin, t_medianfilter, t_noisered)
                        if options.debug_mode == 1:
                            print(f"fillClusterVariables in {t_DBSCAN_4 - t_DBSCAN_3:0.04f} seconds")
                            print()
                        del img_fr_sub,img_fr_satcor,img_fr_zs,img_fr_zs_acc,img_rb_zs
         
                if self.options.pmt_mode:
                    if pmt == True:

                        chs_to_analyse = 4
                        fast_sampling = 1024
                        slow_sampling = 4000

                        ## Fast waveforms
                        if options.debug_mode == 1:
                            print("Number of fast triggers: {}".format(nTriggers_f))

                        for trg in range(nTriggers_f):    

                            insideGE = 0
                            # Method 1 uses a specific signal in ch5 to check if inside GE or not. Now *deprecated*
                            # Method 2 uses the TTTs to check this condition
                            if (TTTs_f[trg] * 8.5/1000/1000) >= 180 and (TTTs_f[trg] * 8.5/1000/1000) <= (camera_exposure*1000):
                                insideGE = 1

                            # Prepare the weighted average waveform 
                            sing_weig_avg_fast_wf = [ [0]*fast_sampling for _ in range(chs_to_analyse)]  
                            fast_wf_weights_snr = [0] * chs_to_analyse
                            weight_average_wf = [0]* fast_sampling

                            for chs in range(chs_to_analyse):

                                ch = chs + 1
                                indx = trg * nChannels_f + ch
                                waveform_info = { 'run' : run, 'event': event, 'channel' : ch, 'trigger' : trg , 'GE' : insideGE, 'sampling' : "fast", 'TTT' : (TTTs_f[trg]*8.5/1000./1000.)}

                                t0 = time.perf_counter()

                                fast_waveform = PMTreco(waveform_info, waveform_f[indx], self.pmt_params)
                                fast_waveform.__repr__()

                                t1 = time.perf_counter()
                                t_waveforms = t1 - t0

                                self.autotree_pmt.fillPMTVariables(fast_waveform) 
                                self.autotree_pmt.fillTimePMTVariables(t_waveforms)
                                self.outTree_pmt.fill()

                                # Weighted averaged waveform (weight = SNR)
                                snr_ratio = fast_waveform.getSignalToNoise()
                                fast_wf_weights_snr[ch-1] = snr_ratio
                                sing_weig_avg_fast_wf[ch-1] = waveform_f[indx]

                                # If one wants to visualize the new weighted waveforms,
                                # meaning how much they actual weight for the final average, 
                                # Ask David for the script changes

                                if ch == 4:

                                    fast_wf_weights_snr = [ (x / max(fast_wf_weights_snr)) for x in fast_wf_weights_snr ]   # Normalization of the weights

                                    for k in range(chs_to_analyse):

                                        for j in range(fast_sampling):

                                            weight_average_wf[j] += sing_weig_avg_fast_wf[k][j]*fast_wf_weights_snr[k]/sum(fast_wf_weights_snr)

                                    waveform_info_fast_wei_avg = { 'run' : run, 'event': event, 'channel' : 9, 'trigger' : trg, 'GE' : 9, 'sampling' : "fast"}

                                    t0 = time.perf_counter()
                                    fast_waveform_wei_avg = PMTreco(waveform_info_fast_wei_avg, weight_average_wf, self.pmt_params)
                                    t1 = time.perf_counter()
                                    t_waveforms = t1 - t0 

                                    self.autotree_pmt_avg.fillPMTVariables_average(fast_waveform_wei_avg)
                                    self.autotree_pmt_avg.fillTimePMTVariables_average(t_waveforms)
                                    #fast_waveform_wei_avg.__repr__()           ## Verbose of averaged waveform
                                    self.outTree_pmt_avg.fill()

                                    del waveform_info_fast_wei_avg
                                    del fast_waveform_wei_avg

                                del waveform_info
                                del fast_waveform
                        # Slow waveforms
                        if options.debug_mode == 1:
                            print("Number of slow triggers: {}".format(nTriggers_s))

                        for trg in range(nTriggers_s):    

                            insideGE = 0
                            if (TTTs_f[trg] * 8.5/1000/1000) >= 180 and (TTTs_f[trg] * 8.5/1000/1000) <= (camera_exposure*1000):
                                insideGE = 1

                            sing_weig_avg_slow_wf = [ [0]* slow_sampling for _ in range(chs_to_analyse)]  
                            slow_wf_weights_snr = [0] * chs_to_analyse
                            weight_average_wf = [0]* slow_sampling

                            for chs in range(chs_to_analyse):

                                ch = chs + 1
                                indx = trg * nChannels_s + ch
                                waveform_info = { 'run' : run, 'event': event, 'channel' : ch, 'trigger' : trg , 'GE' : insideGE , 'sampling' : "slow", 'TTT' : (TTTs_s[trg]*8.5/1000./1000.)}
                                
                                t0 = time.perf_counter()

                                slow_waveform = PMTreco(waveform_info, waveform_s[indx], self.pmt_params)
                                slow_waveform.__repr__()

                                t1 = time.perf_counter()
                                t_waveforms = t1 - t0

                                self.autotree_pmt.fillPMTVariables(slow_waveform) 
                                self.autotree_pmt.fillTimePMTVariables(t_waveforms)
                                self.outTree_pmt.fill()

                                snr_ratio = slow_waveform.getSignalToNoise()
                                slow_wf_weights_snr[ch-1] = snr_ratio
                                sing_weig_avg_slow_wf[ch-1] = waveform_s[indx]

                                if ch == 4:

                                    slow_wf_weights_snr = [ (x / max(slow_wf_weights_snr)) for x in slow_wf_weights_snr ]   # Normalization of the weights

                                    for k in range(chs_to_analyse):

                                        for j in range(slow_sampling):

                                            weight_average_wf[j] += sing_weig_avg_slow_wf[k][j]*slow_wf_weights_snr[k]/sum(slow_wf_weights_snr)

                                    waveform_info_slow_wei_avg = { 'run' : run, 'event': event, 'channel' : 9, 'trigger' : trg, 'GE' : 9, 'sampling' : "slow"}

                                    t0 = time.perf_counter()
                                    slow_waveform_wei_avg = PMTreco(waveform_info_slow_wei_avg, weight_average_wf, self.pmt_params)
                                    t1 = time.perf_counter()
                                    t_waveforms = t1 - t0 

                                    self.autotree_pmt_avg.fillPMTVariables_average(slow_waveform_wei_avg)
                                    self.autotree_pmt_avg.fillTimePMTVariables_average(t_waveforms)
                                    #slow_waveform_wei_avg.__repr__()           ## Verbose of averaged waveform
                                    self.outTree_pmt_avg.fill()

                                    del waveform_info_slow_wei_avg
                                    del slow_waveform_wei_avg

                                del waveform_info
                                del slow_waveform

                ## If single sensor analysis, update event number each loop 
                if (self.options.camera_mode and not self.options.pmt_mode) or (self.options.pmt_mode and not self.options.camera_mode):
                    if camera == True or pmt == True:           ##this is needed otherwise it increases the event number with each bank
                        numev += 1  
                        self.outTree.fill()
                ## If both sensor analysis, update event number only after camera analysis since it comes after the PMT 
                elif (self.options.camera_mode and self.options.pmt_mode):
                    if camera == True:
                        self.outTree.fill()
                        numev += 1

                if camera==True:
                    del obj
                
                if pmt==True:
                    del waveform_f, waveform_s, header

        gc.collect()
                    
        ROOT.gErrorIgnoreLevel = savErrorLevel
                
if __name__ == '__main__':
    from optparse import OptionParser
    
    parser = OptionParser(usage='%prog h5file1,...,h5fileN [opts] ')
    parser.add_option('-r', '--run', dest='run', default='00000', type='string', help='run number with 5 characteres')
    parser.add_option('-j', '--jobs', dest='jobs', default=1, type='int', help='Jobs to be run in parallel (-1 uses all the cores available)')
    parser.add_option(      '--max-entries', dest='maxEntries', default=-1, type='int', help='Process only the first n entries')
    parser.add_option(      '--first-event', dest='firstEvent', default=-1, type='int', help='Skip all the events before this one')
    parser.add_option(      '--pdir', dest='plotDir', default='./', type='string', help='Directory where to put the plots')
    parser.add_option('-t',  '--tmp',  dest='tmpdir', default=None, type='string', help='Directory where to put the input file. If none is given, /tmp/<user> is used')
    parser.add_option(      '--max-hours', dest='maxHours', default=-1, type='float', help='Kill a subprocess if hanging for more than given number of hours.')
    parser.add_option('-o', '--outname', dest='outname', default='reco', type='string', help='prefix for the output file name')
    parser.add_option('-d', '--outdir', dest='outdir', default='./', type='string', help='Directory where to save the output file')
    
    (options, args) = parser.parse_args()
    
    f = open(args[0], "r")
    params = eval(f.read())

    for k,v in params.items():
        setattr(options,k,v)

    run = int(options.run)
    
    if options.debug_mode == 1:
        setattr(options,'outFile','%s_run%d_%s_debug.root' % (options.outname, run, options.tip))
        #if options.ev: options.maxEntries = options.ev + 1
        #if options.daq == 'midas': options.ev +=0.5 
    else:
        setattr(options,'outFile','%s_run%05d_%s.root' % (options.outname, run, options.tip))
    
    patt = re.compile('\S+_(\S+).txt')
    m = patt.match(args[0])
    detector = m.group(1)
    if run > 16798 :
        utilities.setPedestalRun_v2(options,detector)
    else:
        utilities.setPedestalRun(options,detector)
        
        
    try:
        USER = os.environ['USER']
        flag_env = 0
    except:
        flag_env = 1
        USER = os.environ['JUPYTERHUB_USER']
    #tmpdir = '/mnt/ssdcache/' if os.path.exists('/mnt/ssdcache/') else '/tmp/'
    # it seems that ssdcache it is only mounted on cygno-login, not in the batch queues (neither in cygno-custom)
    tmpdir = '/tmp/'
    os.system('mkdir -p {tmpdir}/{user}'.format(tmpdir=tmpdir,user=USER))
    tmpdir = '{tmpdir}/{user}/'.format(tmpdir=tmpdir,user=USER) if not options.tmpdir else options.tmpdir+"/"
    if sw.checkfiletmp(int(options.run),options.rawdata_tier,tmpdir):
        if options.rawdata_tier=='root':
            prefix = 'histograms_Run'
            postfix = 'root'
        elif options.rawdata_tier=='h5':
            prefix = 'histograms_Run'
            postfix = 'h5'
        else:
            prefix = 'run'
            postfix = 'mid.gz'
        options.tmpname = "%s/%s%05d.%s" % (tmpdir,prefix,int(options.run),postfix)
    else:
        if options.rawdata_tier == 'root':
            print ('Downloading file: ' + sw.swift_root_file(options.tag, int(options.run)))
            options.tmpname = sw.swift_download_root_file(sw.swift_root_file(options.tag, int(options.run)),int(options.run),tmpdir)
        else:
            print ('Downloading MIDAS.gz file for run ' + options.run)
    # in case of MIDAS, download function checks the existence and in case it is absent, dowloads it. If present, opens it
    if options.rawdata_tier == 'midas':
        ## need to open it (and create the midas object) in the function, otherwise the async run when multithreaded will confuse events in the two threads
        options.tmpname = [int(options.run),tmpdir,options.tag]
    if options.justPedestal:
        ana = analysis(options)
        print("Pedestals done. Exiting.")
        if options.donotremove == False:
            sw.swift_rm_root_file(options.tmpname)
        sys.exit(0)

    ana = analysis(options)
    nev = ana.getNEvents(options)
    print("This run has ",nev," events.")
    print("Will save plots to ",options.plotDir)
    os.system('cp utils/index.php {od}'.format(od=options.plotDir))
    
    nThreads = 1
    if options.jobs==-1:
        import multiprocessing
        nThreads = multiprocessing.cpu_count()
    else:
        nThreads = options.jobs

    t1 = time.perf_counter()
    firstEvent = 0 if options.firstEvent<0 else options.firstEvent
    lastEvent = nev if options.maxEntries==-1 else min(nev,firstEvent+options.maxEntries)
    print ("Analyzing from event %d to event %d" %(firstEvent,lastEvent))
    base = options.outFile.split('.')[0]
    if nThreads>1:
        print ("RUNNING USING ",nThreads," THREADS.")
        nj = int(nev/nThreads) if options.maxEntries==-1 else max(int((lastEvent-firstEvent)/nThreads),1)
        chunks = [(ichunk,i,min(i+nj-1,nev)) for ichunk,i in enumerate(range(firstEvent,lastEvent,nj))]
        print("Chunks = ",chunks)
        with futures.ProcessPoolExecutor(nThreads) as executor:
            futures_list = [executor.submit(ana,c) for c in chunks]
            for future in futures.as_completed(futures_list):
                # retrieve the result. This is crucial, because result() does not exit until the process is completed.
                future.result()
        print("Now hadding the chunks...")
        if flag_env == 0:
            os.system('hadd -k -f {outdir}/{base}.root {outdir}/{base}_chunk*.root'.format(base=base, outdir=options.outdir))
        else:
            os.system('/usr/bin/hadd -k -f {outdir}/{base}.root {outdir}/{base}_chunk*.root'.format(base=base, outdir=options.outdir))
        os.system('rm {outdir}/{base}_chunk*.root'.format(base=base, outdir=options.outdir))
    else:
        evrange=(-1,firstEvent,lastEvent)
        ana(evrange)
    t2 = time.perf_counter()
    if options.debug_mode == 1:
        print(f'Reconstruction Code Took: {t2 - t1} seconds')

    # now add the git commit hash to track the version in the ROOT file
    tf = ROOT.TFile.Open("{outdir}/{base}.root".format(base=base, outdir=options.outdir),'update')
    githash = ROOT.TNamed("gitHash",str(utilities.get_git_revision_hash()).replace('\n',''))
    githash.Write()
    total_time = ROOT.TNamed("total_time", str(t2-t1))
    total_time.Write()
    tf.Close()
    
    if options.donotremove == False:
        sw.swift_rm_root_file(options.tmpname)
