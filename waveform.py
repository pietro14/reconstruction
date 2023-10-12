#!/usr/bin/env python

import os,math,sys,ctypes
import numpy as np
import ROOT
import random
ROOT.gROOT.SetBatch(True)
from scipy.signal import find_peaks,peak_widths
import pandas as pd

########################################################################################################################

class PMTreco:

    ## Initializes the waveform object with its main properties and operations
    def __init__(self, wf_info, y_array, pmt_params):

        self.run        = wf_info['run']                if 'run' in wf_info else 0
        self.event      = wf_info['event']              if 'event' in wf_info else 0
        self.channel    = wf_info['channel']            if 'channel' in wf_info else 0
        self.trigger    = wf_info['trigger']            if 'trigger' in wf_info else 0
        self.insideGE   = wf_info['GE']                 if 'GE' in wf_info else 0
        self.digitizer  = wf_info['sampling']           if 'sampling' in wf_info else None
        self.TTT        = wf_info['TTT']                if 'TTT' in wf_info else 0
        self.plotname   = 'WF__' + self.digitizer + '_run_' + str(self.run) + '_ev_' + str(self.event) + '_tr_' + str(self.trigger) + '_ch_' + str(self.channel) 

        self.threshold  = pmt_params['threshold']       if 'threshold' in pmt_params else 0
        self.height_RMS = pmt_params['height_RMS']      if 'height_RMS' in pmt_params else 1
        self.minDist    = pmt_params['minPeakDistance'] if 'minPeakDistance' in pmt_params else 1
        self.prominence = pmt_params['prominence']      if 'prominence' in pmt_params else None
        self.fixed_prom = pmt_params['fixed_prom']      if 'fixed_prom' in pmt_params else False
        self.width      = pmt_params['width']           if 'width' in pmt_params else 1
        self.resample   = pmt_params['resample']        if 'resample' in pmt_params else 1
        self.plotpy     = pmt_params['plotpy']          if 'plotpy' in pmt_params else False
        self.wf_in_tree = pmt_params['wf_in_tree']      if 'wf_in_tree' in pmt_params else False
        self.digit_tag  = pmt_params['digit_tag']       if 'digit_tag' in pmt_params else None
        self.pmt_verb   = pmt_params['pmt_verb']        if 'pmt_verb' in pmt_params else 0
        


        if self.digitizer == "fast":
            self.freq       = 1                                                                                               
            # self.freq       = 0.75                                        ## For now we save x_array as the sample array, not converted to ns                        
            self.x_array    = np.linspace(0,(1024-1)/self.freq,1024)
            self.sampling   = 1024

        elif self.digitizer == "slow":
            self.freq       = 1
            self.x_array    = np.linspace(0,(4000-1)/self.freq,4000)
            self.sampling   = 4000

        self.y_array    = y_array

        self.x_array_original = self.x_array
        self.y_array_original = self.y_array

        self.baseline   = self.getBaseline()
        
        # if self.channel in [1,2,3,4]:
        if self.channel in [1,2,3,4,5,9,91,92,93,94]:                   ## For tests
            self.invert_and_center_WF(self.baseline)

            if self.digitizer == "fast":                                    ##  For now let's not average the slow digitizer since hasn't been tested.
                self.moving_average(window_size = self.resample)

        # ## Test normalization
        # if self.channel in [9,99,999,9999,99999]:
        #     self.y_array = [x / max(self.y_array) for x in self.y_array]
        
        self.findPeaks(thr = self.threshold, height = self.height_RMS, 
            mindist = self.minDist, prominence = self.prominence, 
            fixed_prom = self.fixed_prom, width = self.width)
    
        if self.plotpy == True:
            self.plot_and_save( pdir = './waveforms', xlabel = None, ylabel = None, save = True, plot = False)

    ## Displays on the terminal the waveform information
    def __repr__(self):

        style = self.getPMTVerbose()

        GE = ""
        if self.insideGE == 0:
            GE = "No"
        else:
            GE = "Yes"

        generic = "\nWaveform analysis:\n\
            Run: {run:.0f}\
            Event: {event:.0f}\
            Trigger: {trigger:.0f}\
            Channel: {channel:.0f}".format(run=self.run,event=self.event,trigger = self.trigger,  channel=self.channel) 

        specific = "\n\
            Inside waveform: {GE}\n\
            Digitizer: {dgtz}\n\
            RMS: {rms:.2f}\n\
            Tot integral: {ti:.2f}\n\
            N_peaks: {np:.0f}\n\
            Peaks pos: {ppp}".format(GE = GE, dgtz = self.digitizer,rms = self.getRMS(),ti= self.getTotalIntegral(), \
                np =  len(self.getPeaks()) , ppp = str( np.around(self.getPeaksPositions()[:],2) ))

        tot_info = "\n\
            SNR: {snr:.2f}\n\
            TOT_time: {tot:.2f}\n\
            TOT_area: {ar:.2f}".format(snr = self.getSignalToNoise(), tot = self.getTOT('time'), ar = self.getTOT('area'))

        if style == 0 : return 
        if style == 1 : return print(generic)
        if style == 2 : return print(generic + specific)
        if style == 3 : return print(generic + specific + tot_info)

    #################  Operations of the waveform  #################  

    ## Inverts and centers to zero the waveform
    def invert_and_center_WF(self, baseline):

        demo_y = list(self.y_array) 
        for i in range(len(demo_y)):

            demo_y[i] -= baseline 
            demo_y[i] *= (-1.)

        self.y_array = tuple(demo_y)
        self.y_array_original = self.y_array

    ## Applies a low-pass filter, cutting the very high frequencies. Can be tuned with 'window size'
    def moving_average(self, window_size = 7, drop_NaN = True):

        tmp = pd.Series(self.y_array)
        moving_average = (tmp.rolling(window=window_size).mean().to_numpy())
        
        if drop_NaN:
            moving_average = (pd.Series(moving_average).dropna()).to_numpy()
    
        zeros = (window_size - 1) * [0.]                            ## Moving average changes the size of the array, so I fill it until 1024 with 0 to avoid issues
        tmp_f = np.append(moving_average,zeros)       
        
        self.y_array = tuple(tmp_f)

    ## Find the peaks in the waveform and their properties
    def findPeaks(self, height = None, thr = None, mindist = None, prominence = None, fixed_prom = None, width = None):

        height_thr = self.getRMS() * height

        if fixed_prom:
            prominence = max(self.y_array)/4
        else:
            prominence = 0.01

        peaks, properties = find_peaks(self.y_array, height=height_thr, threshold = thr, distance = mindist,  prominence = prominence, width = width)
        self.peaks = peaks
        self.properties = properties
        self.widths_half = peak_widths(self.y_array, self.peaks, rel_height=0.5)        
        self.widths_full = peak_widths(self.y_array, self.peaks, rel_height=0.95)               ## 95% is used due to the large RMS of the waveforms. Maybe can be set to 1.00 after waveform corrections

    #################  Retrieval of values  ###########3######  
    
    ## Get time of threshold of a given waveform
    def getTOT(self, mod):

        threshold_tot = self.getRMS() * 3
        tot_limits = 2*[0]

        tot_time = 0
        tot_area = 0
        
        beginning = False
        ending = False

        begin_x = 0
        end_x = 0

        # Defines how many consectuive samples must be above (below) the threshold to start (end) the signal
        density_start = 10          
        density_finish = 10

        c_up = 0
        c_down = 0

        for i in range(len(self.y_array)):

            elem = self.y_array[i]

            if elem > threshold_tot: c_up += 1
            else: c_up = 0

            if c_up == density_start and beginning == False:

                beginning = True
                begin_x = self.x_array[i] - c_up

            if beginning == True:

                if elem < threshold_tot : c_down += 1
                else: c_down = 0

                if c_down == density_finish and ending == False:

                    ending = True
                    end_x = self.x_array[i] - c_down

        if beginning == False:
            begin_x = 0
            end_x = 0

        if beginning == True and ending == False:
            endind = True
            end_x = self.sampling

        tot_time = end_x - begin_x

        for k in range(int(tot_time)):

            tot_area += self.y_array[int(begin_x) + k]

        tot_limits[0] = begin_x
        tot_limits[1] = end_x

        if   mod == 'time': return tot_time
        elif mod == 'area': return tot_area
        elif mod == 'limits': return tot_limits
        elif mod == 'thr': return threshold_tot
            


    def getSignalToNoise(self):
        
        signal = self.getMaxAmpl()
        noise = self.getRMS()
        S_N_ratio = signal/noise

        return S_N_ratio

    ## Retrieves run number of waveform
    def getRun(self):
        return self.run

    ## Retrieves event number of waveform
    def getEvent(self):
        return self.event

    ## Retrieves trigger number of waveform
    def getTrigger(self):
        return self.trigger

    ## Example for fillling the tree all in the same way
    ## Check 'treeVars.py' for more details
    # def getTrigger2(self):
    #     a = [self.trigger] * len(self.getPeaks())
    #     return a

    ## Retrieves channel number of waveform
    def getChannel(self):
        return self.channel

    def getInGE(self):
        return self.insideGE

    def getWfSaveInfo(self):
        return self.wf_in_tree

    def getSampling(self):
        return self.sampling

    def getDigitizerTag(self):
        return self.digit_tag

    def getPMTVerbose(self):
        return self.pmt_verb

    ## Retrieves run number of waveform
    def getTTT(self):
        return self.TTT

    ## Retrieves baseline of waveform with 100 samples starting from sample 0. Can be tuned
    def getBaseline(self, n_offset = 0):

        if   self.digitizer == "fast": n_samples = 100
        elif self.digitizer == "slow": n_samples = 400

        bl = np.mean(self.y_array[n_offset:n_offset+n_samples])

        return bl

    ## Retrieves RMS of waveform with 100 samples starting from sample 0. Can be tuned
    def getRMS(self, n_offset = 0):

        if   self.digitizer == "fast": n_samples = 100
        elif self.digitizer == "slow": n_samples = 400

        rms = np.std(self.y_array[n_offset:n_offset+n_samples])

        return rms

    ## Used to assign a ID (1,2,3,...) to each peak to more easily identify simultaneous peaks in different channels
    def getPeakIdentifier(self):
        peak_ID = [(ID+1) for ID in range(len(self.peaks))]
        return peak_ID

    ## Allows to save the whole waveform in the tree. Not used for the moment
    def getFullwaveform(self,axis):
        if axis == 'x':
            return self.x_array_original
        if axis == 'y':
            return self.y_array_original
        
    ## Retrieves amplitudes (mV) of identified peaks    
    def getAmplitudes(self):
        return self.properties["peak_heights"]

    ## Retrieves peaks found. Different from peaks_positions
    def getPeaks(self):
        return self.peaks

    ## Retrieves the actual position of the peaks given our x_axis
    def getPeaksPositions(self):
        return self.x_array[self.peaks]

    ## Retrieves the full and half width of the peaks
    def getPeakWidths(self,height):
        if height == 'full':
            return self.widths_full[0]
        if height == 'half':
            return self.widths_half[0]

    ## Retrieves the y_value used for the determination of the widths. Useful for plots
    def getHeightPeakBoundaries(self,height):
        if height == 'full':
            return self.widths_full[1]
        if height == 'half':
            return self.widths_half[1]

    ## Retrieves the x_values used for the determination of the widths (left and right). Useful for plots 
    def getPeakBoundaries(self,side,height):
 
        if height == 'full':   
            if side=='left': 
                return [self.x_array[int(x)] for x in self.widths_full[2]]
            if side == 'right':
                return [self.x_array[int(x)] for x in self.widths_full[3]]

        if height == 'half':   
            if side=='left': 
                return [self.x_array[int(x)] for x in self.widths_half[2]]
            if side == 'right':
                return [self.x_array[int(x)] for x in self.widths_half[3]]

    ## Retrieves the maximum amplitude of a waveform. Could be useful to check the saturation or do a quick particle ID
    def getMaxAmpl(self):
        return max(self.y_array)

    ## Retrieves total integral of the waveform (sum all the points). Could be improved 
    def getTotalIntegral(self, begin=None,end=None):
        if begin is not None:
            return np.sum(self.y_array[begin:end])
        else:
            return np.sum(self.y_array)

    ## Converts value (integral tipically) into charge
    def voltageToCharge(self,vlt):
        charge = vlt * (4./3.) * (1./50.)
        return charge

    ## Plot the waveforms with the peaks  founds and respective widths. Saves them into a folder called 'waveforms'
    def plot_and_save(self,pdir='./waveforms',xlabel='Time (ns)',ylabel='amplitude (mV)', save = True, plot = False):
        import matplotlib.pyplot as plt
        # from matplotlib import rc
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        # rc('text', usetex=True)

        # plot data and the found peaks
        plt.plot(self.x_array, self.y_array, color = 'black')
        plt.plot(self.x_array[self.peaks], self.getAmplitudes(), "x", color = 'orange')
        
        # plot some properties
        plt.hlines(y=self.getHeightPeakBoundaries(height = 'full'), xmin=self.getPeakBoundaries(side = 'left',height = 'full'),
            xmax=self.getPeakBoundaries(side = 'right', height = 'full'), color = "C3")     
        plt.hlines(y=self.getHeightPeakBoundaries(height = 'half'), xmin=self.getPeakBoundaries(side = 'left',height = 'half'),
            xmax=self.getPeakBoundaries(side = 'right', height = 'half'), color = "C2")     

        ###testtt
        plt.hlines(y= self.getTOT('thr'), xmin=self.getTOT(mod = 'limits')[0], xmax=self.getTOT(mod = 'limits')[1], color = "pink")     

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if plot == True:
            plt.show()
        
        if save == True:
            # for ext in ['pdf','png']:                                 ## To be changed if one want to save also in '.png' format
            for ext in ['pdf']:
                plt.savefig('{pdir}/{name}.{ext}'.format(pdir=pdir,name=self.plotname,ext=ext))
            plt.gcf().clear()


    def getWaveformID(self,size = 'peaks'):
        
        IDvec = []
        IDstr = "{}{}{}".format(self.event,self.trigger,self.channel)
        IDnum = int(IDstr)
        
        if size == 'waveforms':
            length = [1]
        elif size == 'peaks':
            length = self.getPeaks()
        elif size == 'fullWF':
            length = self.x_array

        for i in range(len(length)):
            IDvec.append(IDnum)
        
        return IDvec

    ###### Possible missing functions
    ##  function for peak majority (mj2_peak)    
    ##  function for calculate integral of mj2_peak peaks    
    ##  function of SUB-substructure charge? will depend on the relative height choosen
    ##  function for if Saturated(self):


