import numpy as np
from waveform import PMTreco
import ROOT


class AutoFillTreeProducer:
    def __init__(self,tree,eventContent):
        self.outTree = tree
        self.saveKillerVars = False
        self.eventContent = eventContent

    def createEnvVariables(self):
        self.outTree.branch('Lime_pressure', 'F', title="Lime pressure")
        self.outTree.branch('Atm_pressure', 'F', title="Atmosperic pressure")
        self.outTree.branch('Lime_temperature', 'F', title="Lime temperature")
        self.outTree.branch('Atm_temperature', 'F', title="Atmosheric temperature")
        self.outTree.branch('Humidity', 'F', title="Humidity")
        

    def fillEnvVariables(self, dslow):
        self.outTree.fillBranch('Lime_pressure', dslow.P0IIn5)        
        self.outTree.fillBranch('Atm_pressure', dslow.P0IIn3)
        self.outTree.fillBranch('Lime_temperature', dslow.P0IIn0)
        self.outTree.fillBranch('Atm_temperature', dslow.P1UIn1)
        self.outTree.fillBranch('Humidity', dslow.P1UIn5*0.0375*1000-37.7)
    """
    ## This code allows to create a single TTree with all the waveform information of a single picture clustered together.
    ## This means the number of entries is the same as for the camera TTree.
    ## Requires post analysis to desantagle the information. This because the multi-dimensional vector (dimension = ionization event) are rolled down into a single dimension vector.
    ## It is working, but not updated with the more recent variables added.
    ## Check examples of construction to add additional variables.
    def createPMTVariables_singleTree(self,pmt_params):

        ### Waveform variables -- Identifiable with 'pmt_wf_ID'
        self.outTree.branch('pmt_nWaveforms',       'I', title  = 'Number of waveforms per event')
        self.outTree.branch('pmt_wf_ID',            'I', lenVar = 'pmt_nWaveforms', title = 'Waveform identifier')
        self.outTree.branch('pmt_wf_run',           'I', lenVar = 'pmt_nWaveforms', title = 'Waveform run')
        self.outTree.branch('pmt_wf_event',         'I', lenVar = 'pmt_nWaveforms', title = 'Waveform event/picture')
        self.outTree.branch('pmt_wf_trigger',       'I', lenVar = 'pmt_nWaveforms', title = 'Waveform trigger')
        self.outTree.branch('pmt_wf_channel',       'I', lenVar = 'pmt_nWaveforms', title = 'Waveform channel')
        self.outTree.branch('pmt_wf_insideGE',      'I', lenVar = 'pmt_nWaveforms', title = 'Check if waveform inside Global exposure')

        self.outTree.branch('pmt_wf_baseline',      'F', lenVar = 'pmt_nWaveforms', title = 'Waveform baseline')
        self.outTree.branch('pmt_wf_RMS',           'F', lenVar = 'pmt_nWaveforms', title = 'Waveform RMS')
        self.outTree.branch('pmt_wf_tot_integral',  'F', lenVar = 'pmt_nWaveforms', title = 'Whole waveform total integrated amplitude')
        self.outTree.branch('pmt_wf_tot_charge',    'F', lenVar = 'pmt_nWaveforms', title = 'Whole waveform total integrated charge')
        self.outTree.branch('pmt_wf_max_ampl',      'F', lenVar = 'pmt_nWaveforms', title = 'Waveform max voltage')
        self.outTree.branch('pmt_wf_nPeaks',        'I', lenVar = 'pmt_nWaveforms', title = 'Waveform number of peaks')

        ## Peak variables -- Identifiable with 'pmt_wf_ID_peaks'
        self.outTree.branch('pmt_wf_ID_peaks',      'I', lenVar = 'pmt_peaks_in_pic',  title = 'Waveform identifier for peaks')
        self.outTree.branch('pmt_peak_Number',      'I', lenVar = 'pmt_peaks_in_pic',  title = 'Peaks numbers')
        self.outTree.branch('pmt_peak_Position',    'F', lenVar = 'pmt_peaks_in_pic',  title = 'Peaks positions')
        self.outTree.branch('pmt_peak_Height',      'F', lenVar = 'pmt_peaks_in_pic',  title = 'Peaks heights')
        self.outTree.branch('pmt_peak_HalfWidth',   'F', lenVar = 'pmt_peaks_in_pic',  title = 'Peaks half widths')
        self.outTree.branch('pmt_peak_FullWidth',   'F', lenVar = 'pmt_peaks_in_pic',  title = 'Peaks full widths')

        ## Full waveforms -- for analysis -- Identifiable with 'pmt_wf_ID_full'
        if pmt_params['wf_in_tree'] == True:
            self.outTree.branch('pmt_wf_ID_full',       'F', lenVar = 'pmt_full_sized_wf',      title = 'Full waveforms length ID')
            self.outTree.branch('pmt_fullWaveform_X',   'F', lenVar = 'pmt_full_sized_wf',      title = 'Full waveform for in depth analysis') # I know 1024 there is a string, but I don't know how to manually input a number
            self.outTree.branch('pmt_fullWaveform_Y',   'F', lenVar = 'pmt_full_sized_wf',      title = 'Full waveform for in depth analysis') # I know 1024 there is a string, but I don't know how to manually input a number

        ## Example if one wants to fill all the branches with the same 'lenVar'.
        ## The variables are filled N times, where N is the number of peaks in the waveform.
        ## Check "waveform.py" where example with 'pmt_wf_trigger' is given.
         
        # self.outTree.branch('pmt_wf_trigger',       'I', lenVar = 'pmt_peaks_in_pic, title = 'Waveform trigger')
 
    def fillPMTVariables_singleTree(self,wfs):

        ## Waveform variables
        self.outTree.fillBranch('pmt_nWaveforms',        len(wfs))
        self.outTree.fillBranch('pmt_wf_ID',             [wf_id for wf in wfs for wf_id in wf.getWaveformID('waveforms')])
        self.outTree.fillBranch('pmt_wf_run',            [wf.getRun()         for wf in wfs])
        self.outTree.fillBranch('pmt_wf_event',          [wf.getEvent()       for wf in wfs])
        self.outTree.fillBranch('pmt_wf_trigger',        [wf.getTrigger()     for wf in wfs])
        self.outTree.fillBranch('pmt_wf_channel',        [wf.getChannel()     for wf in wfs])
        self.outTree.fillBranch('pmt_wf_insideGE',       [wf.getInGE()        for wf in wfs])

        self.outTree.fillBranch('pmt_wf_baseline',       [wf.getBaseline()       for wf in wfs])
        self.outTree.fillBranch('pmt_wf_RMS',            [wf.getRMS()            for wf in wfs])
        self.outTree.fillBranch('pmt_wf_tot_integral',   [wf.getTotalIntegral()  for wf in wfs])     
        self.outTree.fillBranch('pmt_wf_tot_charge',     [wf.voltageToCharge(wf.getTotalIntegral()) for wf in wfs])
        self.outTree.fillBranch('pmt_wf_max_ampl',       [wf.getMaxAmpl()        for wf in wfs])
        self.outTree.fillBranch('pmt_wf_nPeaks',         [len(wf.getPeaks())     for wf in wfs])     

        ## Peak variables
        self.outTree.fillBranch('pmt_wf_ID_peaks',       [wf_id_p     for wf in wfs for wf_id_p     in wf.getWaveformID('peaks')])
        self.outTree.fillBranch('pmt_peak_Number',       [peakid      for wf in wfs for peakid      in wf.getPeakIdentifier()])
        self.outTree.fillBranch('pmt_peak_Position',     [peakpos     for wf in wfs for peakpos     in wf.getPeaksPositions()])
        self.outTree.fillBranch('pmt_peak_Height',       [peakhei     for wf in wfs for peakhei     in wf.getAmplitudes()])
        self.outTree.fillBranch('pmt_peak_HalfWidth',    [peakhalfwid for wf in wfs for peakhalfwid in wf.getPeakWidths('half')])
        self.outTree.fillBranch('pmt_peak_FullWidth',    [peakfullwid for wf in wfs for peakfullwid in wf.getPeakWidths('full')])

        # # If at some point is necessary to save the whole waveform
        if (wfs[0].getWfSaveInfo() == True):
            self.outTree.fillBranch('pmt_wf_ID_full',       [wf_id_full for wf in wfs for wf_id_full in wf.getWaveformID('fullWF')])
            self.outTree.fillBranch('pmt_fullWaveform_X',   [x_point    for wf in wfs for x_point    in wf.getFullwaveform('x')])
            self.outTree.fillBranch('pmt_fullWaveform_Y',   [y_point    for wf in wfs for y_point    in wf.getFullwaveform('y')])

        ## Example for fillling the TTree all in the same way. Check explanation above
        # self.outTree.fillBranch('pmt_wf_trigger',       [trg for wf in wfs for trg in wf.getTrigger2()])
    """

    ####################################################################################################################################################################################

    ''' ## This code is the one being currently used. Each entry is an individual waveform. This allows to do individual cuts in the features/variables of the single waveforms. '''
    def createPMTVariables_multipleTrees(self,pmt_params):

        ### Waveform variables -- Identifiable with 'pmt_wf_ID'
        # self.outTree.branch('pmt_nWaveforms',       'I', title  = 'Number of waveforms per event')
        self.outTree.branch('pmt_wf_ID',            'I', title = 'Waveform identifier')
        self.outTree.branch('pmt_wf_run',           'I', title = 'Waveform run')
        self.outTree.branch('pmt_wf_event',         'I', title = 'Waveform event/picture')
        self.outTree.branch('pmt_wf_trigger',       'I', title = 'Waveform trigger')
        self.outTree.branch('pmt_wf_channel',       'I', title = 'Waveform channel')
        self.outTree.branch('pmt_wf_insideGE',      'I', title = 'Check if waveform inside Global exposure')
        self.outTree.branch('pmt_wf_sampling',      'I', title = 'Waveform coming from fast or slow digitizer')
        self.outTree.branch('pmt_wf_TTT',           'F', title = 'Waveform/Trigger time of arrival')

        self.outTree.branch('pmt_wf_baseline',      'F', title = 'Waveform baseline')
        self.outTree.branch('pmt_wf_RMS',           'F', title = 'Waveform RMS')
        self.outTree.branch('pmt_wf_tot_integral',  'F', title = 'Whole waveform total integrated amplitude')
        self.outTree.branch('pmt_wf_tot_charge',    'F', title = 'Whole waveform total integrated charge')
        self.outTree.branch('pmt_wf_max_ampl',      'F', title = 'Waveform max voltage')
        self.outTree.branch('pmt_wf_nPeaks',        'I', title = 'Waveform number of peaks')

        ##testttt
        self.outTree.branch('test_tot_time',             'F', title = "Time over threshold")
        self.outTree.branch('test_tot_area',             'F', title = "Area over threshold")

        ## Peak variables -- Identifiable with 'pmt_wf_ID_peaks'
        # self.outTree.branch('pmt_wf_ID_peaks',      'I', lenVar = 'pmt_wf_nPeaks',  title = 'Waveform identifier for peaks')
        self.outTree.branch('pmt_peak_Number',      'I', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks numbers')
        self.outTree.branch('pmt_peak_Position',    'F', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks positions')
        self.outTree.branch('pmt_peak_Height',      'F', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks heights')
        self.outTree.branch('pmt_peak_HalfWidth',   'F', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks half widths')
        self.outTree.branch('pmt_peak_FullWidth',   'F', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks full widths')

        ## Full waveforms -- for analysis -- Identifiable with 'pmt_wf_ID_full'
        if pmt_params['wf_in_tree'] == True:
            self.outTree.branch('pmt_wf_ID_full',       'F', lenVar = 'pmt_full_sized_wf',      title = 'Full waveforms length ID')
            self.outTree.branch('pmt_fullWaveform_X',   'F', lenVar = 'pmt_full_sized_wf',      title = 'Full waveform for in depth analysis') # I know 1024 there is a string, but I don't know how to manually input a number
            self.outTree.branch('pmt_fullWaveform_Y',   'F', lenVar = 'pmt_full_sized_wf',      title = 'Full waveform for in depth analysis') # I know 1024 there is a string, but I don't know how to manually input a number

    ## Working version to fill the tree with singular waveforms. Problem encountered because we fill the tree will the camera variables too many times
    def fillPMTVariables_multipleTrees(self,wf):

        ## Event variables
        self.outTree.fillBranch('pmt_wf_ID',                wf.getWaveformID('waveforms')[0])
        self.outTree.fillBranch('pmt_wf_run',               wf.getRun())
        self.outTree.fillBranch('pmt_wf_event',             wf.getEvent())
        self.outTree.fillBranch('pmt_wf_trigger',           wf.getTrigger())
        self.outTree.fillBranch('pmt_wf_channel',           wf.getChannel())
        self.outTree.fillBranch('pmt_wf_insideGE',          wf.getInGE())
        self.outTree.fillBranch('pmt_wf_sampling',          wf.getSampling())
        self.outTree.fillBranch('pmt_wf_TTT',               wf.getTTT())

        ## Waveform variables
        self.outTree.fillBranch('pmt_wf_baseline',          wf.getBaseline())
        self.outTree.fillBranch('pmt_wf_RMS',               wf.getRMS())
        self.outTree.fillBranch('pmt_wf_tot_integral',      wf.getTotalIntegral())
        self.outTree.fillBranch('pmt_wf_tot_charge',        wf.voltageToCharge(wf.getTotalIntegral()))
        self.outTree.fillBranch('pmt_wf_max_ampl',          wf.getMaxAmpl())
        self.outTree.fillBranch('pmt_wf_nPeaks',            len(wf.getPeaks()))

        ## testttt
        self.outTree.fillBranch('test_tot_time',                 wf.getTOT('time'))
        self.outTree.fillBranch('test_tot_area',                 wf.getTOT('area'))

        # If at some point is necessary to save the whole waveform
        if (wf.getWfSaveInfo() == True):
            self.outTree.fillBranch('pmt_wf_ID_full',       [wf_id_p for wf_id_p in wf.getWaveformID('fullWF')])
            self.outTree.fillBranch('pmt_fullWaveform_X',   [x for x in wf.getFullwaveform('x')])
            self.outTree.fillBranch('pmt_fullWaveform_Y',   [y for y in wf.getFullwaveform('y')])

        ## Peak variables
        self.outTree.fillBranch('pmt_peak_Number',          [pid for pid in wf.getPeakIdentifier()])
        self.outTree.fillBranch('pmt_peak_Position',        [pp for pp in wf.getPeaksPositions()])
        self.outTree.fillBranch('pmt_peak_Height',          [ph for ph in wf.getAmplitudes()])
        self.outTree.fillBranch('pmt_peak_HalfWidth',       [phw for phw in wf.getPeakWidths('half')])
        self.outTree.fillBranch('pmt_peak_FullWidth',       [pfw for pfw in wf.getPeakWidths('full')])


    ####################################################################################################################################################################################

    ''' 
    ## This code is used to create a different tree which is saving the weighted average waveform.
    ## We kept only the revelant variables related to the average waveform to avoid confusions.
    ## Could be further useful in the future once we want to save the PHYSICS EVENT variables mixing automatically Camera and PMT reconstructions
    '''
    def createPMTVariables_average(self,pmt_params):

        self.outTree.branch('pmt_wf_run',           'I', title = 'Waveform run')
        self.outTree.branch('pmt_wf_event',         'I', title = 'Waveform event/picture')
        self.outTree.branch('pmt_wf_trigger',       'I', title = 'Waveform trigger')
        self.outTree.branch('pmt_wf_channel',       'I', title = 'Waveform channel')                          ## I will call the average channel 9
        self.outTree.branch('pmt_wf_insideGE',      'I', title = 'Check if waveform inside Global exposure')
        self.outTree.branch('pmt_wf_sampling',      'I', title = 'Waveform coming from fast or slow digitizer')
        self.outTree.branch('pmt_wf_TTT',           'F', title = 'Waveform/Trigger time of arrival')

        self.outTree.branch('pmt_wf_nPeaks',        'I', title = 'Waveform number of peaks')

        ##testttt
        self.outTree.branch('test_tot_time',             'F', title = "Time over threshold")
        self.outTree.branch('test_tot_area',             'F', title = "Area over threshold")

        ## Peak variables -- Identifiable with 'pmt_wf_ID_peaks'
        self.outTree.branch('pmt_peak_Number',      'I', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks numbers')
        self.outTree.branch('pmt_peak_Position',    'F', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks positions')
        self.outTree.branch('pmt_peak_Height',      'F', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks heights')
        self.outTree.branch('pmt_peak_HalfWidth',   'F', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks half widths')
        self.outTree.branch('pmt_peak_FullWidth',   'F', lenVar = 'pmt_wf_nPeaks',  title = 'Peaks full widths')

    def fillPMTVariables_average(self,wf):

        ## Event variables
        self.outTree.fillBranch('pmt_wf_run',               wf.getRun())
        self.outTree.fillBranch('pmt_wf_event',             wf.getEvent())
        self.outTree.fillBranch('pmt_wf_trigger',           wf.getTrigger())
        self.outTree.fillBranch('pmt_wf_channel',           wf.getChannel())
        self.outTree.fillBranch('pmt_wf_insideGE',          wf.getInGE())
        self.outTree.fillBranch('pmt_wf_sampling',          wf.getSampling())
        self.outTree.fillBranch('pmt_wf_TTT',               wf.getTTT())

        ## Waveform variables
        self.outTree.fillBranch('pmt_wf_nPeaks',            len(wf.getPeaks()))

        ## testttt
        self.outTree.fillBranch('test_tot_time',                 wf.getTOT('time'))
        self.outTree.fillBranch('test_tot_area',                 wf.getTOT('area'))

        ## Peak variables
        self.outTree.fillBranch('pmt_peak_Number',          [pid for pid in wf.getPeakIdentifier()])
        self.outTree.fillBranch('pmt_peak_Position',        [pp for pp in wf.getPeaksPositions()])
        self.outTree.fillBranch('pmt_peak_Height',          [ph for ph in wf.getAmplitudes()])
        self.outTree.fillBranch('pmt_peak_HalfWidth',       [phw for phw in wf.getPeakWidths('half')])
        self.outTree.fillBranch('pmt_peak_FullWidth',       [pfw for pfw in wf.getPeakWidths('full')])

    ####################################################################################################################################################################################

    def createCameraVariables(self):
        self.outTree.branch('cmos_integral', 'F', title="integral counts of the full CMOS sensor")
        self.outTree.branch('cmos_mean',     'F', title="average counts of the full CMOS sensor")
        self.outTree.branch('cmos_rms',      'F', title="RMS of the counts of the full CMOS sensor")

    def createTimeVariables(self):
        self.outTree.branch('t_DBSCAN', 'F', title="DBSCAN time")
        self.outTree.branch('t_variables',     'F', title="Variables time")
        self.outTree.branch('lp_len',     'F', title="# pixel")
        
        self.outTree.branch('t_pedsub',     'F', title="pedestal subtraction")
        self.outTree.branch('t_saturation',     'F', title="saturation correction mode")
        self.outTree.branch('t_zerosup',     'F', title="zero suppression")
        self.outTree.branch('t_xycut',     'F', title="xy acceptance cut")
        self.outTree.branch('t_rebin',     'F', title="rebinning")

        self.outTree.branch('t_medianfilter',     'F', title="median filter")
        self.outTree.branch('t_noisered',     'F', title="noise reductor")


        
    def createClusterVariables(self,name='track'):
        chars = list(name)
        start = chars[0]; rest = chars[1:]
        sizeStr = 'n'+start.upper()+''.join(rest)
        self.outTree.branch('{name}_size'.format(name=name),         'F', lenVar=sizeStr, title="number of pixels of the cluster, without zero-suppression")
        self.outTree.branch('{name}_nhits'.format(name=name),        'F', lenVar=sizeStr, title="number of pixels of the cluster above zero-suppression threshold")
        self.outTree.branch('{name}_integral'.format(name=name),     'F', lenVar=sizeStr, title="uncalibrated integral of counts of all the pixels in the cluster")
        self.outTree.branch('{name}_corrintegral'.format(name=name), 'F', lenVar=sizeStr, title="density-corrected integral of the cluster (LEMON-specific calibration)")
        self.outTree.branch('{name}_rms'.format(name=name),          'F', lenVar=sizeStr, title="RMS of counts of all the pixels in the cluster")
        # filled only for the supercluster
        if name=='sc':
            self.outTree.branch('{name}_energy'.format(name=name),  'F', lenVar=sizeStr, title="calibrated energy of the cluster in keV (LEMON-specific calibration)")
            self.outTree.branch('{name}_pathlength'.format(name=name),    'F', lenVar=sizeStr, title="curved length of the cluster (made with skeletonization)")
            if self.eventContent["scfullinfo"] == True:
                self.outTree.branch('{name}_redpixIdx'.format(name=name),   'F',  lenVar=sizeStr, title="index of the first pixel in the reduced pixel (redpix) collection belonging to the cluster")
                self.outTree.branch('redpix_ix',        'I', lenVar='nRedpix', title="x coordinate of the pixel") 
                self.outTree.branch('redpix_iy',        'I', lenVar='nRedpix', title="y coordinate of the pixel") 
                self.outTree.branch('redpix_iz',        'F', lenVar='nRedpix', title="number of counts of the pixel (after pedestal subtraction)") 
        self.outTree.branch('{name}_theta'.format(name=name),        'F', lenVar=sizeStr, title="polar angle inclination of the major-axis of the cluster")
        self.outTree.branch('{name}_length'.format(name=name),       'F', lenVar=sizeStr, title="length of the major axis of the cluster")
        self.outTree.branch('{name}_width'.format(name=name),        'F', lenVar=sizeStr, title="length of the minor axis of the cluster")
        self.outTree.branch('{name}_longrms'.format(name=name),      'F', lenVar=sizeStr, title="truncated RMS of the cluster along the major axis")
        self.outTree.branch('{name}_latrms'.format(name=name),       'F', lenVar=sizeStr, title="truncated RMS of the cluster along the minor axis")
        self.outTree.branch('{name}_lfullrms'.format(name=name),     'F', lenVar=sizeStr, title="full RMS of the cluster along the major axis")
        self.outTree.branch('{name}_tfullrms'.format(name=name),     'F', lenVar=sizeStr, title="full RMS of the cluster along the minor axis")
        self.outTree.branch('{name}_lp0amplitude'.format(name=name), 'F', lenVar=sizeStr, title="amplitude of the main peak of the longitudinal cluster profile")
        self.outTree.branch('{name}_lp0prominence'.format(name=name),'F', lenVar=sizeStr, title="prominence of the main peak wrt the local baseline along the longitudinal cluster profile")
        self.outTree.branch('{name}_lp0fwhm'.format(name=name),      'F', lenVar=sizeStr, title="full width at half-maximum of the main peak of the longitudinal cluster profile")
        self.outTree.branch('{name}_lp0mean'.format(name=name),      'F', lenVar=sizeStr, title="mean position wrt the start of the cluster of the main peak of the longitudinal cluster profile")
        self.outTree.branch('{name}_tp0fwhm'.format(name=name),      'F', lenVar=sizeStr, title="full width at half-maximum of the main peak of the transverse cluster profile")
        self.outTree.branch('{name}_xmean'.format(name=name),        'F', lenVar=sizeStr, title="x position of the cluster energy baricenter")
        self.outTree.branch('{name}_ymean'.format(name=name),        'F', lenVar=sizeStr, title="y position of the cluster energy baricenter")
        self.outTree.branch('{name}_xmax'.format(name=name),         'F', lenVar=sizeStr, title="x position of the rightmost pixel of the cluster")
        self.outTree.branch('{name}_xmin'.format(name=name),         'F', lenVar=sizeStr, title="x position of the leftmost pixel of the cluster")
        self.outTree.branch('{name}_ymax'.format(name=name),         'F', lenVar=sizeStr, title="y position of the topmost pixel of the cluster")
        self.outTree.branch('{name}_ymin'.format(name=name),         'F', lenVar=sizeStr, title="y position of the bottommost pixel of the cluster")
        self.outTree.branch('{name}_pearson'.format(name=name),      'F', lenVar=sizeStr, title="Pearson coefficient of the cluster")
        self.outTree.branch('{name}_tgaussamp'.format(name=name),    'F', lenVar=sizeStr, title="amplitude of the Gaussian transverse profile")
        self.outTree.branch('{name}_tgaussmean'.format(name=name),   'F', lenVar=sizeStr, title="mean position of the Gaussian transverse profile")
        self.outTree.branch('{name}_tgausssigma'.format(name=name),  'F', lenVar=sizeStr, title="standard deviation of the Gaussian transverse profile")
        self.outTree.branch('{name}_tchi2'.format(name=name),        'F', lenVar=sizeStr, title="chi-squared of the Gaussian fit to the transverse profile")
        self.outTree.branch('{name}_tstatus'.format(name=name),      'F', lenVar=sizeStr, title="status of the Gaussian fit to the transverse profile")
        self.outTree.branch('{name}_lgaussamp'.format(name=name),    'F', lenVar=sizeStr, title="amplitude of the Gaussian longitudinal profile")
        self.outTree.branch('{name}_lgaussmean'.format(name=name),   'F', lenVar=sizeStr, title="mean position of the Gaussian longitudinal profile")
        self.outTree.branch('{name}_lgausssigma'.format(name=name),  'F', lenVar=sizeStr, title="standard deviation of the Gaussian longitudinal profile")
        self.outTree.branch('{name}_lchi2'.format(name=name),        'F', lenVar=sizeStr, title="chi-squared of the Gaussian fit to the longitudinal profile")
        self.outTree.branch('{name}_lstatus'.format(name=name),      'F', lenVar=sizeStr, title="status of the Gaussian fit to the longitudinal profile")

    def addCosmicKillerVariables(self,name='track'):
        self.saveKillerVars = True
        if name=='sc':
            self.outTree.branch('{name}_mindist'.format(name=name),       'F', lenVar=sizeStr, title="minimal distance of the cluster to a long track")
            self.outTree.branch('{name}_nmatchweak'.format(name=name),    'F', lenVar=sizeStr, title="number of pixels of the cluster matching a long track extrapolation (tight selection)")
            self.outTree.branch('{name}_nmatchrobust'.format(name=name),  'F', lenVar=sizeStr, title="number of pixels of the cluster matching a long track extrapolation (loose selection)")


    def fillCameraVariables(self,pic):
        self.outTree.fillBranch('cmos_integral',np.sum(pic))
        self.outTree.fillBranch('cmos_mean',np.mean(pic))
        self.outTree.fillBranch('cmos_rms',np.std(pic))

    def fillTimeVariables(self, t_variables, t_DBSCAN, lp, t_pedsub, t_saturation, t_zerosup, t_xycut, t_rebin, t_medianfilter, t_noisered):
        self.outTree.fillBranch('t_DBSCAN', t_DBSCAN)
        self.outTree.fillBranch('t_variables',t_variables)
        self.outTree.fillBranch('lp_len',lp)

        self.outTree.fillBranch('t_pedsub', t_pedsub)
        self.outTree.fillBranch('t_saturation', t_saturation)
        self.outTree.fillBranch('t_zerosup', t_zerosup)
        self.outTree.fillBranch('t_xycut', t_xycut)
        self.outTree.fillBranch('t_rebin', t_rebin)

        self.outTree.fillBranch('t_medianfilter', t_medianfilter)
        self.outTree.fillBranch('t_noisered', t_noisered) 

    def fillClusterVariables(self,clusters,name='track'):
        chars = list(name)
        start = chars[0]; rest = chars[1:]
        sizeStr = 'n'+start.upper()+''.join(rest)
        self.outTree.fillBranch('{name}_size'.format(name=name),     [cl.size() for cl in clusters])
        self.outTree.fillBranch('{name}_nhits'.format(name=name),    [cl.sizeActive() for cl in clusters])
        self.outTree.fillBranch('{name}_integral'.format(name=name), [cl.integral() for cl in clusters])
        self.outTree.fillBranch('{name}_corrintegral'.format(name=name), [cl.corr_integral() for cl in clusters])
        self.outTree.fillBranch('{name}_rms'.format(name=name),      [cl.rms() for cl in clusters])
        # filled only for the supercluster
        if name=='sc':
            self.outTree.fillBranch('{name}_energy'.format(name=name), [cl.calibratedEnergy for cl in clusters])
            self.outTree.fillBranch('{name}_pathlength'.format(name=name),    [cl.pathlength for cl in clusters])
            if self.saveKillerVars == True:
                self.outTree.fillBranch('{name}_mindist'.format(name=name), [cl.minDistKiller for cl in clusters])
                self.outTree.fillBranch('{name}_nmatchweak'.format(name=name), [cl.nMatchKillerWeak for cl in clusters])
                self.outTree.fillBranch('{name}_nmatchrobust'.format(name=name), [cl.nMatchKiller for cl in clusters])
            if self.eventContent["scfullinfo"] == True:
                selection = self.eventContent["scpixels_sel"]
                redPixIdxs = []
                ix = []; iy = []; iz = []
                for cl in clusters:
                    if cl.shapes['long_width'] < float(selection["max_len"]) and cl.integral() > float(selection["min_integral"]):
                        redPixIdxs.append(len(ix))
                        ix = ix + [round(cl.xallpixelcoord[i]) for i in range(cl.nallintpixels)]
                        iy = iy + [round(cl.yallpixelcoord[i]) for i in range(cl.nallintpixels)]
                        iz = iz + [round(cl.zallpixel[i]*10)/10. for i in range(cl.nallintpixels)]
                    else:
                        redPixIdxs.append(-1)
                self.outTree.fillBranch('redpix_ix', ix)
                self.outTree.fillBranch('redpix_iy', iy)
                self.outTree.fillBranch('redpix_iz', iz)
                self.outTree.fillBranch('{name}_redpixIdx'.format(name=name),   redPixIdxs)
        self.outTree.fillBranch('{name}_theta'.format(name=name),    [cl.shapes['theta'] for cl in clusters])
        self.outTree.fillBranch('{name}_length'.format(name=name),   [cl.shapes['long_width'] for cl in clusters])
        self.outTree.fillBranch('{name}_width'.format(name=name),    [cl.shapes['lat_width'] for cl in clusters])
        self.outTree.fillBranch('{name}_longrms'.format(name=name),  [cl.shapes['longrms'] for cl in clusters])
        self.outTree.fillBranch('{name}_latrms'.format(name=name),   [cl.shapes['latrms'] for cl in clusters])
        self.outTree.fillBranch('{name}_lfullrms'.format(name=name), [cl.shapes['long_fullrms'] for cl in clusters])
        self.outTree.fillBranch('{name}_tfullrms'.format(name=name), [cl.shapes['lat_fullrms'] for cl in clusters])
        # self.outTree.fillBranch('{name}_lp0amplitude'.format(name=name), [cl.shapes['long_p0amplitude'] for cl in clusters])
        # self.outTree.fillBranch('{name}_lp0prominence'.format(name=name), [cl.shapes['long_p0prominence'] for cl in clusters])
        # self.outTree.fillBranch('{name}_lp0fwhm'.format(name=name),   [cl.shapes['long_p0fwhm'] for cl in clusters])
        # self.outTree.fillBranch('{name}_lp0mean'.format(name=name),   [cl.shapes['long_p0mean'] for cl in clusters])
        # self.outTree.fillBranch('{name}_tp0fwhm'.format(name=name),   [cl.shapes['lat_p0fwhm'] for cl in clusters])
        self.outTree.fillBranch('{name}_xmean'.format(name=name),     [cl.shapes['xmean'] for cl in clusters])
        self.outTree.fillBranch('{name}_ymean'.format(name=name),     [cl.shapes['ymean'] for cl in clusters])
        self.outTree.fillBranch('{name}_xmax'.format(name=name),      [cl.shapes['xmax'] for cl in clusters])
        self.outTree.fillBranch('{name}_xmin'.format(name=name),      [cl.shapes['xmin'] for cl in clusters])
        self.outTree.fillBranch('{name}_ymax'.format(name=name),      [cl.shapes['ymax'] for cl in clusters])
        self.outTree.fillBranch('{name}_ymin'.format(name=name),      [cl.shapes['ymin'] for cl in clusters])
        self.outTree.fillBranch('{name}_pearson'.format(name=name),   [cl.getPearson() for cl in clusters])
        self.outTree.fillBranch('{name}_tgaussamp'.format(name=name),   [cl.shapes['tgaussamp'] for cl in clusters])
        self.outTree.fillBranch('{name}_tgaussmean'.format(name=name),  [cl.shapes['tgaussmean'] for cl in clusters])
        self.outTree.fillBranch('{name}_tgausssigma'.format(name=name), [cl.shapes['tgausssigma'] for cl in clusters])
        self.outTree.fillBranch('{name}_tchi2'.format(name=name),       [cl.shapes['tchi2'] for cl in clusters]) 
        self.outTree.fillBranch('{name}_tstatus'.format(name=name),     [cl.shapes['tstatus'] for cl in clusters])
        self.outTree.fillBranch('{name}_lgaussamp'.format(name=name),   [cl.shapes['lgaussamp'] for cl in clusters])
        self.outTree.fillBranch('{name}_lgaussmean'.format(name=name),  [cl.shapes['lgaussmean'] for cl in clusters])
        self.outTree.fillBranch('{name}_lgausssigma'.format(name=name), [cl.shapes['lgausssigma'] for cl in clusters])
        self.outTree.fillBranch('{name}_lchi2'.format(name=name),       [cl.shapes['lchi2'] for cl in clusters]) 
        self.outTree.fillBranch('{name}_lstatus'.format(name=name),     [cl.shapes['lstatus'] for cl in clusters])
