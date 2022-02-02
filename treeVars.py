import numpy as np

class AutoFillTreeProducer:
    def __init__(self,tree,eventContent):
        self.outTree = tree
        self.saveKillerVars = False
        self.eventContent = eventContent
        
    def createPMTVariables(self):
        self.outTree.branch('pmt_integral', 'F', title="integral of the PMT waveform")
        self.outTree.branch('pmt_tot', 'F', title="time over threshold of the PMT waveform")
        self.outTree.branch('pmt_amplitude', 'F', lenVar='nPeak', title="amplitude of the main peak of the PMT waveform")
        self.outTree.branch('pmt_time', 'F', lenVar='nPeak', title="time of the raising edge of the PMT waveform")
        self.outTree.branch('pmt_prominence', 'F', lenVar='nPeak', title="amplitude of the main peak of the PMT waveform wrt the surrounding baseline")
        self.outTree.branch('pmt_fwhm', 'F', lenVar='nPeak', title="full width at half-maximum of the main peak of the PMT waveform")
        self.outTree.branch('pmt_hm', 'F', lenVar='nPeak', title="half-maximum of the PMT waveform")
        self.outTree.branch('pmt_risetime', 'F', lenVar='nPeak', title="length of the PMT waveform rise-time")
        self.outTree.branch('pmt_falltime', 'F', lenVar='nPeak', title="length of the PMT waveform fall-time")

    def fillPMTVariables(self,peakFinder,sampleSize):
        self.outTree.fillBranch('pmt_integral',peakFinder.getIntegral()*sampleSize)        
        self.outTree.fillBranch('pmt_tot',peakFinder.getTot())
        self.outTree.fillBranch('pmt_amplitude',peakFinder.getAmplitudes())
        self.outTree.fillBranch('pmt_time',peakFinder.getPeakTimes())
        self.outTree.fillBranch('pmt_prominence',peakFinder.getProminences())
        self.outTree.fillBranch('pmt_fwhm',peakFinder.getFWHMs())
        self.outTree.fillBranch('pmt_hm',peakFinder.getHMs())
        self.outTree.fillBranch('pmt_risetime',peakFinder.getTimes('rise'))
        self.outTree.fillBranch('pmt_falltime',peakFinder.getTimes('fall'))

    def createCameraVariables(self):
        self.outTree.branch('cmos_integral', 'F', title="integral counts of the full CMOS sensor")
        self.outTree.branch('cmos_mean',     'F', title="average counts of the full CMOS sensor")
        self.outTree.branch('cmos_rms',      'F', title="RMS of the counts of the full CMOS sensor")

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
        self.outTree.fillBranch('{name}_lp0amplitude'.format(name=name), [cl.shapes['long_p0amplitude'] for cl in clusters])
        self.outTree.fillBranch('{name}_lp0prominence'.format(name=name), [cl.shapes['long_p0prominence'] for cl in clusters])
        self.outTree.fillBranch('{name}_lp0fwhm'.format(name=name),   [cl.shapes['long_p0fwhm'] for cl in clusters])
        self.outTree.fillBranch('{name}_lp0mean'.format(name=name),   [cl.shapes['long_p0mean'] for cl in clusters])
        self.outTree.fillBranch('{name}_tp0fwhm'.format(name=name),   [cl.shapes['lat_p0fwhm'] for cl in clusters])
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
