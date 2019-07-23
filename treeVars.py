class AutoFillTreeProducer:
    def __init__(self,tree):
        self.outTree = tree

    def createPMTVariables(self):
        self.outTree.branch('pmt_integral', 'F')
        self.outTree.branch('pmt_tot', 'F')
        self.outTree.branch('pmt_amplitude', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_time', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_prominence', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_fwhm', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_hm', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_risetime', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_falltime', 'F', lenVar='nPeak')

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
        self.outTree.branch('cmos_integral', 'F')
        self.outTree.branch('cmos_mean', 'F')
        self.outTree.branch('cmos_rms', 'F')
        self.outTree.branch('track_nhits', 'F', lenVar='nTrack')
        self.outTree.branch('track_integral', 'F', lenVar='nTrack')
        self.outTree.branch('track_length', 'F', lenVar='nTrack')
        self.outTree.branch('track_width', 'F', lenVar='nTrack')
        self.outTree.branch('track_lfullrms', 'F', lenVar='nTrack')
        self.outTree.branch('track_tfullrms', 'F', lenVar='nTrack')
        self.outTree.branch('track_lp0amplitude', 'F', lenVar='nTrack')
        self.outTree.branch('track_lp0prominence', 'F', lenVar='nTrack')
        self.outTree.branch('track_lp0fwhm', 'F', lenVar='nTrack')
        self.outTree.branch('track_lp0mean', 'F', lenVar='nTrack')
        self.outTree.branch('track_tp0fwhm', 'F', lenVar='nTrack')
        self.outTree.branch('track_iteration', 'F', lenVar='nTrack')
        self.outTree.branch('track_xmax', 'F', lenVar='nTrack')
        self.outTree.branch('track_xmin', 'F', lenVar='nTrack')
        self.outTree.branch('track_ymax', 'F', lenVar='nTrack')
        self.outTree.branch('track_ymin', 'F', lenVar='nTrack')

    def fillCameraVariables(self,pic,clusters):
        self.outTree.fillBranch('cmos_integral',pic.Integral())
        self.outTree.fillBranch('cmos_mean',pic.GetMean())
        self.outTree.fillBranch('cmos_rms',pic.GetRMS())
        self.outTree.fillBranch('track_nhits', [cl.size() for cl in clusters])
        self.outTree.fillBranch('track_integral', [cl.integral() for cl in clusters])
        self.outTree.fillBranch('track_length', [cl.shapes['long_width'] for cl in clusters])
        self.outTree.fillBranch('track_width', [cl.shapes['lat_width'] for cl in clusters])
        self.outTree.fillBranch('track_lfullrms', [cl.shapes['long_fullrms'] for cl in clusters])
        self.outTree.fillBranch('track_tfullrms', [cl.shapes['lat_fullrms'] for cl in clusters])
        self.outTree.fillBranch('track_lp0amplitude', [cl.shapes['long_p0amplitude'] for cl in clusters])
        self.outTree.fillBranch('track_lp0prominence', [cl.shapes['long_p0prominence'] for cl in clusters])
        self.outTree.fillBranch('track_lp0fwhm', [cl.shapes['long_p0fwhm'] for cl in clusters])
        self.outTree.fillBranch('track_lp0mean', [cl.shapes['long_p0mean'] for cl in clusters])
        self.outTree.fillBranch('track_tp0fwhm', [cl.shapes['lat_p0fwhm'] for cl in clusters])
        self.outTree.fillBranch('track_iteration', [cl.iterations() for cl in clusters])
        self.outTree.fillBranch('track_xmax', [cl.getXmax() for cl in clusters])
        self.outTree.fillBranch('track_xmin', [cl.getXmin() for cl in clusters])
        self.outTree.fillBranch('track_ymax', [cl.getYmax() for cl in clusters])
        self.outTree.fillBranch('track_ymin', [cl.getYmin() for cl in clusters])




        
