class AutoFillTreeProducer:
    def __init__(self,tree):
        self.outTree = tree

    def createPMTVariables(self):
        self.outTree.branch('pmt_amplitude', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_time', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_prominence', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_fwhm', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_hm', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_raisetime', 'F', lenVar='nPeak')
        self.outTree.branch('pmt_falltime', 'F', lenVar='nPeak')

    def fillPMTVariables(self,peakFinder):
        self.outTree.fillBranch('pmt_amplitude',peakFinder.getAmplitudes())
        self.outTree.fillBranch('pmt_time',peakFinder.getPeakTimes())
        self.outTree.fillBranch('pmt_prominence',peakFinder.getProminences())
        self.outTree.fillBranch('pmt_fwhm',peakFinder.getFWHMs())
        self.outTree.fillBranch('pmt_hm',peakFinder.getHMs())
        self.outTree.fillBranch('pmt_raisetime',peakFinder.getTimes('rise'))
        self.outTree.fillBranch('pmt_falltime',peakFinder.getTimes('fall'))

    def createCameraVariables(self):
        self.outTree.branch('track_integral', 'F', lenVar='nTrack')
        self.outTree.branch('track_length', 'F', lenVar='nTrack')
        self.outTree.branch('track_width', 'F', lenVar='nTrack')
        self.outTree.branch('track_p0amplitude', 'F', lenVar='nTrack')
        self.outTree.branch('track_p0prominence', 'F', lenVar='nTrack')
        self.outTree.branch('track_p0fwhm', 'F', lenVar='nTrack')
        self.outTree.branch('track_p0mean', 'F', lenVar='nTrack')

    def fillCameraVariables(self,clusters):
        self.outTree.fillBranch('track_integral', [cl.integral() for cl in clusters])
        self.outTree.fillBranch('track_length', [cl.shapes['long_width'] for cl in clusters])
        self.outTree.fillBranch('track_width', [cl.shapes['lat_width'] for cl in clusters])
        self.outTree.fillBranch('track_p0amplitude', [cl.shapes['long_p0amplitude'] for cl in clusters])
        self.outTree.fillBranch('track_p0prominence', [cl.shapes['long_p0prominence'] for cl in clusters])
        self.outTree.fillBranch('track_p0fwhm', [cl.shapes['long_p0fwhm'] for cl in clusters])
        self.outTree.fillBranch('track_p0mean', [cl.shapes['long_p0mean'] for cl in clusters])



        
