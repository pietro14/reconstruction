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

    def fillCameraVariables(self,clusters):
        self.outTree.fillBranch('track_integral', [cl.integral() for cl in clusters])
        self.outTree.fillBranch('track_length', [cl.getSize('long') for cl in clusters])
        self.outTree.fillBranch('track_width', [cl.getSize('lat') for cl in clusters])
        
