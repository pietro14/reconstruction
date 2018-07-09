class AutoFillTreeProducer:
    def __init__(self,tree):
        self.outTree = tree

    def createPMTVariables(self):
        self.outTree.branch('pmt_amplitude', 'F', lenVar='nPeaks')
        self.outTree.branch('pmt_time', 'F', lenVar='nPeaks')
        self.outTree.branch('pmt_prominence', 'F', lenVar='nPeaks')
        self.outTree.branch('pmt_fwhm', 'F', lenVar='nPeaks')
        self.outTree.branch('pmt_hm', 'F', lenVar='nPeaks')

    def fillPMTVariables(self,peakFinder):
        self.outTree.fillBranch('pmt_amplitude',peakFinder.getAmplitudes())
        self.outTree.fillBranch('pmt_time',peakFinder.getPeakTimes())
        self.outTree.fillBranch('pmt_prominence',peakFinder.getProminences())
        self.outTree.fillBranch('pmt_fwhm',peakFinder.getFWHMs())
        self.outTree.fillBranch('pmt_hm',peakFinder.getHMs())
