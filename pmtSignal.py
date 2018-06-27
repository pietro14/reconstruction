#!/usr/bin/env python

import ROOT,math,os
from cameraChannel import cameraGeometry

class PMTSignal:
    def __init__(self,tgraph,clusters,options):
        self.waveform = tgraph
        self.clusters = clusters
        self.options = options
        
    def plotNice(self):
        sig_width = 150 #ns
        sig_min = 6150 # at least at FNG with DAQ

        canv = ROOT.TCanvas("cfr","",600,600)
        canv.SetLeftMargin(0.20)
        canv.SetBottomMargin(0.15)
        self.waveform.Draw('AL')
        self.waveform.GetXaxis().SetRangeUser(sig_min,sig_min+sig_width)
        self.waveform.GetXaxis().SetTitle('Time (ns)')
        self.waveform.GetYaxis().SetTitle('Amplitude (mV)')

        maxwidth = max([cl.widths['long'] for cl in self.clusters]) # mm
        title = 'N clusters = {nclu}, max length = {maxl:.1f}mm'.format(nclu=len(self.clusters), maxl=maxwidth)
        self.waveform.SetTitle(title)

        for ext in ['png','pdf']:
            canv.SaveAs('{od}/{name}.{ext}'.format(od=self.options.plotDir,name=self.waveform.GetName(),ext=ext))
        
