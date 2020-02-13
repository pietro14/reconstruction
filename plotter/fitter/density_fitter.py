### python specific import
import argparse
import os
import sys
import shutil

from libPython.classUtils import fitSample
import libPython.fitUtils as fitRoot


parser = argparse.ArgumentParser(description='simple fitter')
parser.add_argument('--doFit'      , action='store_true'  , help = 'fit sample (sample should be defined in settings.py)')
parser.add_argument('--doPlot'     , action='store_true'  , help = 'plotting')

args = parser.parse_args()

samplesDef = {
    'bkg_sample'  : fitSample('bkg_sample','/Users/emanuele/Work/cygnus/analysis/plotter/fitter/templates/nosource_density.root'),
    ## NB both bkg (no source) and signal (source) are in the same file, just with two different names
    'data_sample' : fitSample('data_sample','/Users/emanuele/Work/cygnus/analysis/plotter/fitter/templates/nosource_density.root')
}

setattr( samplesDef['bkg_sample'], 'hist', 'fe_density')
setattr( samplesDef['data_sample'], 'hist', 'ambe_density')

#############################################################
########## fitting params to tune fit by hand if necessary
#############################################################
fitParNomFit = [
    "meanS[19,17,20]","sigmaS[1,0.5,3.0]",
    ]

sampleToFit         = samplesDef['data_sample']
sampleToFit.pathBkg = samplesDef['bkg_sample'].path
sampleToFit.histBkg = samplesDef['bkg_sample'].hist

if  args.doFit:
    fitRoot.histFitterNominal( sampleToFit, fitParNomFit )

    args.doPlot = False

if   args.doPlot:
    plottingDir = "fitplots"
    if os.path.exists( plottingDir ):
            shutil.rmtree( plottingDir )
    os.makedirs( plottingDir )
    fitRoot.histPlotter( 'fitresults.root', plottingDir )
