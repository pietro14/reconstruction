# analysis code to reconstruct and analyze track/clusters from CYGNUS camera

## Checkout instructions:
git@github.com:CYGNUS-RD/analysis.git
git checkout tune_sel

## Convert the H5 files in ROOT files with TH2D
This could be avoided reading directly the H5 files, but for now it is like this...
See instructions in:

https://github.com/CYGNUS-RD/hdf2root

## Run the code for track/cluster analysis:
1. First calculate the pedestals with 1-pixel width and store it for later (slow...)
`./analysis.py -r 1 --numPedEvents 100 --max-entries 0 runXXX.root`

2. The real analysis can be done with whatever rebin of the image on the fly:
`./analysis.py --pdir plots ~/cernbox/CYGNUS/Run742.root`


