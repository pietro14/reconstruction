# analysis code to reconstruct and analyze track/clusters from CYGNUS camera

## Checkout instructions:
git@github.com:CYGNUS-RD/analysis.git
git checkout tune_sel

## Convert the H5 files in ROOT files with TH2D
This could be avoided reading directly the H5 files, but for now it is like this...
See instructions in:

https://github.com/CYGNUS-RD/hdf2root

The usual name of the runs after the conversion is *histogram_Run00494.root*
or *histograms_Run01515.root* if it was already a root file.


# Updated HOW-TO-RUN
## Running the analysis code:

`python3 reconstruction.py configFile.txt --pdir plots --max-entries X -jX`

- *configFile.txt* is the configuration file with all the settings.
- *pdir* is the directory where the plots will be saved.
- *max-entries* is the number of images you want to analyse.
- *j* is the number of cores you want to use.


# Dependences
- Python 3.X
- Root 6.X
- root-numpy
- Numpy
- Matplotlib
and a few other common python libraries

