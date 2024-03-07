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
## Running the analysis code (in general):

`python3.X reconstruction.py configFile.txt --pdir plots --max-entries X -jX -r R -t /inputdir/path/`

- *configFile.txt* is the configuration file with all the settings.
- *pdir* is the directory where the plots will be saved.
- *max-entries* is the number of images you want to analyse.
- *j* is the number of cores you want to use.
- *R* is the run number 
- */inputdir/path/* is the path to the directory where to put the input root files (for MC)



## Running the analysis code on MC data (updated on June 2023):
Firstly, create pedestal map (example with pedestal run 4904 of LNGS):
1. set `justPedestal` to `True` in `configFile_LNGS.txt`
2. create pedmap with: `python3 reconstruction.py configFile_LNGS.txt -r 4904`
 
Now, you should have the file `pedestals/pedmap_run4904_rebin1.root`

Then, reconstruct digitized (simulated) images:
1. set `justPedestal` to `False` (and `debug_mode` to `0`) in `configFile_MC.txt`
2. if needed, create the file `pedestals/runlog_MC.csv` with the script `scripts/make_runlog_tmp.py`. If the pedestal run is 4904, the simulated runs must be named with higher run numbers (4905, 4906 ...)
3. move the simulated runs in the input directory `/path/to/inputdir/` and recostruct them with:

`python3 reconstruction.py configFile_MC.txt -r 4905 -t /path/to/inputdir/`



# Prerequisite to run:
- Python 3.X (X>6)
- Root 6.X

## Python libraries:
- cycler>=0.10.0
- numpy >= 1.20
- cython >= 3.0.2
- cygnolib (see [Install cygno-lib](https://github.com/CYGNUS-RD/cygno?tab=readme-ov-file#install-the-cygno-library) )
  (which requires oidc-agent==4.2.6, boto3sts and midas== 0.0.1 )
- scipy>=1.3.1
- root-numpy==4.8.0 (only for older versions of the code before winter23 branch and with python < python3.10)
- uproot==5.2.2
- h5py==3.10.0
- scikit_image==0.22.0
- scikit-learn>=0.21.3
- mahotas==1.4.13
 

(Beware that some dependent packages will be installed automatically as requirements for some of these packages like:
 - kiwisolver==1.1.0
 - matplotlib==3.1.1
 - networkx==2.4
 - pyparsing==2.4.2
 - python-dateutil==2.8.0
 - six==1.12.0


...so proceed in order)

# Small guide to PMT parameters

Read [here](https://github.com/CYGNUS-RD/reconstruction/blob/winter23/ReadMe_PMT.md)

## Example

**Download the code from github:**

`git clone git@github.com:CYGNUS-RD/analysis.git`
or
`git clone https://github.com/CYGNUS-RD/analysis.git`

`cd analysis`


**Get a file for a specific run taken with the DAQ (eg. run 2113):**

`wget https://swift.cloud.infn.it:8080/v1/AUTH_1e60fe39fba04701aa5ffc0b97871ed8/Cygnus/Data/LAB/histograms_Run02113.root`


**Change the run number in the config file (Line 34)**

**https://github.com/CYGNUS-RD/analysis/blob/fng_18/configFile.txt#L34**

`emacs -nw configFile.txt`


**Then run the code on all the events**

`python reconstruction.py configFile.txt`


**If your computer has X cores**

(check  on linux with `cat /proc/cpuinfo | awk '/^processor/{print $3}’`)


**You can speed up the processing by parallelizing it:**

`python reconstruction.py configFile.txt -j X`


**You can now look at the output ROOT file with a tree containing 1 event/image with:**

`root -l reco_run02113_3D.root`
