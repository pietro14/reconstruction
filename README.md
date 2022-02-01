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

`python3 reconstruction.py configFile.txt --pdir plots --max-entries X -jX`

- *configFile.txt* is the configuration file with all the settings.
- *pdir* is the directory where the plots will be saved.
- *max-entries* is the number of images you want to analyse.
- *j* is the number of cores you want to use.



## Running the analysis code on MC data (updated on December 2021):
Firstly, create pedestal map (example with pedestal run 4159):
1. set `justPedestal` to `True` in `configFile_MC.txt`
2. add this line to `pedestals/pedruns_MC.txt`: `(4159, 4159) :   4159,`
3. copy pedestal root file (`histograms_Run04159.root`) to `/tmp/<username>`  
4. create pedmap with: `python reconstruction.py configFile_MC.txt -r 4159`
 
Now, you should have the file `pedestals/pedmap_run4159_rebin1.root`

Secondly, reconstruct digitized images (example with 3 root input files):
1. set `justPedestal` to `False` and `debug_mode` to `0` in `configFile_MC.txt`
2. copy pedestal root file (`histograms_Run04159.root`) in `pedestals/`
3. rename the 3 root input files: `histograms_Run00001.root`, `histograms_Run00002.root` and `histograms_Run00003.root` 
4. create the directory `./input` and copy the 3 root input files in it
5. add to `configFile_MC.txt`: `'tmpdir'  : './input/',`
6. set right pedestal for run 1 to 3 in `pedestals/pedruns_MC.txt`: `(1,3) :   4159,`
7. run in batch: `python scripts/submit_batch.py --outdir out/ --cfg configFile_MC.txt $(pwd) '[1-3]'`




# Prerequisite to run:
- Python 2.X (Stable version) or Python 3.X
- Root 6.X

## Python libraries:
- cycler==0.10.0
- decorator==4.4.0
- imageio==2.6.1
- joblib==0.14.0
- kiwisolver==1.1.0
- matplotlib==3.1.1
- networkx==2.4
- numpy==1.17.3
- Pillow==6.2.0
- pyparsing==2.4.2
- python-dateutil==2.8.0
- PyWavelets==1.1.0
- root-numpy==4.8.0
- scikit-image==0.16.1
- scikit-learn==0.21.3
- scipy==1.3.1
- six==1.12.0

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
