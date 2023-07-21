#!/bin/bash

python3 condor_folder_creator.py
python3 pmt_reconstruction.py $1 -j$2 --max-entries $3 -r $4 --outdir ./output_files
#python3 pmt_reconstruction.py $1 -j$2 --max-entries $3 -r $4 --tmp $5 --outdir output_files/
#python3 pmt_reconstruction.py configFile_LNGS.txt -j1 --max-entries 20 -r 12170 --tmp ../data/ --outdir ./output_files/

##the numbers will then be the options for the reconstructions, where their values are put in the other file
