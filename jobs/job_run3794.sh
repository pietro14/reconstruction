#!/bin/bash
ulimit -c 0 -S
ulimit -c 0 -H
set -e
cd /media/giorgio/DATA/Giorgio/Documenti/Uni/Dottorato/Lavoro/Analisis_algorithm/reconstruction_simple_tris/
source scripts/activate_cygno_lngs.sh
python3.8 reconstruction.py configFile.txt -r 3794 -j 24 
