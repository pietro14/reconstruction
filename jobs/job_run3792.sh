#!/bin/bash
ulimit -c 0 -S
ulimit -c 0 -H
set -e
cd /media/giorgio/DATA/Giorgio/Documenti/Uni/Dottorato/Lavoro/Analisis_algorithm/reconstruction_simple_tris/Size_analyzed/
source scripts/activate_cygno_lngs.sh
./../After_reco/davidesize.exe ../reco_runs/reco_run3792-3792.root
