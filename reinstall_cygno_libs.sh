## This reinstalls the latest cygno libs required for the PMT_reco

#!/bin/bash
python3 -m pip uninstall cygno
python3 -m pip install git+https://github.com/CYGNUS-RD/cygno.git -U 
