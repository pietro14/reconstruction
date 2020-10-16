echo "Entering the CYGNO virtual environment software..."
export LD_LYBRARY_PATH="/nfs/cygno/software/python3.8/pyenv/lib"
echo "Setting up ROOT..."
export PYTHONPATH=$PYTHONPATH:/home/$USER/.local/lib/python3.8/site-packages/
source /nfs/cygno/software/root-v6-22-00-py36-build/bin/thisroot.sh 
alias python="python3.8"
echo "DONE."
