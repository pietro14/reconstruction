echo "Entering the CYGNO virtual environment software..."
export PATH="/nfs/cygno/software/python3.8/pyenv/bin:$PATH"
export LD_LYBRARY_PATH="/nfs/cygno/software/python3.8/pyenv/lib"
echo "Setting up ROOT..."
export PYTHONPATH="/nfs/cygno/software/python3.8/env/lib/python3.8/site-packages:/nfs/cygno/software/python3.8/env/lib/python3.8/site-packages-pipinstalled"
source /nfs/cygno/software/root-v6-22-00-py36-build/bin/thisroot.sh 
alias python="python3.8"
echo "DONE."
