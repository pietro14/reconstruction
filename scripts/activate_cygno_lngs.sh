echo "Entering the CYGNO virtual environment software..."
#source /users/d/dimarcoe/python/env/bin/activate
export PATH=/users/d/dimarcoe/python/pyenv/bin:$PATH
export LD_LYBRARY_PATH=/users/d/dimarcoe/python/pyenv/lib
echo "Setting up ROOT..."
export PYTHONPATH="/afs/lngs.infn.it/user/d/dimarcoe/python/env/lib/python3.8/site-packages:/afs/lngs.infn.it/user/d/dimarcoe/.local/lib/python3.8/site-packages"
source /nfs/cygno/software/root-v6-22-00-py36-build/bin/thisroot.sh 
echo "DONE."
