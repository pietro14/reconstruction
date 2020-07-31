#!/bin/env python

PYTHON3_DIR = '/nfs/cygno/software/python3.8'
PYTHON3     = PYTHON3_DIR+'/pyenv/bin/python3.8'
PYTHON_PATH = '{pdir}/env/lib/python3.8/site-packages:{pdir}/env/lib/python3.8/site-packages-pipinstalled'.format(pdir=PYTHON3_DIR)

jobstring  = '''#!/bin/bash
ulimit -c 0 -S
ulimit -c 0 -H
set -e
export CYGNO_BASE="CYGNOBASE"
cd $CYGNO_BASE
echo "Entering the CYGNO virtual environment software..."
export PATH="{pdir}/pyenv/bin:$PATH"
export LD_LIBRARY_PATH="{pdir}/pyenv/lib:$LD_LIBRARY_PATH"
echo "Setting up ROOT..."
export PYTHONPATH="{ppath}"
. /nfs/cygno/software/root-v6-22-00-py36-build/bin/thisroot.sh 
alias python="python3.8"
echo "DONE."
RECOSTRING
'''.format(pdir=PYTHON3_DIR,ppath=PYTHON_PATH)

import os, sys, re

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage='%prog workdir runs [options] ')
    parser.add_option(        '--dry-run'       , dest='dryRun'        , action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option(        '--outdir', dest='outdir', type="string", default=None, help='outdirectory');
    (options, args) = parser.parse_args()

    if len(args)<2:
        parser.print_help()
        exit(1)

    abswpath  = os.path.abspath(args[0]) 
    if not os.path.isdir(abswpath):
        raise RuntimeError('ERROR: {p} is not a valid directory. This is the base dir where the jobs run'.format(p=abswpath))

    runr = args[1]
    p = re.compile('\[(\d+)-(\d+)\]')
    m = p.match(runr)
    runs = []
    if m:
        minr=m.group(1); maxr=m.group(2)
        if not str.isdigit(minr) or not str.isdigit(maxr):
            raise RuntimeError ("not a good formatted run range: {rr}. Aborting.".format(rr=runr))
        minr=int(minr); maxr=int(maxr)
        if minr>maxr:
            raise RuntimeError ("Warning: first run > end run: {rr}. Aborting.".format(rr=runr))
        runs = list(range(minr,maxr+1))
        print("Will run reconstruction on run range ",runs)
    else:
        raise RuntimeError ("not a good formatted run range string: {rr}. Aborting.".format(rr=runr))

    if not options.outdir:
        raise RuntimeError('ERROR: give at least an output directory. This is where the jobs configs and logs are stored')
    else:
        absopath  = os.path.abspath(options.outdir)
        if not os.path.isdir(absopath):
            print ('making a directory ',absopath,' and running in it')
            os.system('mkdir -p {od}'.format(od=absopath))

    jobdir = absopath+'/jobs/'
    if not os.path.isdir(jobdir):
        os.system('mkdir {od}'.format(od=jobdir))
    logdir = absopath+'/logs/'
    if not os.path.isdir(logdir):
        os.system('mkdir {od}'.format(od=logdir))

    commands = []
    for run in runs:
        job_file_name = jobdir+'/job_run{r}.sh'.format(r=run)
        log_file_name = logdir+'/job_run{r}.log'.format(r=run)
        tmp_file = open(job_file_name, 'w')

        tmp_filecont = jobstring
        cmd = '{py3} reconstruction.py configFileNeutrons7030.txt -r {r}'.format(py3=PYTHON3,r=run)
        tmp_filecont = tmp_filecont.replace('RECOSTRING',cmd)
        tmp_filecont = tmp_filecont.replace('CYGNOBASE',abswpath+'/')
        tmp_file.write(tmp_filecont)
        tmp_file.close()
        
        sub_cmd = 'qsub -q cygno -d {dpath} -e localhost:{logf} -o localhost:{logf} -v LD_LIBRARY_PATH="{pdir}/pyenv/lib:$LD_LIBRARY_PATH",PYTHONPATH="{ppath}",PATH="{pdir}/pyenv/bin:$PATH" -l mem=8000mb,ncpus=8 {jobf}'.format(dpath=abswpath,logf=log_file_name,jobf=job_file_name,pdir=PYTHON3_DIR,ppath=PYTHON_PATH)
        commands.append(sub_cmd)

    if options.dryRun:
        for c in commands:
            print (c)
    else:
        for c in commands:
            os.system(c)

    print ("DONE")


        
