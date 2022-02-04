#!/bin/env python
# USAGE: python3.8 scripts/submit_batch.py $PWD "[3792-3796]" --outdir cosmics_1stset_261120 --dry-run
# USAGE 2.0: python3.8 scripts/submit_batch.py $PWD "[4455-4456]" --outdir test_batch --mh 8 --nev 200 24 --dry-run
# reading run list from a filem: python3.8 scripts/submit_batch.py $PWD datasets/lime_April2021.txt --outdir test_batch --mh 8 --nev 200 24 --dry-run
# remove --dry-run to submit for real (otherwise only the scripts are created and commands are printed)

jobstring  = '''#!/bin/bash
ulimit -c 0 -S
ulimit -c 0 -H
set -e
cd CYGNOBASE
source scripts/activate_cygno_lngs.sh
RECOSTRING
'''

import os, sys, re

def prepare_jobpack(jobdir,logdir,workdir,cmd,ijob=0):
    job_file_name = jobdir+'/job{ij}_run{r}.sh'.format(ij=ijob,r=run)
    log_file_name = logdir+'/job{ij}_run{r}.log'.format(ij=ijob,r=run)
    tmp_file = open(job_file_name, 'w')
    tmp_filecont = jobstring

    tmp_filecont = tmp_filecont.replace('RECOSTRING',cmd)
    tmp_filecont = tmp_filecont.replace('CYGNOBASE',workdir+'/')
    tmp_file.write(tmp_filecont)
    tmp_file.close()
    sub_cmd = 'qsub -q {queue} -l {ssd}ncpus={nt},mem={ram}mb -d {dpath} -e localhost:{logf} -o localhost:{logf} -j oe {jobf}'.format(dpath=workdir,logf=log_file_name,jobf=job_file_name,nt=nThreads,ram=RAM,ssd=ssdcache_opt,queue=options.queue)
    return sub_cmd

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage='%prog workdir runs [options] ')
    parser.add_option(        '--dry-run',  dest='dryRun',   action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option(        '--outdir',   dest='outdir',   type="string", default=None, help='outdirectory');
    parser.add_option(        '--nthreads', dest='nthreads', type="string", default=None, help='number of threads / job. If not given, it will be decided based on the queue');
    parser.add_option('-q',   '--queue',    dest='queue',    type="string", default='cygno-custom', help='queue to be used for the jobs');
    parser.add_option('--mh'  '--max-hours',dest='maxHours', default=-1, type='float', help='Kill a subprocess if hanging for more than given number of hours.')
    parser.add_option('--nev' '--event-chunks',dest='eventChunks', default=[], nargs=2,  type='int', help='T C: Total number of events to process and events per job')
    parser.add_option('--cfg' '--config-file',dest='configFile', default="configFile_LNF.txt",  type='string', help='the config file to be run')
    (options, args) = parser.parse_args()

    if len(args)<2:
        parser.print_help()
        exit(1)

    abswpath  = os.path.abspath(args[0]) 
    if not os.path.isdir(abswpath):
        raise RuntimeError('ERROR: {p} is not a valid directory. This is the base dir where the jobs run'.format(p=abswpath))

    runr = args[1]
    runs = []
    if os.path.exists(runr):
        print("Opening the text file %s with run ranges to process "%runr)
        f = open(runr, "r")
        runs.append(exec(f.read()))
    # if a simple run range with the format [minr-maxr] is given
    else:
        p = re.compile('\[(\d+)-(\d+)\]')
        m = p.match(runr)
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
    print ("SUBMITTING JOBS ON RUN RANGE = ",runs)

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

    # typically for LIME images (~4MP) MEM ~1GB/core
    if options.nthreads:
        nThreads = int(options.nthreads)
        RAM = nThreads * 1000
    else:
        if options.queue=='cygno':
            nThreads = 7
            RAM = nThreads * 1250 # max is 9 GB
        elif options.queue=='cygno-custom':
            nThreads = 64
            RAM = nThreads * 1500 # max is 110 GB
        else:
            print("WARNING: queue ",options.queue," not foreseen. Using 4 threads and 1 GB/thread by default")
            nThreads = 4
            RAM = nThreads * 1000
    ssdcache_opt = 'nodes=1:disk10,' if options.queue=='cygno-custom' else ''
    # cygno-custom mounts /mnt/ssdcache, not the other queues. It seems that sometimes testing its presence works, sometimes not,
    # so when not in cygno-custom, force the usage to /tmp
    #tmpdir_opt = '' if  options.queue=='cygno-custom' else ' --tmp /tmp/'
    # IT SEEMS THAT /mnt/ssdcache is not mounted even in cygno-custom
    tmpdir_opt = ' --tmp /tmp/'
    maxtime_opt = '' if options.maxHours < 0 else '--max-hours {hr}'.format(hr=options.maxHours)
    commands = []
    for run in runs:
        if len(options.eventChunks)>0:
            (totEv,evPerJob) = options.eventChunks
            print ("Preparing jobs for run {r}. The task subdivides a total of {nT} events in chunks of {nJ} events per job.".format(r=run,nT=totEv,nJ=evPerJob))
            for ij,firstEvent in enumerate(range(0,totEv,evPerJob)):
                print ("Will submit job #{ij}, processing event range: [{fev}-{lev}]".format(ij=ij,fev=firstEvent,lev=min(firstEvent+evPerJob,totEv)))
                cmd = 'python3.8 reconstruction.py {cfg} -r {r} -o reco_job{ijob} --first-event {fev} --max-entries {me} -j {nt} {tmpopt} {maxtimeopt}'.format(r=run,nt=nThreads,tmpopt=tmpdir_opt,maxtimeopt=maxtime_opt,fev=firstEvent,me=evPerJob,ijob=ij,cfg=options.configFile)
                sub_cmd = prepare_jobpack(jobdir,logdir,abswpath,cmd,ij)
                commands.append(sub_cmd)
        else:
            cmd = 'python3.8 reconstruction.py {cfg} -r {r} -j {nt} {tmpopt} {maxtimeopt}'.format(r=run,nt=nThreads,tmpopt=tmpdir_opt,maxtimeopt=maxtime_opt,cfg=options.configFile)
            prepare_jobpack(jobdir,logdir,abswpath,cmd)
            sub_cmd = prepare_jobpack(jobdir,logdir,abswpath,cmd)
            commands.append(sub_cmd)

    if options.dryRun:
        for c in commands:
            print (c)
    else:
        for c in commands:
            os.system(c)

    print ("DONE")


        
