#!/usr/bin/env python 
import os

jobstring = '''#!/bin/sh
source /ua9/user/dimarcoe/cygno/MY_PYTHON_ENV/bin/activate
cd WORKDIR
RECONSTRUCTIONSTRING
'''



if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog txtfileWithRuns [opts] ')
    parser.add_option("--dry-run", dest="dryRun",    action="store_true", default=False, help="Do not run the job, only print the command");
    parser.add_option("--configFile", dest="configFile", default="configFileNeutrons7030.txt", type='string', help="configuration file to be used for the reconstruction. Default: configFileNeutrons7030.txt")
    parser.add_option('--outdir', dest='outdir', type="string", default=None, help='outdirectory')

    (options, args) = parser.parse_args()

    f = open(args[0], "r")
    runs = eval(f.read())
    print "RUNNING ON THE FOLLOWING RUNS: ",runs," with config: ",options.configFile

    absopath  = os.path.abspath(options.outdir)
    if not options.outdir:
        raise RuntimeError, 'ERROR: give at least an output directory. there will be a HUGE number of jobs!'
    else:
        if not os.path.isdir(absopath):
            print 'making a directory and running in it'
            os.system('mkdir -p {od}'.format(od=absopath))

    jobdir = absopath+'/jobs/'
    if not os.path.isdir(jobdir):
        os.system('mkdir {od}'.format(od=jobdir))
    logdir = absopath+'/logs/'
    if not os.path.isdir(logdir):
        os.system('mkdir {od}'.format(od=logdir))
    errdir = absopath+'/errs/'
    if not os.path.isdir(errdir):
        os.system('mkdir {od}'.format(od=errdir))

    srcfiles = []
    for j,r in enumerate(runs):
        cmd = 'python reconstruction.py {config} -r {run} -j -1'.format(config=options.configFile,run=r)
        job_file_name = jobdir+'/job_{j}.sh'.format(j=j)
        log_file_name = logdir+'/log_{j}.log'.format(j=j)
        tmp_file = open(job_file_name, 'w')
        tmp_filecont = jobstring
        tmp_filecont = tmp_filecont.replace('WORKDIR', os.environ['PWD'])
        tmp_filecont = tmp_filecont.replace('RECONSTRUCTIONSTRING', cmd)
        tmp_file.write(tmp_filecont)
        tmp_file.close()
        os.system("chmod a+x "+job_file_name)
        bsub_cmd = 'bsub -o {logf} -J job{j} {src}'.format(logf=log_file_name,j=j,src=job_file_name)
        if options.dryRun:
            print bsub_cmd
        else:
            os.system(bsub_cmd)

