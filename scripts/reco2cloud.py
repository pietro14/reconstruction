import os, glob, sys, re, csv, subprocess

# this trick of a bash file to be sourced is a trick, because we need to use the old python3.6 to use the cygno_libs
jobstring = '''#!/bin/bash
ulimit -c 0 -S
ulimit -c 0 -H
set -e
echo "Setting python to older 3.6 to use the cygno_libs"
alias python=python3.6
unset PYTHONPATH
cd FILESDIR
echo moved to FILESDIR
echo Starting to upload...
'''

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage='%prog workdir runs [options] ')
    parser.add_option(        '--dry-run',  dest='dryRun',   action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option('-o',   '--outdir',   dest='outdir',   type="string", default=None, help='outdirectory');
    parser.add_option('-b',   '--bucket',   dest='bucket',   type="string", default="cygno-analysis", help='bucket of the cloud');
    parser.add_option('-t',   '--tag',      dest='tag',      type="string", default="test", help='tag of the reconstruction code (corresponds to CLOUD directory under <bucket>/RECO/<tag>');
    (options, args) = parser.parse_args()

    outdir = os.path.abspath(options.outdir)
    allfiles = os.listdir(outdir)
    files = [f for f in allfiles if re.match('reco_run\d+_3D.root',f)]

    shfname = "reco2cloud_tmp.sh"
    shellfile = open(shfname, 'w')
    header = jobstring
    header = header.replace("FILESDIR",str(outdir))
    shellfile.write(header)
    for f in files:
        shellfile.write('cygno_repo put {bucket} {thefile} -t RECO/{tag} -v\n'.format(bucket=options.bucket,thefile=f,tag=options.tag))
    shellfile.write("echo DONE upload.")
    shellfile.close()

    print ("Ready to upload %d ROOT files to cloud in bucket '%s' and tag '%s'" % (len(files),options.bucket,options.tag))
    if options.dryRun:
        print ("Written shell file ",shfname)
        print ("To upload for real, remove --dry-run option")
    else:
        os.system("bash %s" % shfname)
        #os.system("rm %s" % shfname)

    print ("DONE")
    
        

    

    
    
