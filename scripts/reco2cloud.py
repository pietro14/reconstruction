import os, glob, sys, re, csv, subprocess

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

    cmd = "cd {odir}; {cyput}; cd -" # the cd is needed, the cygno lib only checks the file name, the abspath confuses it
    cmds = []
    for f in files:
        cmds.append(cmd.format(odir=outdir,cyput='cygno_repo put {bucket} {thefile} -t RECO/{tag} -v'.format(bucket=options.bucket,thefile=f,tag=options.tag)))

    cmds = sorted(cmds)
    if options.dryRun:
        for c in cmds:
            print(c)
        print ("To upload for real, remove --dry-run option")
    else:
        for c in cmds:
            os.system(c)
    print ("Ready to upload %d ROOT files to cloud in bucket '%s' and tag '%s'" % (len(files),options.bucket,options.tag))

    
        

    

    
    
