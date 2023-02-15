#!/usr/bin/env python
from os import listdir,popen,access,F_OK,getcwd,system
from sys import argv

from optparse import OptionParser
parser = OptionParser()
parser.add_option('--targetString', default=None, help='String to match to include')
parser.add_option('--skipString', default=None, help='Strong to match to skip')
parser.add_option('--runRange', default=(-1,-1), nargs=2, type='int', help='Run range (min, max). If not given, try to merge everything')
parser.add_option('--runList', default=None, type="string", help='file with runlist to hadd. Need to be interpreted as python list')
parser.add_option('--outputname', default="merged", type="string", help='output file name')
(opts,args) = parser.parse_args()


targetstring = ""
if opts.targetString is not None: 
    targetstring = opts.targetString

skipstring = ""
if opts.skipString is not None: 
    skipstring = opts.skipString

def printAndExec(cmd):
    print(cmd)
    result = popen(cmd).read()
    print(result)

filelist = []
selectedruns = []

if opts.runList and access(opts.runList,F_OK):
    runlist = open(opts.runList)
    selectedruns = eval(runlist.read())
    print ("Runs to be hadded taken from runlist: ",opts.runList)
    print ("Selected runs = ",selectedruns)
    print

system("mkdir -p Chunks")
for fn in listdir("."):
    if fn.count(".root") and fn.count("reco_run") and fn.count(targetstring) and (not skipstring or not fn.count(skipstring)):
        run = int(fn.replace("reco_run","").replace("_3D.root",""))
        if opts.runRange[0]>0 and run<opts.runRange[0]: continue
        if opts.runRange[1]>0 and run>opts.runRange[1]: continue
        if run not in selectedruns: continue
        filelist.append(fn)

outputfn = opts.outputname
if opts.runRange[0]>0: outputfn += "_rmin%i" % opts.runRange[0]
if opts.runRange[1]>0: outputfn += "_rmax%i" % opts.runRange[1]
outputfn += ".root"

result = sorted(filelist)
print("===> List of ",len(result)," files to be hadded:")
print(result)
print("================================")

if access(outputfn,F_OK):
    print("skipping",outputfn)
    exit(0)
if len(result) > 8:
    filesperintermediate = int((len(result))**0.5)
    subres = []
    nextone = 0
    while nextone < len(result):
        subres.append(result[nextone:nextone+filesperintermediate])
        nextone += filesperintermediate
    mediumlist = []    
    for i in range(len(subres)):
        mediumfile = outputfn.replace(".root","intermediate%i.root"%i)
        mediumlist.append(mediumfile)
        if access(mediumfile,F_OK):
            print("skipping",mediumfile)
            continue
        cmd = "hadd %s %s" % (mediumfile," ".join([fnn for fnn in subres[i]]))
        printAndExec(cmd)
    cmd = "hadd %s %s" % (outputfn," ".join(mediumlist)) 
    printAndExec(cmd)
else:    
    cmd = "hadd %s %s" % (outputfn," ".join([fnn for fnn in result]))
    printAndExec(cmd)
print("Now moving all the chunks files in Chunks and remove the intermediate files...")
system("mv %s Chunks" % " ".join([fnn for fnn in result]))
system("rm *intermediate*root")

