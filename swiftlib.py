#!/usr/bin/env python
import uproot
import midas.file_reader
import h5py
import cygno as cy
import os

def swift_root_file(tag, run):    
    BASE_URL = "https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/"
    if 'MC' in tag:
        tag,path=tag.split('$')
        bucket= 'cygno-sim'
        sel=path
    if tag=='LNGS':
        if run>138:
            bucket='cygno-data'
            sel=tag
        else:
            print("WARNING: Runs prior 138 not present on cloud storage cygno-data/LNGS/")
            exit()
    if tag=='LNF':
        if (run>5000 and run<6633) or (run>10037):
            bucket='cygno-data'
            if run>10037:
               sel=tag
            else:
               sel='LAB'
        else:
            print("WARNING: Runs Data prior 5000 or in the range [6633-10036] are not present on cloud storage cygno-data/LNF/ or cygno-data/LAB")
            exit()
    if tag=='MAN':
        bucket='cygno-data'
        sel=tag
    
    file_root = ('/histograms_Run%05d.root' % run)
    return BASE_URL+bucket+'/'+sel+file_root

    
def reporthook(blocknum, blocksize, totalsize):
    import sys
    readsofar = blocknum * blocksize
    if totalsize > 0:
        percent = readsofar * 1e2 / totalsize
        s = "\r%5.1f%% %*d / %d" % (
            percent, len(str(totalsize)), readsofar, totalsize)
        sys.stderr.write(s)
        if readsofar >= totalsize: # near the end
            sys.stderr.write("\n")
    else: # total size is unknown
        sys.stderr.write("read %d\n" % (readsofar,))

def swift_download_root_file(url,run,tmp=None,justName=False):
    import ROOT
    from urllib.request import urlretrieve
    try:
        USER = os.environ['USER']
    except:
        try:
          USER = os.environ['JUPYTERHUB_USER']
        except:
          USER = "autoreco"
    tmpdir = tmp if tmp else '/tmp/'
    if tmpdir == '/tmp/':
         os.system('mkdir -p {tmpdir}/{user}'.format(tmpdir=tmpdir,user=USER))
         tmpname = ("%s/%s/histograms_Run%05d.root" % (tmpdir,USER,run))
    else:
         tmpname = ("%s/histograms_Run%05d.root" % (tmpdir,run))
    if not justName:
        urlretrieve(url, tmpname, reporthook)
    return tmpname 

def swift_read_root_file(tmpname):
    f  = uproot.open(tmpname)
    return f

def swift_read_h5_file(tmpname):
    f  = h5py.File(tmpname, 'r')
    return f

def swift_rm_root_file(tmpname):
    import os
    os.remove(tmpname)
    print("tmp file removed")

def checkfiletmp(run,tier,tmp=None):
    import os.path
    try:
        USER = os.environ['USER']
    except:
        try:
          USER = os.environ['JUPYTERHUB_USER']
        except:
          USER = "autoreco"
    tmpdir = tmp if tmp else '/tmp/'
         
    
    if tmpdir=='/tmp/':
         tmpdir= '{tmpdir}/{user}'.format(tmpdir=tmpdir,user=USER)
     
    os.system('mkdir -p {tmpdir}'.format(tmpdir=tmpdir))
    
    if tier=='root':
       prefix = 'histograms_Run'
       postfix = 'root'
    elif tier=='h5':
       prefix = 'histograms_Run'
       postfix = 'h5'
    else:
       prefix = 'run'
       postfix = 'mid.gz'
    return os.path.isfile("%s/%s%05d.%s" % (tmpdir,prefix,run,postfix))

def swift_download_midas_file(run,tmpdir,tag='LNGS'):
    print("download or open midas file for run ",int(run))
    mfile = cy.open_mid(int(run), path=tmpdir, cloud=True, tag=tag, verbose=True)
    return mfile
    
def root_TH2_name(root_file):
    pic = []
    wfm = []
    for i,e in enumerate(root_file.GetListOfKeys()):
        che = e.GetName()
        if ('pic_run' in str(che)):
            pic.append(che)
        elif ('wfm_run' in str(che)):
            wfm.append(che)
    return pic, wfm

def swift_pedestal_file(run):
    pedrun = selectPedestal(run)    
    
    BASE_URL = "https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/cygnus/Pedestals/"
    file_root = ('pedmap_run%05d_rebin1.root' % pedrun)
    return BASE_URL+file_root

def selectPedestal(run):
    
    f = open('runvspedmap.txt', "r")
    params = eval(f.read())
    
    for k,v in params.items():
        setattr(options,k,v)
        
    options.pedavailable
    
 
    return sel

