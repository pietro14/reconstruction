#!/usr/bin/env python

def s3_root_file(tag, run):
    sel = rootlocation(tag,run)    
    
    BASE_URL  = "https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/cygnus/"
    file_root = (sel+'/histograms_Run%05d.root' % run)
    return BASE_URL+file_root

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

def s3_download_root_file(url,run):
    import ROOT
    import os
    from urllib.request import urlretrieve
    USER = os.environ['USER']
    tmpdir = '/mnt/ssdcache/' if os.path.exists('/mnt/ssdcache/') else '/tmp/'
    os.system('mkdir -p {tmpdir}/{user}'.format(tmpdir=tmpdir,user=USER))
    tmpname = ("%s/%s/histograms_Run%05d.root" % (tmpdir,USER,run))
    urlretrieve(url, tmpname, reporthook)
    return tmpname 

def rootlocation(tag,run):
    
    if tag == 'Data':
        if (run>=936) and (run<=1601):
            sel = 'Data/LTD/Data_Camera/ROOT'
        elif (run>=1632) and (run<=4000):
            sel = 'Data/LAB'
        else:
            print("WARNING: Data taken with another DAQ or not yet uploaded to the cloud")
            exit()
    
    elif tag == 'DataMango':
            sel= 'Data/MAN'
            
    elif tag == 'MC':
        sel = 'Simulation'
        print("WARNING: automatic download for Simulated data not implemented yet")
        exit()
 
    return sel

