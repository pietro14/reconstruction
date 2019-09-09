import ROOT
import matplotlib.pyplot as plt
import numpy as np
from root_numpy import hist2array

def showPic(filen,run,event):
     tf = ROOT.TFile(filen)
     th2name = "pic_run{run:05}_ev{event}".format(run=run,event=event)
     print "--> ",th2name
     th2 = tf.Get(th2name)
     arr = hist2array(th2)
     fig = plt.figure(figsize=(10, 10))
     plt.imshow(arr.T,cmap='viridis', vmin=85, vmax=110, origin='lower' )
     plt.show()
     
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog -f file -r run -e event ')
    parser.add_option('-f', '--file',  dest='filename', type='string', help='input filename (ROOT)')
    parser.add_option('-r', '--run',   dest='run',      type='int',    help='run number')
    parser.add_option('-e', '--event', dest='event',    type='int',    help='event number')
    (options, args) = parser.parse_args()

    showPic(options.filename,options.run,options.event)
    
