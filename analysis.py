import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.SetBatch(True)

from skimage import measure
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter
from skimage import data
from skimage import img_as_float
from skimage.morphology import reconstruction

def zs(th2):
    ped = ROOT.TH1F('ped','',1000,0,20000)
    nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
    [ped.Fill(th2.GetBinContent(ix+1,iy+1)) for ix in xrange(nx) for iy in xrange(ny)]
    pedmean = ped.GetMean()
    pedrms = ped.GetRMS()

    print pedmean, "  ,  ",pedrms
    
    th2_zs = ROOT.TH2D(th2.GetName()+'_zs',th2.GetName()+'_zs',nx,0,nx,ny,0,ny)
    th2_zs.SetDirectory(None)
    for ix in xrange(nx):
        for iy in xrange(ny):
            z = max(th2.GetBinContent(ix+1,iy+1) - pedmean,0)
            if z>10*pedrms: th2_zs.SetBinContent(ix+1,iy+1,z)
            th2_zs.GetZaxis().SetRangeUser(0.1,10000)
    return th2_zs
    
def analyze(rfile):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    tf = ROOT.TFile.Open(rfile)
    c1 = ROOT.TCanvas('c1','',600,600)
    for e in tf.GetListOfKeys() :
        name=e.GetName()
        obj=e.ReadObj()
        if not obj.InheritsFrom('TH2'): continue
        print "Processing histogram: ",name
        obj.RebinX(5); obj.RebinY(5)
        h2zs = zs(obj)

        print "Analyzing its contours..."
        x_bins = h2zs.GetNbinsX()
        y_bins = h2zs.GetNbinsY()
        bins = np.zeros((x_bins,y_bins))
        for y_bin in xrange(y_bins): 
            for x_bin in xrange(x_bins): 
                z = h2zs.GetBinContent(x_bin + 1,y_bin + 1)
                if z>0:
                    bins[y_bin,x_bin] = h2zs.GetBinContent(x_bin + 1,y_bin + 1)
                
        h2zs.Draw('colz')
        for ext in ['png']: #,'pdf']:
            c1.SaveAs('{name}.{ext}'.format(name=name,ext=ext))



        import matplotlib.pyplot as plt
        from skimage import data, filters

        fig, ax = plt.subplots()

        image = bins
        edges = filters.sobel(image)

        low = 0.20
        high = 0.35

        lowt = (edges > low).astype(int)
        hight = (edges > high).astype(int)

        ax.imshow(lowt, cmap='magma')
        ax.set_title('Low threshold')
        ax.axis('off')

        plt.tight_layout()

        plt.show()



            
        # from skimage import feature

        # # Compute the Canny filter for two values of sigma
        # edges1 = feature.canny(bins)
        # edges2 = feature.canny(bins, sigma=5)

        # # display results
        # fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(8, 3),
        #                                     sharex=True, sharey=True)

        # ax1.imshow(bins, cmap=plt.cm.gray)
        # ax1.axis('off')
        # ax1.set_title('noisy image', fontsize=20)

        # ax2.imshow(edges1, cmap=plt.cm.gray)
        # ax2.axis('off')
        # ax2.set_title('Canny filter, $\sigma=1$', fontsize=20)
        
        # ax3.imshow(edges2, cmap=plt.cm.gray)
        # ax3.axis('off')
        # ax3.set_title('Canny filter, $\sigma=3$', fontsize=20)
        
        # fig.tight_layout()
        
        # plt.show()
            
        # from math import sqrt
        # from skimage.feature import blob_dog, blob_log, blob_doh
        # from skimage.color import rgb2gray

        # import matplotlib.pyplot as plt

        # image_gray = bins
        # blobs_doh = blob_doh(image_gray, min_sigma=1, max_sigma=10, threshold=100,overlap=0.5)

        # fig, ax = plt.subplots()
        # ax.imshow(image_gray, interpolation='nearest')
        
        # for blob in blobs_doh:
        #     y, x, r = blob
        #     c = plt.Circle((x, y), r, color='lime', linewidth=2, fill=False)
        #     ax.add_patch(c)
            
        # plt.tight_layout()
        # plt.show()


if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog h5file1,...,h5fileN [opts] ')
    (options, args) = parser.parse_args()

    inputf = args[0]
    analyze(inputf)
