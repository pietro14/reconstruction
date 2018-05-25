import os,math
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

class analysis:

    def __init__(self,rfile,rebin,options):
        self.rebin = rebin
        self.rfile = rfile
        self.pedfile_name = '{base}_ped_rebin{rb}.root'.format(base=os.path.splitext(self.rfile)[0],rb=self.rebin)
        if not options.calcPedestals and not os.path.exists(self.pedfile_name):
            print "WARNING: pedestal file ",self.pedfile_name, " not existing. First calculate them..."
            self.calcPedestal()
        print "Pulling pedestals..."
        pedrf = ROOT.TFile.Open(self.pedfile_name)
        self.pedmap = pedrf.Get('pedmap').Clone()
        self.pedmap.SetDirectory(None)
        self.pedmean = pedrf.Get('pedmean').GetMean()
        self.pedrms = pedrf.Get('pedrms').GetMean()
        pedrf.Close()
        
    def zs(self,th2):
        nx = th2.GetNbinsX(); ny = th2.GetNbinsY();
        th2_zs = ROOT.TH2D(th2.GetName()+'_zs',th2.GetName()+'_zs',nx,0,nx,ny,0,ny)
        th2_zs.SetDirectory(None)
        for ix in xrange(1,nx+1):
            for iy in xrange(1,ny+1):
                if not self.isGoodChannel(ix,iy): continue
                ped = self.pedmap.GetBinContent(ix,iy)
                noise = self.pedmap.GetBinError(ix,iy)
                z = max(th2.GetBinContent(ix,iy)-ped,0)
                if z>10*noise: th2_zs.SetBinContent(ix,iy,z)
                #print "x,y,z=",ix," ",iy," ",z,"   3*noise = ",3*noise
        th2_zs.GetZaxis().SetRangeUser(0,1000)
        return th2_zs

    def calcPedestal(self,maxImages=-1):
        nx=ny=2048
        nx=int(nx/self.rebin); ny=int(ny/self.rebin); 
        pedfile = ROOT.TFile.Open(self.pedfile_name,'recreate')
        pedmap = ROOT.TProfile2D('pedmap','pedmap',nx,0,nx,ny,0,ny,'s')
        tf = ROOT.TFile.Open(self.rfile)
        for i,e in enumerate(tf.GetListOfKeys()):
            if maxImages>-1 and i==maxImages: break
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print "Processing histogram: ",name
            obj.RebinX(self.rebin); obj.RebinY(self.rebin); 
            for ix in xrange(nx):
                for iy in xrange(ny):
                    pedmap.Fill(ix,iy,obj.GetBinContent(ix,iy)/float(math.pow(self.rebin,2)))

        tf.Close()
        pedfile.cd()
        pedmap.Write()
        pedmean = ROOT.TH1D('pedmean','pedestal mean',500,97,103)
        pedrms = ROOT.TH1D('pedrms','pedestal RMS',500,0,5)
        for ix in xrange(nx):
            for iy in xrange(ny):
               pedmean.Fill(pedmap.GetBinContent(ix,iy)) 
               pedrms.Fill(pedmap.GetBinError(ix,iy)) 
        pedmean.Write()
        pedrms.Write()
        pedfile.Close()

    def isGoodChannel(self,ix,iy):
        pedval = self.pedmap.GetBinContent(ix,iy)
        pedrms = self.pedmap.GetBinError(ix,iy)
        if pedval > 110: return False
        if pedrms < 0.2: return False
        if pedrms > 5: return False
        return True

    def store_evolution_in(self,lst):
        """Returns a callback function to store the evolution of the level sets in
        the given list.
        """
        def _store(x):
            lst.append(np.copy(x))
        return _store
    
    def reconstruct(self):
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        tf = ROOT.TFile.Open(self.rfile)
        c1 = ROOT.TCanvas('c1','',600,600)
        for e in tf.GetListOfKeys() :
            name=e.GetName()
            obj=e.ReadObj()
            if not obj.InheritsFrom('TH2'): continue
            print "Processing histogram: ",name
            obj.RebinX(self.rebin); obj.RebinY(self.rebin)
            obj.Scale(1./float(math.pow(self.rebin,2)))
            h2zs = self.zs(obj)
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
     
     
     
            # import matplotlib.pyplot as plt
            # from skimage import data, io, filters
     
            # fig, ax = plt.subplots()
     
            # image = bins
            # edges = filters.sobel(image)

            # low = 1
     
            # lowt = (edges > low).astype(int)
     
            # ax.imshow(lowt)
            # ax.set_title('Low threshold')     
            # plt.tight_layout()
     
            # plt.show()
     
     
     
                
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
            # import matplotlib.pyplot as plt
     
            # image_gray = bins
            # blobs_doh = blob_doh(image_gray, min_sigma=2, max_sigma=100, threshold=.1,overlap=0.1)
     
            # fig, ax = plt.subplots()
            # ax.imshow(image_gray, interpolation='nearest')
            
            # for blob in blobs_doh:
            #     y, x, r = blob
            #     c = plt.Circle((x, y), r, color='lime', linewidth=2, fill=False)
            #     ax.add_patch(c)
              
            # plt.tight_layout()
            # plt.show()

            from morphsnakes import(morphological_chan_vese,
                                    morphological_geodesic_active_contour,
                                    inverse_gaussian_gradient,
                                    checkerboard_level_set)

            fig, axes = plt.subplots(1, 2, figsize=(8, 4))
            ax = axes.flatten()
            
            # Morphological GAC
            image = img_as_float(bins)
            gimage = inverse_gaussian_gradient(image)

            # Initial level set
            init_ls = np.zeros(image.shape, dtype=np.int8)
            init_ls[10:-10, 10:-10] = 1
            # List with intermediate results for plotting the evolution
            evolution = []
            callback = self.store_evolution_in(evolution)
            ls = morphological_geodesic_active_contour(gimage, 230, init_ls,
                                                       smoothing=1, balloon=-1,
                                                       threshold=0.69,
                                                       iter_callback=callback)


            ax[0].imshow(image, cmap="gray")
            ax[0].set_axis_off()
            ax[0].contour(ls, [0.5], colors='r')
            ax[0].set_title("Morphological GAC segmentation", fontsize=12)
            
            ax[1].imshow(ls, cmap="gray")
            ax[1].set_axis_off()
            contour = ax[1].contour(evolution[0], [0.5], colors='g')
            contour.collections[0].set_label("Iteration 0")
            contour = ax[1].contour(evolution[100], [0.5], colors='y')
            contour.collections[0].set_label("Iteration 100")
            contour = ax[1].contour(evolution[-1], [0.5], colors='r')
            contour.collections[0].set_label("Iteration 230")
            ax[1].legend(loc="upper right")
            title = "Morphological GAC evolution"
            ax[1].set_title(title, fontsize=12)
            
            fig.tight_layout()
            plt.show()
            
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage='%prog h5file1,...,h5fileN [opts] ')
    parser.add_option('-r', '--rebin', dest='rebin', default=1, type='int', help='Rebin factor (same in x and y)')
    parser.add_option('-p', '--pedestal', dest='calcPedestals', default=False, action='store_true', help='First calculate the pedestals')
    (options, args) = parser.parse_args()

    inputf = args[0]
    ana = analysis(inputf,10,options)
    if options.calcPedestals:
        ana.calcPedestal()
    ana.reconstruct()
