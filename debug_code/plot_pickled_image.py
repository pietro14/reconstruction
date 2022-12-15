#!/usr/bin/env python
import optparse
import matplotlib.pyplot as plt
import pickle as pl


if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog run event ', version='%prog 1.0')
    parser.add_option('-c', '--cmap' , type='string'       , default='gray_r'      , help='palette for 2D image')
    parser.add_option('-s', '--step' , type='string'       , default='raw'         , help='step of the reconstruction: raw,1st,2nd,all,sc')

    (options, args) = parser.parse_args()

    font = {'family': 'arial',
            'color':  'black',
            'weight': 'normal',
            'size': 24,
    }

    suff = {'raw': 'oriIma',
            '0th': '0th',
            '1st': '1st_3D',
            '2nd': '2nd_3D',
            'all': 'all_3D',
            'sc' : 'sc_3D',
            }
            
    
    fig_handle = pl.load(open('CAM0_{step}.pkl'.format(step=suff[options.step]),'rb'))
    plt.set_cmap(options.cmap)
    if options.step=='raw':
        plt.title('Image after zero suppression', font, pad=40)
        plt.xlabel('x (pixels)', font, labelpad=20)
        plt.ylabel('y (pixels)', font, labelpad=20)
        plt.clim(vmin=-5,vmax=10)
        plt.clim(vmin=0,vmax=25)
        csize = 160
        plt.ylim(250,2000)
        #plt.xlim(1140,1275)
        #plt.ylim(1500,1640)


    else:
        plt.title('Rebinned image', font, pad=40)
        plt.xlabel('x (macro-pixels)', font, labelpad=20)
        plt.ylabel('y (macro-pixels)', font, labelpad=20)
        plt.clim(vmin=0,vmax=25)
        plt.ylim(62,500)
        plt.gca().invert_yaxis()
#        cmap = plt.cm.get_cmap("binary")
#        reversed_color_map = cmap.reversed()
        
        csize = 160
        #plt.xlim(1140,1275)
        #plt.ylim(1500,1640)

            
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    #plt.show()

    for ext in ['pdf','png']:
        plt.savefig('pic_{step}_paper.{ext}'.format(step=suff[options.step],ext=ext))
    plt.gcf().clear()
    plt.close('all')
