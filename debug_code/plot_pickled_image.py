#!/usr/bin/env python
import optparse
import matplotlib.pyplot as plt
import pickle as pl


if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog run event ', version='%prog 1.0')
    parser.add_option('-c', '--cmap' , type='string'       , default='gray_r'      , help='palette for 2D image')

    (options, args) = parser.parse_args()

    font = {'family': 'arial',
            'color':  'black',
            'weight': 'normal',
            'size': 24,
    }
    
    fig_handle = pl.load(open('pic_run0{run}_ev{ev}_oriIma.pkl'.format(run=args[0],ev=args[1]),'rb'))
    plt.set_cmap(options.cmap)
    plt.clim(vmin=-5,vmax=15)
    plt.xlabel('x (pixels)', font, labelpad=20)
    plt.xticks(fontsize=14)
    plt.ylabel('y (pixels)', font, labelpad=20)
    plt.yticks(fontsize=14)
    plt.title('Image after zero suppression', font, pad=40)
    #plt.show()

    plt.savefig('pic_run0{run}_ev{ev}_oriIma_paper.pdf'.format(run=args[0],ev=args[1]))
    plt.gcf().clear()
    plt.close('all')
