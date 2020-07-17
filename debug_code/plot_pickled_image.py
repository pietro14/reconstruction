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
            '1st': '1st_3D',
            '2nd': '2nd_3D',
            'all': 'all_3D',
            'sc' : 'sc_3D',
            }
            
    
    fig_handle = pl.load(open('pic_run0{run}_ev{ev}_{step}.pkl'.format(run=args[0],ev=args[1],step=suff[options.step]),'rb'))
    plt.set_cmap(options.cmap)
    if options.step=='raw':
        plt.title('Image after zero suppression', font, pad=40)
        plt.xlabel('x (pixels)', font, labelpad=20)
        plt.ylabel('y (pixels)', font, labelpad=20)
        plt.clim(vmin=-5,vmax=10)
        if int(args[0])==2317 and int(args[1])==8: ## example of split track
            plt.clim(vmin=0,vmax=25)
            csize = 160
            plt.xlim(4*240,4*(240+csize))
            plt.ylim(4*70,4*(70+csize))
        elif int(args[0])==2097 and int(args[1])==317: # ambe 60/40 (6 keV NR candidate)
            plt.clim(vmin=0,vmax=40)
            csize = 100
            plt.xlim(1200,1200+csize)
            plt.ylim(880,880+csize)
        elif int(args[0])==2097 and int(args[1])==59: # ambe 60/40 (6 keV NR candidate)
            plt.clim(vmin=0,vmax=40)
            csize = 100
            plt.xlim(660,660+csize)
            plt.ylim(1020,1020+csize)

    else:
        plt.title('Rebinned image', font, pad=40)
        plt.xlabel('x (macro-pixels)', font, labelpad=20)
        plt.ylabel('y (macro-pixels)', font, labelpad=20)
        if int(args[0])==1843: ## iron
            plt.clim(vmin=98,vmax=120)
            csize = 60
            plt.xlim(240,240+csize)
            plt.ylim(170,170+csize)
        elif int(args[0])==2317 and int(args[1])==8: ## example of split track
            plt.clim(vmin=100,vmax=110)
            csize = 160
            plt.xlim(240,240+csize)
            plt.ylim(70,70+csize)
        elif int(args[0])==2156: # cosmics
            plt.clim(vmin=98,vmax=110)            
            
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    #plt.show()

    for ext in ['pdf','png']:
        plt.savefig('pic_run0{run}_ev{ev}_{step}_paper.{ext}'.format(run=args[0],ev=args[1],step=suff[options.step],ext=ext))
    plt.gcf().clear()
    plt.close('all')
