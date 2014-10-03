#!/usr/bin/python

import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
#import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import AxesGrid
import bopy as bp

cdict3 = {'red':  ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75,1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25,1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

plt.register_cmap(name='BlueRed3', data=cdict3)

if __name__ == '__main__':

    try:
        trgdir = sys.argv[1]
        refdir = sys.argv[2]
        trg = sys.argv[3]
        ref = sys.argv[4]
        param = sys.argv[5]
        compare = sys.argv[6]
        root = sys.argv[7]
        cutoff = float(sys.argv[8]) 
        wavelength = sys.argv[9]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), \
               "<trgdir> <refdir> <trg> <ref> <meas> <compare> <root> <cutoff> <wavelength>"
        sys.exit(0)

    file1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, trg, wavelength, param)
    file2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, ref, wavelength, param)
    filex1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, trg, wavelength, "A")
    filex2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, ref, wavelength, "A")

    nrow = 11
    ncol = 19

    mx = -np.inf
    images = []
    fig = plt.figure()
    grid = AxesGrid(fig, 111, # similar to subplot(131)
                    nrows_ncols = (nrow, ncol),
                    axes_pad = 0.0,
                    cbar_mode = 'single',
                    cbar_pad = 0.1,
                    #cbar_location = 'top',
                    share_all=True,
                    label_mode = 'L'
                    )

    #fig, ax = plt.subplots(nrow, ncol)
    #gs = gridspec.GridSpec(nrow, ncol+1)
    for r in range(nrow):
        for c in range(ncol):
            i = c + ncol*r + 1
            d1 = bp.io.read_frame(os.path.join(trgdir, file1%(i)))
            d2 = bp.io.read_frame(os.path.join(refdir, file2%(i)))
            dmask = bp.utils.mask2_maxcutoff(bp.io.read_frame(os.path.join(trgdir, filex1%(i))), 
                                     bp.io.read_frame(os.path.join(refdir, filex2%(i))), 
                                     cutoff, toprowsmask=14) 
            if compare == '-':
                d = d1 - d2
            elif compare == '/':
                d = np.log(d1/d2)
            d = ma.masked_array(d, mask=dmask)
            mx = np.max((np.abs(np.min(d)), np.abs(np.max(d)), mx))

            im = grid[i-1].imshow(d, interpolation='nearest')
            grid[i-1].get_xaxis().set_ticks([])
            grid[i-1].get_yaxis().set_ticks([])
            grid[i-1].set_ylabel(i)
            images.append(im)

    #plt.subplots_adjust(wspace=0.1, hspace=0.0)
    vmin = -mx
    vmax = mx
    print "vmin:", vmin
    print "vmax:", vmax
    for im in images:
        im.set_clim(vmin=vmin)
        im.set_clim(vmax=vmax)
        im.set_cmap('BlueRed3')
    #plt.colorbar(images[0], cax=plt.subplot(gs[:, ncol]))
    plt.colorbar(im, cax = grid.cbar_axes[0])
    plt.tight_layout()
    plt.show()
