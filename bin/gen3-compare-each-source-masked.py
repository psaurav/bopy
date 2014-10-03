#!/usr/bin/python

import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
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
        vmin = sys.argv[10]
        vmax = sys.argv[11]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), \
               "<trgdir> <refdir> <trg> <ref> <meas> <compare> <root> <cutoff> <wavelength> [<savefilename>]"
        sys.exit(0)

    file1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, trg, wavelength, param)
    file2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, ref, wavelength, param)
    filex1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, trg, wavelength, "A")
    filex2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, ref, wavelength, "A")

    nrow = 11
    ncol = 19

    images = []
    mx = -np.inf
    #fig, ax = plt.subplots(nrow, ncol)
    #fig = plt.figure()
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
            fig = plt.figure()
            f, (ax1, ax2) = fig.subplots(1, 2)
            ax1.imshow(d, interpolation='nearest', vmin=vmin, vmax=vmax)
            ax1.set_cmap('BlueRed3')
            ax1.colorbar()
            ax2.hist(d, 50, (vmin, vmax))
            mean = np.mean(d)
            ax2.vxline(mean, color='r')
            ax2.title('mean = {0}'.format(mean))
            f.show()
            
            mx = np.max((np.abs(np.min(d)), np.abs(np.max(d)), mx))
            ax = plt.subplot(gs[r, c])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title('%d'%(i), size='x-small')
            im = ax.imshow(d, interpolation='nearest')
            images.append(im)

