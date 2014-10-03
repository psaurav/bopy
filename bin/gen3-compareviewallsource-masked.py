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
        toprowsmask = int(sys.argv[10])
        genmask = sys.argv[11] 
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), \
               "<trgdir> <refdir> <trg> <ref> <meas> <compare> <root> <cutoff> <wavelength> <toprowsmask> <genmask> [<savefilename>]"
        sys.exit(0)

    if genmask is not 'None':
        genmask = np.load(genmask)
        
    try:
        save_file_name = sys.argv[12]
    except:
        save_file_name = None

    file1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, trg, wavelength, param)
    file2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, ref, wavelength, param)
    filex1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, trg, wavelength, "A")
    filex2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, ref, wavelength, "A")

    nrow = 11
    ncol = 19

    images = []
    mx = -np.inf
    #fig, ax = plt.subplots(nrow, ncol)
    fig = plt.figure()
    gs = gridspec.GridSpec(nrow, ncol+1)
    for r in range(nrow):
        for c in range(ncol):
            i = c + ncol*r + 1
            d1 = bp.io.read_frame(os.path.join(trgdir, file1%(i)))
            d2 = bp.io.read_frame(os.path.join(refdir, file2%(i)))
            dmask = bp.utils.mask2_maxcutoff(bp.io.read_frame(os.path.join(trgdir, filex1%(i))), 
                                     bp.io.read_frame(os.path.join(refdir, filex2%(i))), 
                                     cutoff, toprowsmask=toprowsmask, genmask=genmask) 
            if compare == '-':
                d = d1 - d2
            elif compare == '/':
                d = np.log(d1/d2)
            elif compare == 'p':
                d = d1/d2
                d = d - 1.
            d = ma.masked_array(d, mask=dmask)
            mx = np.max((np.abs(np.min(d)), np.abs(np.max(d)), mx))
            #if compare == 'p':
            #    d = d + 1.
            ax = plt.subplot(gs[r, c])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title('%d'%(i), size='x-small')
            im = ax.imshow(d, interpolation='nearest')
            images.append(im)

            d1 = ma.masked_array(d1, mask=dmask)
            d2 = ma.masked_array(d2, mask=dmask)
            d1_dynamrange = np.max(d1) - np.min(d1)
            d2_dynamrange = np.max(d2) - np.min(d2)
            print param, wavelength, i, d1_dynamrange, d2_dynamrange

    #plt.subplots_adjust(wspace=0.1, hspace=0.0)
    vmin = -mx
    vmax = mx
    #if compare == 'p':
    #    vmin = vmin + 1.
    #    vmax = vmax + 1.
    print "vmin:", vmin
    print "vmax:", vmax
    for im in images:
        im.set_clim(vmin=vmin)
        im.set_clim(vmax=vmax)
        im.set_cmap('BlueRed3')
    plt.colorbar(images[0], cax=plt.subplot(gs[:, ncol]))
    if save_file_name:
        fig.savefig(save_file_name, bbox_inches='tight', dpi=150)
    else:
        plt.show()
