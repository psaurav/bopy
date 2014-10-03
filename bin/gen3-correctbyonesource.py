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
        src = sys.argv[10]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), \
               "<trgdir> <refdir> <trg> <ref> <meas> <compare> <root> <cutoff> <wavelength> <src> [<savefilename>]"
        sys.exit(0)
    try:
        save_file_name = sys.argv[11]
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

    i = int(src)
    d1 = bp.io.read_frame(os.path.join(trgdir, file1%(i)))
    d2 = bp.io.read_frame(os.path.join(refdir, file2%(i)))
    dmask = bp.utils.mask2_maxcutoff(bp.io.read_frame(os.path.join(trgdir, filex1%(i))), 
                                     bp.io.read_frame(os.path.join(refdir, filex2%(i))), 
                                     cutoff, toprowsmask=14) 
    if compare == '-':
                d = d1 - d2
    elif compare == '/':
                d = np.log(d1/d2)
    elif compare == 'p':
                d = d1/d2
                d = d - 1.
    d = ma.masked_array(d, mask=dmask)
    mx = np.max((np.abs(np.min(d)), np.abs(np.max(d)), mx))
    d = d.compressed() 
    dmedian = np.median(d)
    print dmedian
    plt.hist(d, 100)
    plt.axvline(dmedian)
    #ax = plt.subplot(gs[r, c])
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)
    #ax.set_title('%d'%(i), size='x-small')
    #im = ax.imshow(d, interpolation='nearest')
    #images.append(im)

    #plt.subplots_adjust(wspace=0.1, hspace=0.0)
    #vmin = -mx
    #vmax = mx
    #print "vmin:", vmin
    #print "vmax:", vmax
    #for im in images:
    #    im.set_clim(vmin=vmin)
    #    im.set_clim(vmax=vmax)
    #    im.set_cmap('BlueRed3')
    #plt.colorbar(images[0], cax=plt.subplot(gs[:, ncol]))
    if save_file_name:
        fig.savefig(save_file_name, bbox_inches='tight', dpi=150)
    else:
        plt.show()
