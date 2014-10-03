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
        source = sys.argv[10]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), \
               "<trgdir> <refdir> <trg> <ref> <meas> <compare> <root> <cutoff> <wavelength> source"
        sys.exit(0)

    file1 = "{0}_{1}_wl{2}_s{3}_{4}.npy".format(root, trg, wavelength, source, param)
    file2 = "{0}_{1}_wl{2}_s{3}_{4}.npy".format(root, ref, wavelength, source, param)
    filex1 = "{0}_{1}_wl{2}_s{3}_{4}.npy".format(root, trg, wavelength, source, "A")
    filex2 = "{0}_{1}_wl{2}_s{3}_{4}.npy".format(root, ref, wavelength, source, "A")

    d1 = bp.io.read_frame(os.path.join(trgdir, file1))
    d2 = bp.io.read_frame(os.path.join(refdir, file2))

    if compare == '-':
        d = d1 - d2
    elif compare == '/':
        d = np.log(d1/d2)

    dmask = bp.utils.mask2_maxcutoff(bp.io.read_frame(os.path.join(trgdir, filex1)), 
                                     bp.io.read_frame(os.path.join(refdir, filex2)), 
                                     cutoff, toprowsmask=14) 
    d = ma.masked_array(d, mask=dmask)
    mx = np.max((np.abs(np.min(d)), np.abs(np.max(d))))
    vmin = -mx
    vmax = mx

    #gs = gridspec.GridSpec(1, 30)
    #ax_hist = plt.subplot(gs[0, 0:8])
    #ax_hist.hist(ma.compressed(d), 65)

    #ax_imag = plt.subplot(gs[0, 8:26]) 
    #img = ax_imag.imshow(d, interpolation='nearest', vmin=vmin, vmax=vmax)
    #img.set_cmap('BlueRed3')
    #plt.tight_layout()
    #plt.colorbar(img, cax=plt.subplot(gs[0, 26:29]))

    fig = plt.figure()
    img = plt.imshow(d, interpolation='nearest', vmin=vmin, vmax=vmax)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    img.set_cmap('BlueRed3')
    plt.colorbar(img, cax=cax)
    plt.tight_layout()

    plt.show()
