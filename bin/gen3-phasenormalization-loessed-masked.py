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
from bopy.gen3pp.gen3pp_utils import phase_normalization
from bopy.utils.rfuncs import loess2d

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
        maskdir = sys.argv[5]
        studyname = sys.argv[6]
        wavelength = sys.argv[7]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), \
               "<trgdir> <refdir> <trg> <ref> <maskdir> <studyname> <wavelength> [<savefilename>]"
        sys.exit(0)
    try:
        save_file_name = sys.argv[8]
    except:
        save_file_name = None

    file1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(studyname, trg, wavelength, "phi")
    file2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(studyname, ref, wavelength, "phi")
    filemask = "{0}_wl{1}_s%d_{2}.npy".format(studyname, wavelength, "mask")

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
            pt = bp.io.read_frame(os.path.join(trgdir, file1%(i)))
            pr = bp.io.read_frame(os.path.join(refdir, file2%(i)))
            pt = loess2d(pt, 0.02)
            pr = loess2d(pr, 0.02)
            dmask = np.load(os.path.join(maskdir, filemask%(i))) 
            pt = phase_normalization(pt, pr, dmask)
            d = pt - pr
            d = ma.masked_array(d, mask=dmask)
            mx = np.max((np.abs(np.min(d)), np.abs(np.max(d)), mx))
            ax = plt.subplot(gs[r, c])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title('%d'%(i), size='x-small')
            im = ax.imshow(d, interpolation='nearest')
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
    plt.colorbar(images[0], cax=plt.subplot(gs[:, ncol]))
    if save_file_name:
        fig.savefig(save_file_name, bbox_inches='tight', dpi=150)
    else:
        plt.show()
