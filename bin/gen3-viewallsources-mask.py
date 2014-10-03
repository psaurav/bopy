#!/usr/bin/python

import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import bopy as bp

if __name__ == '__main__':

    try:
        fileformatstring = sys.argv[1]
        fileformatstring2 = sys.argv[2]
        cutoff = float(sys.argv[3])
        cmap = sys.argv[4]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), "<fileformat> <fileformat2> <cutoff> <cmap>"
        sys.exit(0)

    if cmap == 'grey':
        cmap = cm.Greys_r
    else:
        cmap = cm.jet

    nrow = 11
    ncol = 19

    fig, ax = plt.subplots(nrow, ncol)
    for r in range(nrow):
        for c in range(ncol):
            i = c + ncol*r + 1
            d = np.load(fileformatstring%(i))
            dmask = np.load(fileformatstring2%(i))
            dmask = bp.utils.mask_maxcutoff(dmask, cutoff, toprowsmask=14)
            d = ma.masked_array(d, mask=dmask)
            ax[r,c].get_xaxis().set_visible(False)
            ax[r,c].get_yaxis().set_visible(False)
            ax[r,c].set_title('%d'%(i), size='x-small')
            #im = ax[r,c].imshow(d, interpolation='nearest', cmap=cm.Greys_r)
            ax[r,c].imshow(d, interpolation='nearest', cmap=cmap)
    plt.subplots_adjust(wspace=0.1, hspace=0.0)
    #cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])
    #fig.colorbar(im, cax)
    plt.show()
