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
        file1 = sys.argv[1]
        file2 = sys.argv[2]
        cutoff = float(sys.argv[3])        
        cmap = sys.argv[4]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), "<fileformat> <fileformat> <cutoff> <cmap>"
        sys.exit(0)

    if cmap == 'grey':
        cmap = cm.Greys_r
    else:
        cmap = cm.jet

    d = np.load(file1)
    dmask = np.load(file2)
    dmask = bp.utils.mask_maxcutoff(dmask, cutoff, toprowsmask=14)
    d = ma.masked_array(d, mask=dmask)
    plt.imshow(d, interpolation='nearest', cmap=cmap)
    plt.colorbar()
    plt.show()
