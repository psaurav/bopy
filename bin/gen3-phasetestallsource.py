#!/usr/bin/python

import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import bopy as bp
from bopy.utils.rfuncs import loess2d

if __name__ == '__main__':

    try:
        trgdir = sys.argv[1]
        refdir = sys.argv[2]
        trg = sys.argv[3]
        ref = sys.argv[4]
        root = sys.argv[5]
        wavelength = sys.argv[7]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), \
               "<trgdir> <refdir> <trg> <ref> <root> <wavelength>"
        sys.exit(0)

    file1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, trg, wavelength, "phi")
    file2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, ref, wavelength, "phi")
    filex1 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, trg, wavelength, "A")
    filex2 = "{0}_{1}_wl{2}_s%d_{3}.npy".format(root, ref, wavelength, "A")

    cutoff = 0.1

    nrow = 11
    ncol = 19

    im = np.zeros((11, 19))
    images = []
    fig = plt.figure()
    for r in range(nrow):
        for c in range(ncol):
            i = c + ncol*r + 1
            d1 = bp.io.read_frame(os.path.join(trgdir, file1%(i)))
            d2 = bp.io.read_frame(os.path.join(refdir, file2%(i)))
            dmask = bp.utils.mask2_maxcutoff(bp.io.read_frame(os.path.join(trgdir, filex1%(i))), 
                                     bp.io.read_frame(os.path.join(refdir, filex2%(i))), 
                                     cutoff, toprowsmask=14) 
            d1 = loess2d(d1, span=0.02) 
            d1 = np.min(ma.compressed(ma.masked_array(d1, mask=dmask)))
            d2 = loess2d(d2, span=0.02) 
            d2 = np.min(ma.compressed(ma.masked_array(d2, mask=dmask)))
            d = d1 - d2
            im[r, c] = d

    plt.imshow(im, interpolation='nearest')
    plt.colorbar()
    plt.show()
