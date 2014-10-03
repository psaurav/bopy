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

    me = os.path.basename(sys.argv[0])
    try:
        maskdir = sys.argv[1]
    except:
        print >>sys.stderr, 'Usage: {0} <maskdir>'.format(me)
        sys.exit(1)

    nrow = 11
    ncol = 19

    images = []

    mx = -np.inf
    mn = np.inf

    #fig, ax = plt.subplots(nrow, ncol)
    gs = gridspec.GridSpec(nrow, ncol+1)
    for r in range(nrow):
        for c in range(ncol):
            i = c + ncol*r + 1
            fmask = os.path.join(maskdir, 'mask_s{0}.npy'.format(i))
            dmask = np.load(fmask)
            ax = plt.subplot(gs[r, c])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title('%d'%(i), size='x-small')
            im = ax.imshow(dmask, interpolation='nearest')
            images.append(im)

    #plt.subplots_adjust(wspace=0.1, hspace=0.0)

    #vmin = mn
    #vmax = mx

    #maximum = np.max(np.abs(mx), np.abs(mn))
    #vmin = - maximum
    #vmax = maximum

    #print "vmin:", vmin
    #print "vmax:", vmax

    for im in images:
        im.set_clim(vmin=0)
        im.set_clim(vmax=1)
        #im.set_cmap('BlueRed3')
        im.set_cmap(cm.gray)

    plt.colorbar(images[0], cax=plt.subplot(gs[:, ncol]))

    plt.show()
