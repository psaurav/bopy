#!/usr/bin/python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib
import numpy.ma as ma




if __name__ == '__main__':

    matplotlib.rcParams.update({'font.size': 6})

    try:
        fileformatstring = sys.argv[1]
        cmap = sys.argv[2]
    except:
        print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), "<fileformat> <cmap>"
        sys.exit(0)

    try:
        save_file_name = sys.argv[3]
        print >>sys.stderr, "Output file:", save_file_name
    except:
        save_file_name = None

    if cmap == 'grey':
        cmap = cm.Greys_r
    else:
        cmap = cm.jet

    nrow = 11
    ncol = 19

    #fig, ax = plt.subplots(nrow, ncol)
    images = []
    fig = plt.figure(figsize=(10, 8))
    grid = AxesGrid(fig, 111, # similar to subplot(131)
                    nrows_ncols = (nrow, ncol),
                    axes_pad = 0.0,
                    cbar_mode = 'single',
                    cbar_pad = 0.1,
                    #cbar_location = 'top',
                    share_all=True,
                    label_mode = 'L'
                    )
    vmin = np.inf
    vmax = -np.inf
    for r in range(nrow):
        for c in range(ncol):
            i = c + ncol*r + 1
            d = np.load(fileformatstring.format(i))
            #ax[r,c].get_xaxis().set_visible(False)
            #ax[r,c].get_yaxis().set_visible(False)
            #ax[r,c].set_title('{0}'.format(i), size='x-small')
            #im = ax[r,c].imshow(d, interpolation='nearest', cmap=cm.Greys_r)
            #ax[r,c].imshow(d, interpolation='nearest', cmap=cmap)
            im = grid[i-1].imshow(d, interpolation='nearest')
            grid[i-1].get_xaxis().set_ticks([])
            grid[i-1].get_yaxis().set_ticks([])
            grid[i-1].set_ylabel(i)
            images.append(im)
            m = np.zeros(d.shape).astype(bool)
            m[:] = False
            m[0:1,:] = True
            m[:,0:1] = True
            m[:-2,:] = True
            m[:,:-2] = True
            d = ma.masked_array(d, mask=m) 
            
            mx = np.max(d)
            if vmax < mx:
                vmax = mx
            mn = np.min(d)
            if vmin > mn:
                vmin = mn

    print "vmin:", vmin
    print "vmax:", vmax
    for im in images:
        im.set_clim(vmin=vmin)
        im.set_clim(vmax=vmax)

    
    #plt.subplots_adjust(wspace=0.1, hspace=0.0)
    #cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])
    #fig.colorbar(im, cax)
    plt.colorbar(im, cax = grid.cbar_axes[0])
    plt.tight_layout()
    if save_file_name:
        fig.savefig(save_file_name, bbox_inches='tight', dpi=300)
    else:
        plt.show()

