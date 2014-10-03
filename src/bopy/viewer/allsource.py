import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

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

def display_raw(arr, nrow=11, ncol=19, mask=None, ofile=None):
    images = []
    mx = -np.inf
    fig = plt.figure()
    gs = gridspec.GridSpec(nrow, ncol+1)
    for r in range(nrow):
        for c in range(ncol):
            i0 = c + ncol*r
            i = i0 + 1
            d = arr[i0]
            ax = plt.subplot(gs[r, c])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title('%d'%(i), size='x-small')
            im = ax.imshow(d, interpolation='nearest')
            images.append(im)
            mx = np.max((np.abs(np.min(d)), np.abs(np.max(d)), mx))
            
    vmax = np.max(arr)
    vmin = np.min(arr)
    for im in images:
        im.set_clim(vmin=vmin)
        im.set_clim(vmax=vmax)
        #im.set_cmap('BlueRed3')
    plt.colorbar(images[0], cax=plt.subplot(gs[:, ncol]))
    if ofile is None:
        plt.show()
    else:
        plt.savefig(ofile)


def display(arr, mask=None, ofile=None, nrow=11, ncol=19, kind='phi'):

    images = []
    if mask is not None:
        arr = ma.masked_array(arr, mask=mask)

    mx = -np.inf
    mx = np.max((np.abs(np.min(arr)), np.abs(np.max(arr)), mx))

    fig = plt.figure()
    gs = gridspec.GridSpec(nrow, ncol+1)
    for r in range(nrow):
        for c in range(ncol):
            i0 = c + ncol*r
            i = i0 + 1
            d = arr[i0]
            ax = plt.subplot(gs[r, c])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title('%d'%(i), size='x-small')
            im = ax.imshow(d, interpolation='nearest')
            images.append(im)
            #mx = np.max((np.abs(np.min(d)), np.abs(np.max(d)), mx))
            
    #vmin = -mx
    #vmax = mx
    if kind is 'phi':
        vmin = -mx
        vmax = mx
    else:
        vmax = mx
        vmin = 2.0 - mx

    for im in images:
        im.set_clim(vmin=vmin)
        im.set_clim(vmax=vmax)
        im.set_cmap('BlueRed3')
    plt.colorbar(images[0], cax=plt.subplot(gs[:, ncol]))
    if ofile is None:
        plt.show()
    else:
        plt.savefig(ofile)

def read(ddir, root, wl, kind='A'):

    for i in np.arange(209)+1:
        d = np.load(os.path.join(ddir, '{0}_wl{1}_s{2}_{3}.npy'.format(root, wl, i, kind)))
        try:
            dd = np.vstack((dd, d[np.newaxis,:]))     # stack them together
        except:
            dd = d[np.newaxis,:]
    return dd

def read_op(tdir, troot, rdir, rroot, wl, kind='A'):
    t = read(tdir, troot, wl, kind=kind)
    r = read(rdir, rroot, wl, kind=kind)
    if kind is 'phi':
        return t - r 
    else:
        return t/r

def read_masks(ddir, root):

    for i in np.arange(209)+1:
        d = np.load(os.path.join(ddir, '{0}_s{1}.npy'.format(root, i)))
        try:
            dd = np.vstack((dd, d[np.newaxis,:]))     # stack them together
        except:
            dd = d[np.newaxis,:]
    return dd

def display_op(ddir, root, ddirr, rootr, wl, kind='A'):
    d = read_op(ddir, root, ddirr, rootr, wl, kind=kind)
    display(d)
    
