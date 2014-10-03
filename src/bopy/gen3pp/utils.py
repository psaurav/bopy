import os
import numpy as np
import h5py
import matplotlib.pyplot as plt

def mult(x, A):
    return np.dot(A, x)

def xz2cr(xz, pl):
    dp = float(pl.get_val("CALIBRATED_CCD_DPIXEL"))
    ct = float(pl.get_val("CALIBRATED_CCD_COS_THETA"))
    st = float(pl.get_val("CALIBRATED_CCD_SIN_THETA"))
    x0 = float(pl.get_val("CALIBRATED_CCD_X0"))
    z0 = float(pl.get_val("CALIBRATED_CCD_Z0"))
    v0 = np.array([x0, z0])

    R = np.array([ct, st, -st, ct]).reshape((2,2))

    ns = len(xz)

    xz = xz/dp              # scale
    cr = np.zeros(ns*2).reshape((ns, 2))
    for i in xrange(ns):
        cr[i] = mult(xz[i], R)
    cr = cr + v0
    return cr


def create_mask(shape):
    """
    Create a mask of a desired shape (all pixels unmasked)
    """
    return np.zeros(shape).astype(bool)

def mask_annular(m, bound):
    """
    Add an annular mask, given by (rmin, rmax, cmin, cmax) 
    """
    rmin, rmax, cmin, cmax = bound
    m[:rmin,:] = True
    m[rmax:,:] = True
    m[:,:cmin] = True
    m[:,cmax:] = True
    return m

def mask_bottomleft(m, toprightpoint):
    """
    Add a rectangular mask, defined by the top right point (r, c)
    """
    r, c = toprightpoint
    m[r:, :c] = True
    return m

def point_in_poly(x, y, poly):
    """
    Taken from the internet
    """ 
    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def get_transill(ddir, rootext, wl, kind='A', collapse='mean', nrow=11, ncol=19, srcs=None):
    import os
    import numpy as np
    ddir = ddir.rstrip('/')

    if srcs is None:
        srcs = np.arange(nrow*ncol) + 1
    
    #for r in range(nrow):
    #    for c in range(ncol):
    for i in srcs:
        #i = c + ncol*r + 1
        fpath = '{0}_wl{1}_s{2}_{3}.npy'.format(rootext, wl, i, kind)
        fpath = os.path.join(ddir, fpath)
        d = np.load(fpath)
        try:
            dd = np.vstack((dd, d[np.newaxis,:]))     # stack them together
        except:
            dd = d[np.newaxis,:]

    if collapse is 'median':
        dd = np.median(dd, axis=0)
    elif collapse is 'min':
        dd = np.min(dd, axis=0)
    elif collapse is 'max':
        dd = np.max(dd, axis=0)
    else:
        dd = np.mean(dd, axis=0)
    return dd

def get_transill2(ddir, rootext, ddirr, rootextr, wl, kind='A', collapse='mean', mask=None, srcs=None):
    t = get_transill(ddir, rootext, wl, kind=kind, collapse=collapse, srcs=srcs)
    r = get_transill(ddirr, rootextr, wl, kind=kind, collapse=collapse, srcs=srcs)
    print 't max:', np.max(t)
    print 't min:', np.min(t)
    print 'r max:', np.max(r)
    print 'r min:', np.min(r)
    if mask is not None:
        return np.log(ma.masked_array(t, mask=mask)/ma.masked_array(r, mask=mask))
    else:
        return np.log(np.abs(t)/np.abs(r))

def npy2h5(fname, root, exs, wls, dkinds=['A', 'phi', 'DC'], ns=209):
    srcs = np.arange(ns) + 1
    f = h5py.File(fname, 'w')
    for dk in dkinds:
        for ex in exs:
            for w in wls:
                for s in srcs:
                    gname = '/{0}/{1}/{2}/{3}'.format(dk, ex, w, s)
                    fl = os.path.join('.', ex, '{0}_{1}_wl{2}_s{3}_{4}.npy'.format(root, ex, w, s, dk)) 
                    d = np.load(fl)
                    dset = f.create_dataset(gname, shape=d.shape, 
                                            #dtype=d.dtype, data=d) 
                                            #dtype=d.dtype, data=d, 
                                            dtype=d.dtype, 
                                            #chunks=d.shape, compression='gzip',compression_opts=9)
                                            chunks=True, compression='lzf')
                    dset[...] = d
    f.close()

def hottest_pixels(ddir, exs, root, wls, zupto=None, ns=209):

    import matplotlib.pyplot as plt
    import bopy.io as bio
    from scipy.ndimage.filters import gaussian_filter as gaussian_filter

    colrs = {660: 'b', 690: 'g', 785: 'r', 808: 'c', 830: 'm'}

    ddir = ddir.rstrip()
    srcs = np.arange(ns) + 1
    rowbnds = np.arange(11)*19
    fig = plt.figure(figsize=(15, 5), dpi=100)
    ax = plt.subplot(111)
    for ex in exs:
        for w in wls:
            max_arr = []
            for s in srcs:
                fpath =  os.path.join(ddir, ex, '{0}_{1}_wl{2}_s{3}.fits'.format(root, ex, w, s))
                d = bio.read_frame(fpath)
                if zupto is not None:
                    mx = np.max(d[:, :zupto, :])
                    #mx = np.argmax(d[:, :zupto, :])
                    #mx = np.unravel_index(mx, d[:, :zupto, :].shape)
                    #mx = np.max(gaussian_filter(d[mx[0], :zupto, :], sigma=3))
                    
                else:
                    mx = np.max(d)
                    #mx = np.argmax(d)
                    #mx = np.unravel_index(mx, d.shape)
                    #mx = np.max(gaussian_filter(d[mx[0], :, :], sigma=3))
                max_arr.append(mx)
            max_arr = np.array(max_arr)
            colr = colrs[w]
            if ex is 'REF':
                plt.plot(srcs, max_arr, 'o', color=colr, label='{0}-{1}'.format(ex, w))
            else:
                plt.plot(srcs, max_arr, 'v', color=colr, label='{0}-{1}'.format(ex, w))
    
    plt.axhline(65535, color='0.5')
    ylim = ax.get_ylim() 
    ylim_new = (0, ylim[1])
    ax.set_ylim(ylim_new) 
    plt.vlines(rowbnds, ymin=ylim_new[0], ymax=ylim_new[1], color='0.25')

    # Shink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_title(root)
    ax.set_xlabel('Source')
    ax.set_ylabel('Max CCD reading')
    ax.set_xlim((0, np.max(srcs)+1)) 
    #plt.legend()
    plt.show()


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

def compareviewallsource_masked(trgdir, refdir, trg, ref, param, compare, root, cutoff, wavelength, toprowmask, genmask, save_file_name=None):

    if genmask is not 'None':
        genmask = np.load(genmask)
        
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
