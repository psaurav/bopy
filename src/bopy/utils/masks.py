import os
import numpy as np
import numpy.ma as ma

from bopy.utils.utils import bkmax_bottommasked

def add(m1, m2):
    """
    Add two masks
    """
    return m1 | m2

def create(shape):
    """
    Create a mask of a desired shape (all pixels unmasked)
    """
    return np.zeros(shape).astype(bool)

def annular(m, bound):
    """
    Add an annular mask, given by (rmin, rmax, cmin, cmax) 
    """
    rmin, rmax, cmin, cmax = bound
    m[:rmin,:] = True
    m[rmax:,:] = True
    m[:,:cmin] = True
    m[:,cmax:] = True
    return m

def bottomleft(m, toprightpoint):
    """
    Add a rectangular mask, defined by the top right point (r, c)
    """
    r, c = toprightpoint
    m[r:, :c] = True
    return m

def bottomband(m, rmax):
    """
    Add a band mask, defined by the top of the band rmax 
    """
    m[rmax:, :] = True
    return m

def genmask(shape, boundary, righttop):
    mask = create(shape)
    mask = annular(mask, boundary)
    mask = bottomleft(mask, righttop)
    return mask

def sample_mask(shape, roffset, rstride, coffset, cstride):
    """
    sample_mask: returns the sampled mask
    """
    rindex, cindex = np.indices(shape)
    rmask = (rindex - roffset) % rstride == 0
    cmask = (cindex - coffset) % cstride == 0
    return rmask & cmask


def bkmaxcutoff(d, cutoff, toprowsmask=0, uppercutoff=None):
    """
    Returns the mask based on bkmax of an array
    """
    mask = np.ones_like(d, dtype=bool)
    mask[d > cutoff*bkmax(d)] = False
    if not uppercutoff == None:
        if uppercutoff > mask_cutoff:
            mask[d > uppercutoff*bkmax(d)] = True
    if toprowsmask > 0:
        mask[0:toprowsmask,:] = True
    return mask

def bkmaxcutoff_bottommasked(d, cutoff, toprbb=0):
    """
    Returns the mask based on bkmax 
    In finding the max, a lower band is occluded.
    """
    mask = np.zeros_like(d, dtype=bool)
    maxval = bkmax_bottommasked(d, rmax=toprbb)
    mask[d < cutoff*maxval] = True
    return mask


def get_masked_array(d, mask):
    return np.ma.compressed((np.ma.masked_array(d, mask)))

def mask_toprows(m, ntoprows):
    if ntoprows> 0:
        m[:ntoprows,:] = True
    return m
    
def mask_bottomrows(m, fromrow):
    if fromrow > 0:
        m[fromrow:,:] = True
    return m
    

def mask_maxcutoff(d, cutoff, toprowsmask=0, uppercutoff=None, genmask=None):
    """
    Returns the mask based on np.max of an array
    False - signal pixel
    True - non-signal pixel
    """

    # Init mask to all True
    #print d.shape
    mask = np.ones_like(d, dtype=bool) # init to all True

    # Mark the signal as False
    mask[d > cutoff*np.max(d)] = False

    # If we want to mask the central portion (for an annular choice)
    if not uppercutoff == None:
        if uppercutoff > mask_cutoff:
            mask[d > uppercutoff*np.max(d)] = True

    # Mask toprows 
    if toprowsmask > 0:
        #toprowsmask = toprowsmask -1
        mask[0:toprowsmask,:] = True

    # General mask
    if genmask is not None:
        mask = mask | genmask

    return mask


def mask2_maxcutoff(d1, d2, cutoff, toprowsmask=0, uppercutoff=None, genmask=None):
    """
    Returns the common signal mask for two data-arrays
    True - non-signal pixel
    False - signal pixel
    """
    mask1 = mask_maxcutoff(d1, cutoff, toprowsmask=toprowsmask, uppercutoff=uppercutoff, genmask=genmask)
    mask2 = mask_maxcutoff(d2, cutoff, toprowsmask=toprowsmask, uppercutoff=uppercutoff, genmask=genmask)
    return mask1 | mask2

def unmask_noise_sea(d, pmin=0.02, pmax=0.022):
    shp = d.shape
    import scipy.ndimage
    cnt = scipy.ndimage.measurements.maximum_position(d)
    m1 = mask_maxcutoff(d, pmin, cnt[0])
    m2 = mask_maxcutoff(d, pmax, cnt[0])
    m = m1 != m2 
    m_col_indices = np.ma.nonzero(m)[1]
    col_1fourth = (np.max(m_col_indices) - np.min(m_col_indices))/4
    col_1eighth = (np.max(m_col_indices) - np.min(m_col_indices))/8
    if cnt[1] >= shp[1]/2:
        mx = cnt[1] - col_1fourth
        mn = mx - col_1eighth
    else:
        mn = cnt[1] + col_1fourth
        mx = mn + col_1eighth
    m[:,0:mn] = False
    m[:,mx::] = False
    return m == False

def create_gen3maskspo(wls, ppdir, trgdir, trg, refdir, ref, root, shape, boundary, rpo, cpo, cutoff):
    righttop = (rpo, cpo)
    gmask = genmask(shape, boundary, righttop)
    np.save('genmask2.npy', gmask)
    print 'saved genmask2.npy'

    # Creating and saving masks for each source
    if not os.path.exists('masks'):
        print 'masks does not exist'
        os.makedirs('masks')
    else:
        print 'masks exists'
        if os.path.isfile('masks'):
            print 'masks is file'
            os.remove('masks')
            os.makedirs('mask')

    print '  wls:', wls

    srcs = np.arange(209)+1
    for measdir, meas in zip([refdir, trgdir], [ref, trg]):
        print measdir, meas

    master = create(shape)
    master = np.logical_not(master)
    for s in srcs:
        m = create(shape)
        for w in wls:
            for measdir, meas in zip([refdir, trgdir], [ref, trg]):
            #for measdir, meas in zip([refdir, ], [ref, ]):
                fname = '{0}_{1}_wl{2}_s{3}_A.npy'.format(root, meas, w, s)
                fname = os.path.join(ppdir, measdir, fname)
                d = np.load(fname)
                m1 = bkmaxcutoff_bottommasked(d, cutoff, toprbb=rpo)
                m = m | m1
        m = m | gmask
        ofile = 'mask_s{0}.npy'.format(s)
        np.save(os.path.join('masks', ofile), m)
        master = master & m
    ofile = 'mask_master.npy'
    np.save(os.path.join('masks', ofile), master)
