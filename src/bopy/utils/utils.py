import numpy as np
import numpy.ma as ma
from scipy.interpolate import interp1d

def rc_bkmax(dd):
    d = np.median(dd, axis=0)
    return d.sum(1).argmax(), d.sum(0).argmax()

def bkmax(d, bb_rmax=0):
    """
    bkmax: returns the maximum of an array by picking the index from 
    a sum over rows and columns.  this prevents picking up a spurious
    maximum.
    """
    if bb_rmax == 0 or bb_rmax == d.shape[0]:
        dp = d[:bb_rmax, :]
        return d[dp.sum(1).argmax(), dp.sum(0).argmax()]
    else:
        return d[d.sum(1).argmax(), d.sum(0).argmax()]

def bkmax_bottommasked(d, rmax=0):
    if rmax == 0 or rmax == d.shape[0]:
        return bkmax(d)
    dp = d[:rmax, :]
    return d[dp.sum(1).argmax(), dp.sum(0).argmax()]

def mask_add(m1, m2):
    return m1 | m2


def mask_create(shape):
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


def sample_mask(shape, roffset, rstride, coffset, cstride):
    """
    sample_mask: returns the sampled mask
    """
    rindex, cindex = np.indices(shape)
    rmask = (rindex - roffset) % rstride == 0
    cmask = (cindex - coffset) % cstride == 0
    return rmask & cmask


def mask_bkmaxcutoff(d, cutoff, toprowsmask=0, uppercutoff=None):
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


def get_masked_array(d, mask):
    return np.ma.compressed((np.ma.masked_array(d, mask)))


def mask_toprows(m, toprowsmask):
    if toprowsmask > 0:
        m[0:toprowsmask,:] = True
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


def get_min_phase(a, p, cutoff, toprowsmask=0, uppercutoff=None, loess_span=0.02):
    mask = mask_maxcutoff(a, cutoff, toprowsmask=toprowsmask, uppercutoff=uppercutoff)
    from bopy.utils.rloess import loess2d
    p_new = loess2d(p, span=loess_span)
    p_new = get_masked_array(p_new, mask)
    return np.min(p_new)
    

def get_gen3_index(d):
    rows, columns = np.indices(d.shape)
    return d.shape[1]*rows + columns


def get_gen3_oldindex(d, mask):
    idx = get_gen3_index(d)
    return ma.compressed(ma.array(idx, mask=mask))


def interpolate1d(x, y, xs, kind='linear'):
    f = interp1d(x, y, kind=kind)
    return f(xs) 

def abfit2musp(wl_nm, musp_mm):
    import scipy
    fit_func = lambda p, x: p[0]*(x**(-p[1]))
    err_func = lambda p, x, y: (y - fit_func(p,x))
    p = np.array([1000.0, 1.0])
    p1, success = scipy.optimize.leastsq(err_func, p[:], args=(wl_nm, musp_mm))
    return p1
