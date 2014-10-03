import os
import sys
import numpy as np
import statsmodels.api as sm
import bopy as bp
import bopy.io as bio
from bopy.gen3fit.gen3fit import read_data, get_bkmax_phi, read_darkmean
import bopy.cgen3fit.gen3fit as g3f


def pinv_extended(X, rcond=1e-15):
    """
    Return the pinv of an array X as well as the singular values
    used in computation.

    Code adapted from statmodels.
    """
    X = np.asarray(X)
    X = X.conjugate()
    u, s, vt = np.linalg.svd(X, 0)
    s_orig = np.copy(s)
    m = u.shape[0]
    n = vt.shape[1]
    cutoff = rcond * np.maximum.reduce(s)
    for i in range(min(n, m)):
        if s[i] > cutoff:
            s[i] = 1./s[i]
        else:
            s[i] = 0.
    pinv = np.dot(np.transpose(vt), np.multiply(s[:, np.core.newaxis], np.transpose(u)))
    return pinv, s_orig

def fit(pinv, y):
    return np.dot(pinv, y)

def fit_apd(pinv, y):
    res = fit(pinv, y)
    d = res[0]
    a = np.sqrt(np.sum(res[1:]**2))
    p = np.arctan(res[2]/res[1])
    return a, p, d

def get_gen3data(filename):
    data = bio.read_frame(filename)
    shp = data.shape
    return shp, data.reshape((shp[0], shp[1]*shp[2])).T

def get_X(f, dt, n):
    """
    Get the fitting matrix

    Args:
        f -- frequency in Hz
        dt -- data rate in sec
        n -- number of data points

    Returns:
        shp -- data shape (n, nrow, ncol)
        data -- the (n, nrow*ncol) shaped data
    """
    wt = 2.*np.pi*f*np.arange(n)*dt
    cwt = np.cos(wt)
    swt = np.sin(wt)
    return np.column_stack((np.ones(n), cwt, -swt))

def get_gen3pinv(f, dt, n, rcond=1e-15):
    X = get_X(f, dt, n)
    pinv, s_orig =  pinv_extended(X, rcond=rcond)
    return pinv

def get_apd(y, X):
    model = sm.OLS(y, X)
    res = model.fit()
    d = res.params[0]
    a = np.sqrt(np.sum(res.params[1:]**2))
    p = np.arctan(res.params[2]/res.params[1])
    return a, p, d

class LinFit:
    def __init__(self, f, dt, n, rcond=1e-15):
        X = get_X(f, dt, n)
        self.pinv, self.s = pinv_extended(X, rcond=rcond) 

    def apd(self, y):
        res = np.dot(self.pinv, y)
        d = res[0]
        a = np.sqrt(np.sum(res[1:]**2))
        p = np.arctan(res[2]/res[1])
        return a, p, d

def alloc_d_forLinAPhiDC(phi_init, shape, delt, freq, nOfData):
    d = {}
    d['phi_init'] = phi_init
    d['nOfData'] = nOfData

    d['A'] = np.zeros(shape)
    d['phi'] = np.zeros(shape)
    d['DC'] = np.zeros(shape)
    d['delA'] = np.zeros(shape)
    d['delPhi'] = np.zeros(shape)
    d['delDC'] = np.zeros(shape)
    d['chi2ByDOF'] = np.zeros(shape)
    d['work'] = np.zeros((nOfData))
    d['pinv'] = get_gen3pinv(freq, delt, nOfData, rcond=1e-15)
    return d

def aphidc_flat_andro_udrot2_pickoff(datafilelist, darkfile, cfgfile, p, phi_init=None, flatfile=None, flatfiledark=None):
    """
    """
    #
    # If no dark frame, then just read the first data file 
    # and create a dark of zeros from the shape
    if darkfile is None:
        dark = read_darkmean(datafilelist[0], p)
        dark = np.zeros_like(dark, dtype=float)
    else:
        dark = read_darkmean(darkfile, p)
    np.save("darkmean", dark)

    pl = bp.utils.ParamList(cfgfile)

    #
    # flat file
    #
    if not flatfile is None:
        flat = read_darkmean(flatfile, p)
        flat_dark = read_darkmean(flatfiledark, p)
        np.save("flatdarkmean", flat_dark)
        flat = flat - flat_dark

    # slice dark for 
    #dark = dark[p.pickoff[0]:p.pickoff[1], p.pickoff[2]:p.pickoff[3]]
    #
    # Get d
    #
    d = alloc_d_forLinAPhiDC(0., dark.shape, p.delt, p.freq, p.nOfData)
    #
    # Loop over data file list
    #
    print >>sys.stderr, 'Length datafile list : ', len(datafilelist)
    for f in datafilelist:
        print >>sys.stderr, ""
        print >>sys.stderr, "--------------------------"
        print >>sys.stderr, "Data file:", f
        basename = os.path.basename(f)
        basename = os.path.splitext(basename)[0]
        data = read_data(f, p)
        # dark field correction
        data = data - dark
        # flat field correction
        if not flatfile is None:
            data = data/flat
        # start fitting
        data = data.astype(float)   # needed, because c routines are for standard double data
        #if phi_init is None:
        #    if p.pickoff_mask_rows_upto is not None:
        #        print >>sys.stderr, 'p.pickoff_mask_rows_upto:', p.pickoff_mask_rows_upto
        #        pinit = get_bkmax_phi(data[:, :p.pickoff_mask_rows_upto, :], p.delt, p.freq, 0)
        #    if p.pickoff_mask_cols_from is not None:
        #        print >>sys.stderr, 'p.pickoff_mask_cols_from:', p.pickoff_mask_cols_from
        #        pinit = get_bkmax_phi(data[:, :, p.pickoff_mask_cols_from:], p.delt, p.freq, 0, p.pickoff_mask_cols_from)
        #    else:
        #        pinit = get_bkmax_phi(data, p.delt, p.freq, 0)
        #else:
        #    pinit = phi_init
        pinit = 0.
        d['phi_init'] = pinit
        print >>sys.stderr, 'phi_init:', pinit
        g3f.fitLinAPhiDC(data, d)
        for m in p.outputs:
            output = bp.utils.data_reorient_udrot(d[m], pl)
            np.save(basename + "_" + m, output)
