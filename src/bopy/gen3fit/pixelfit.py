import sys
import numpy as np
from scipy.optimize import leastsq

#import rpy2.robjects as ro
#rpchisq = ro.r['pchisq']
#r_pchisq = ro.globalenv['pchisq']

import bopy as bp


#def pixel_cosfit(s, Pinit, delt_t, w):

def get_pcts(o, e):
    return np.sum((o-e)**2/e)

#def get_pvalue(o, e, np):
#    pcts = get_pcts(o, e)
#    dof = len(o) - 1 - np
#    return 2.*r_pchisq(pcts, dof)

def get_wt(n, dt, freq):
    return 2.*np.pi*freq*np.arange(n).astype(float)*dt

def cos_func(v, wt):
    """
    Cosines from amplitude given as square-root
    Used for fitting with amplitude constrained to be positive
    """
    return v[1]*v[1]*np.cos(wt + v[2]) + v[0]

def cos_from_aphidc(v, wt):
    """
    Cosines from regular amplitude
    """
    return v[1]*np.cos(wt + v[2]) + v[0]


def cos_dfunc(v, wt, s):
    return [np.ones_like(wt), 2.*v[1]*np.cos(wt + v[2]), -v[1]*v[1]*np.sin(wt + v[2])]

def cos_res(v, wt, s):
    return cos_func(v, wt) - s 

def pixel_cosfit2(s, wt, Pinit):
    mx = np.max(s)
    mn = np.min(s)
    Dinit = (mx+mn)/2.0
    Cinit = np.sqrt((mx-mn)/2.0)
    v0 = [Dinit, Cinit, Pinit]
    v, succ = leastsq(cos_res, v0, args=(wt, s), Dfun=cos_dfunc, col_deriv=1)
    return v[0], v[1]*v[1], v[2]

def pixel_cosfit3((i, s, wt, Pinit)):
    mx = np.max(s)
    mn = np.min(s)
    Dinit = (mx+mn)/2.0
    Cinit = np.sqrt((mx-mn)/2.0)
    v0 = [Dinit, Cinit, Pinit]
    v, succ = leastsq(cos_res, v0, args=(wt, s), Dfun=cos_dfunc, col_deriv=1)
    return i, v[0], v[1]*v[1], v[2]

def pixel_cosfit(s, Pinit, delt_t, freq):
    """ Perform cosfit on gen3 kinetic series
    """
    n = len(s)
    wt = get_wt(n, delt_t, freq)

    dfnc = lambda v: [np.ones(n), 2.*v[1]*np.cos(wt + v[2]), -v[1]*v[1]*np.sin(wt + v[2])] 
    fnc = lambda v: v[1]*v[1]*np.cos(wt + v[2]) + v[0]
    err = lambda v: fnc(v) - s

    mx = np.max(s)
    mn = np.min(s)
    Dinit = (mx+mn)/2.0
    Cinit = np.sqrt((mx-mn)/2.0)
    v0 = [Dinit, Cinit, Pinit]

    d = (wt, s)

    v, succ = leastsq(err, v0, Dfun=dfnc, col_deriv=1)
    print >>sys.stderr, "pixel_cosfit success:", succ
    return v

def pixel_fftfit(s, delt_t, freq=1.0):
    n = len(s)
    tt = np.arange(n)*delt_t
    if n < 20:
        y = s[0:10]
    if n >= 20:
        y = s[0:20]
    f, a, p = bp.utils.get_fft_fap(y, delt_t)
    idx = np.argmin(np.abs(f-freq)) 
    return a[0], a[idx], p[idx]
