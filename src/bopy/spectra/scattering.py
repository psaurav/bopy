import os
import sys
import numpy as np

def fit_Ab(wls_nm, musp_mm):
    fit_func = lambda p, x: p[0]*(x**(-p[1]))
    err_func = lambda p, x, y: (y - fit_func(p,x))
    p = np.array([1000.0, 1.0])
    from scipy.optimize import leastsq
    ans, success = leastsq(err_func, p[:], args=(wls_nm, musp_mm))
    print "fitting success:", success
    return ans
