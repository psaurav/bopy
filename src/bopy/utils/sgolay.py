import math
import numpy as np
import scipy as sp
import scipy.ndimage as ndi

def sg_coeff2d(k_size=5, p_order=3):
    """
    Creates coefficients for 2D Savitzky-Golay filter
    Author: Jaka
    """
    ks = (k_size-1)/2
    y_grid, x_grid = np.mgrid[-ks:(ks+1),-ks:(ks+1)]
    xgf = x_grid.flatten()
    ygf = y_grid.flatten()

    if p_order == 1:
        X = np.vstack((xgf**0, xgf, ygf)).T
    elif p_order == 2:
        X = np.vstack((xgf**0, xgf, ygf, xgf**2, xgf*ygf, ygf**2)).T
    elif p_order == 3:
        X = np.vstack((xgf**0, xgf, ygf, xgf**2, xgf*ygf, ygf**2, \
            xgf**3, xgf**2*ygf, xgf*ygf**2, ygf**3)).T
    elif p_order == 4:
        X = np.vstack((xgf**0, xgf, ygf, xgf**2, xgf*ygf, ygf**2, \
            xgf**3, xgf**2*ygf, xgf*ygf**2, ygf**3, \
            xgf**4, xgf**3*ygf, xgf*ygf**3, ygf**4)).T
    else:
        raise Error('Incorrect p_order')

    return np.linalg.solve(np.dot(X.T,X), X.T)

def sgolay2d(img, k_size=5, p_order=3, out_noise=None):
    """
    Performs 2D Savitzky-Golay filter
    Author: Jaka
    """
    C = sg_coeff2d(k_size, p_order)
    C = C[0, :]
    C = C.reshape(k_size, k_size).copy()
    smoothed = ndi.filters.convolve(img, C, mode='wrap')
    if out_noise:
        return smoothed, noise_estimation(img, smoothed)
    else:
        return smoothed
