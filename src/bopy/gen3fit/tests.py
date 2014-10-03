import os
import numpy as np
import bopy.io as bio

def degrees_of_freedom(data, nfits=3):
    return data.shape[0] - nfits

def get_wt(f, dt, n):
    return np.arange(n)*2.*np.pi*f*dt

def norm_data_expected(ddir, dfile, fdir, wt):

    # data
    ff = os.path.join(ddir, dfile)
    data = bio.read_frame(ff)
    # fitted values
    dark = bio.read_frame(os.path.join(fdir, 'darkmean.npy'))
    a = bio.read_frame(os.path.join(fdir, dfile.replace('.fits', '_A.npy')))
    p = bio.read_frame(os.path.join(fdir, dfile.replace('.fits', '_phi.npy')))
    d = bio.read_frame(os.path.join(fdir, dfile.replace('.fits', '_DC.npy')))
    # normalize data
    data = (data - d - dark)/a
    # expected [cos(wt+phi)]
    e = np.ones(data.shape).transpose()*wt
    e = np.cos(e.transpose() + p)
    return data, e

def norm_rss(ddir, dfile, fdir, wt):

    data, e = norm_data_expected(ddir, dfile, fdir, wt)
    data = data - e
    return np.sum(np.power(data, 2), axis=0)

def norm_tss(ddir, dfile, fdir, wt):

    data, e = norm_data_expected(ddir, dfile, fdir, wt)
    return np.sum(np.power(data, 2), axis=0)

def norm_r2(ddir, dfile, fdir, wt):
    rss = norm_rss(ddir, dfile, fdir, wt)
    tss = norm_tss(ddir, dfile, fdir, wt)
    return 1.0 - rss/tss

def norm_residual(ddir, dfile, fdir, wt):

    data, e = norm_data_expected(ddir, dfile, fdir, wt)
    return data - e

def norm_chisqr(ddir, dfile, fdir, wt):

    data, e = norm_data_expected(ddir, dfile, fdir, wt)
    return np.mean(np.power(data - e, 2)/e, axis=0)

def norm_chisqr(ddir, dfile, fdir, wt):

    data, e = norm_data_expected(ddir, dfile, fdir, wt)
    return np.mean(np.power(data - e, 2)/e, axis=0)

#def norm_r2(ddir, dfile, fdir, wt):
#
#    data, e = norm_data_expected(ddir, dfile, fdir, wt)
#    e_mean = np.mean(e, axis=0) 
#    return 1.0 - np.sum(data - e, axis=0)/np.sum(data - e_mean, axis=0)

def norm_se(ddir, dfile, fdir, wt, nfits=3):

    data, e = norm_data_expected(ddir, dfile, fdir, wt)
    dof = data.shape[0] - nfits
    return np.sqrt(np.sum(np.power(data - e, 2), axis=0)/dof)
