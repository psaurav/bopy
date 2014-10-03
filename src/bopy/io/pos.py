import os
import numpy as np

def read_sdpos2d(filename, path=None):
    if path is not None:
        d = np.loadtxt(os.path.join(path, filename))
    else:
        d = np.loadtxt(filename)
    nrows = int(np.max(d[:,2])) + 1
    ncols = int(np.max(d[:,1])) + 1
    x = d[:3].reshape((nrows, ncols))
    y = d[:4].reshape((nrows, ncols))
    z = d[:5].reshape((nrows, ncols))
    return x, y, z

#def read_sdpos(filename, path=None):
#    if path is not None:
#        d = np.loadtxt(os.path.join(path, filename))
#    else:
#        d = np.loadtxt(filename)
#    return d[:,3], d[:,4], d[:,5]

def read_sdpos(*args):
    d = np.loadtxt(os.path.join(*args))
    return np.vstack((d[:,3], d[:,4], d[:,5])).T

def read_sdpos_withindex(*args):
    d = np.loadtxt(os.path.join(*args))
    return np.vstack((d[:,0], d[:,3], d[:,4], d[:,5])).T
