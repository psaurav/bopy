import numpy as np
import os
import sys

def read_gim(filename, dtype=float):
    """
    Read gim file and return a numpy array 

    Arguments:
        filename (str) -- name of gim file
        dtype -- data type (default float)

    Returns:
        d (n-dim numpy array)
    """
    if not os.path.exists(filename):
        sys.exit("IOError: "+__name__+": File does not exist: "+filename)

    f = open(filename, 'rt')
    n = np.array([int(v) for v in f.readline().split()])
    #d = np.array([float(v) for v in f.readline().split()])
    d = np.array([v for v in f.readline().split()]).astype(dtype)
    f.close()
    d = np.reshape(d, n[::-1])
    return d

def write_gim(filename, d):
    """
    Write numpy array to a gim file (n-dim)

    Arguments:
        filename (str) -- name of gim file
        d (numpy array) -- data
    """
    f = open(filename, 'w')
    n = np.shape(d)
    s = ' '.join(map(str, n[::-1])) + '\n'
    f.write(s)
    d = d.flatten()
    d = ' '.join(map(str, d))
    f.write(d)
    f.close()

