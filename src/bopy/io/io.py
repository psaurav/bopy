import os
import sys
import re
import numpy as np
import pyfits

def read_fem(filename):
    """
    """
    f = open(filename, 'r')
    s = f.readline()
    s = s.replace('[', '')
    s = s.replace(']', '')
    return np.array([float(v) for v in s.split()])

def read_gim(filename, dtype=float):
    """
    Read gim file and return a numpy array 

    Arguments:
        filename (str) -- name of gim file
        dtype -- data type (default float)

    Returns:
        d (n-dim numpy array)
    """
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

def read_fits(filename):
    """
    Read fits file

    Arguments:
        filename (str) -- name of FITS file

    Returns:
        d (numpy array) -- data
    """
    hdulist = pyfits.open(filename)
    n = len(hdulist) - 1
    d = hdulist[n].data.astype(float)
    if 'DOFFSET' in hdulist[n].header:
        doffset = float(hdulist[n].header['DOFFSET'])
        d = d + doffset
    return d

def read_frame(filename):
    """
    Read fits, npy, gim files
    """
    if re.search(r'\.gim$', filename, re.M|re.I):
        print "gim file"
        d = read_gim(filename)
    elif re.search(r'\.npy$', filename, re.M|re.I):
        print "npy file"
        d = np.load(filename)
    elif re.search(r'\.fits$', filename, re.M|re.I):
        print "fits file"
        d = read_fits(filename)


def fits_info(filename):
    """
    """
    hdulist = pyfits.open(filename)
    print >>sys.stderr, "Fits info (hdulist.info())"
    print >>sys.stderr, "-------------------------"
    hdulist.info()
    print >>sys.stderr, "-------------------------"
    print >>sys.stderr, "hdulist length:", len(hdulist)
    n = len(hdulist) - 1
    #print >>sys.stderr, hdulist[n].header
    for k,v in hdulist[n].header.items():
        print >>sys.stderr, k, '=', v

def write_bin(filename, arr):
    np.ndarray.tofile(filename, arr)

def read_slabfit_results(filename):
    """
    Read slabfit result files
    returns a self-evident dict based data-structure
    """
    f = open(filename, 'r')
    d = {}
    for l in f:
        a = [v for v in l.split()]
        src = int(a[0])
        wl = int(a[1])
        mua = float(a[3])
        musp = float(a[4])
        if src not in d:
            d[src] = {}
        d[src][wl] = {'mua': mua, 'musp': musp}
    return d
