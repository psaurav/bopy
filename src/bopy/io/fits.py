import numpy as np
import pyfits
import os
import sys

def read_fits(filename):
    """
    Read fits file

    Arguments:
        filename (str) -- name of FITS file

    Returns:
        d (numpy array) -- data
    """
    if not os.path.exists(filename):
        sys.exit("IOError: "+__name__+": File does not exist: "+filename)

    hdulist = pyfits.open(filename)
    n = len(hdulist) - 1
    d = hdulist[n].data.astype(float)
    if 'DOFFSET' in hdulist[n].header:
        doffset = float(hdulist[n].header['DOFFSET'])
        d = d + doffset
    return d

def read_fits_pars(filename):
    """
    Read fits file

    Arguments:
        filename (str) -- name of FITS file

    Returns:
        d (numpy array) -- data
    """
    if not os.path.exists(filename):
        sys.exit("IOError: "+__name__+": File does not exist: "+filename)

    hdulist = pyfits.open(filename)
    n = len(hdulist) - 1
    d = hdulist[n].data.astype(float)
    if 'DOFFSET' in hdulist[n].header:
        doffset = float(hdulist[n].header['DOFFSET'])
        d = d + doffset
    return d

def fits_info(filename):
    """
    Display FITS file info

    Arguments:
        filename (str) -- Fits file name
    
    Returns:
        None
    """
    if not os.path.exists(filename):
        sys.exit("IOError: "+__name__+": File does not exist: "+filename)

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

