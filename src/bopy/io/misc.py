import os
import sys
import numpy as np
import bopy as bp

def read_frame(filename):
    """
    Read fits, npy, gim files
    """
    import numpy as np
    import re

    import os

    if re.search(r'\.gim$', filename, re.M|re.I):
        #print >>sys.stderr, "gim file"
        d = bp.io.read_gim(filename)
    elif re.search(r'\.npy$', filename, re.M|re.I):
        if not os.path.exists(filename):
            sys.exit("IOError: "+__name__+": file does not exist:"+filename)
        #print >>sys.stderr, "npy file"
        d = np.load(filename)
    elif re.search(r'\.fits$', filename, re.M|re.I):
        #print >>sys.stderr, "fits file"
        d = bp.io.read_fits(filename)
    elif re.search(r'\.fem$', filename, re.M|re.I):
        #print >>sys.stderr, "fits file"
        d = bp.io.read_fem(filename)
    return d

def write_frame(filename, d, path=None):
    if path is not None:
        np.save(os.path.join(path, filename), d)
    else:
        np.save(filename, d)

def read_mask(filename, path=None):
    if path is not None:
        return read_frame(os.path.join(path, filename)).astype(bool)
    else:
        return read_frame(filename).astype(bool)

def read_ppfile(ddir, meas, studyname, filestring, src, wl, ptype):
    fpath = os.path.join(ddir, meas, filestring%(studyname, meas, wl, src, ptype))
    return bp.io.read_frame(fpath)

def read_ppfile2(ddir, measdir, meas, studyname, filestring, src, wl, ptype):
    fpath = os.path.join(ddir, measdir, filestring.format(studyname, meas, str(wl), str(src), ptype))
    return bp.io.read_frame(fpath)

def read_float_table(filename, comments=None, skiprows=0, colstart=None, colend=None, colskip=None):
    """
    Reads tabular text data files
    """
    if colskip == 0:
        colskip = None
    if not os.path.exists(filename):
        sys.exit("IOError: "+__name__+": file does not exist:"+filename)
    return np.loadtxt(filename, comments=comments, skiprows=skiprows)[:,colstart:colend+1:colskip]

def write_bin(filename, arr):
    np.ndarray.tofile(filename, arr)
