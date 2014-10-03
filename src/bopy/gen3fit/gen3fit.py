import os
import sys
import re
import numpy as np
import bopy as bp
import bopy.io as io
import bopy.cgen3fit.gen3fit as g3f
from bopy.utils import struct, rc_bkmax
from bopy.utils import rebindata3d
from . import pixelfit as pxfit

def get_bkmax_phi(data, delt, freq, ch=0, cols_from=0):
    r,c = rc_bkmax(data)
    rr = range(r-ch, r+ch+1)
    cr = range(c-ch, c+ch+1)
    print "bkmax_phi r, c, ch:", rr, cr, ch
    phi_init = []
    for r in rr:
        for c in cr: 
            fit = pxfit.pixel_fftfit(data[:,r, c], delt, freq)
            phi_init.append(fit[2])
    return np.mean(phi_init)

def read_data(filename, p):
    """
    Read the dark file and create the mean frame
    """
    print >>sys.stderr, __name__, ": reading data..."
    d = io.read_frame(filename)
    print >>sys.stderr, __name__, ": shape(data):", d.shape

    #
    # drop first frames
    #
    if not p.drop_first_frames == None:
        d = d[p.drop_first_frames:,:,:]
    #
    # rebin the data
    #
    if p.rebin is not None:
        if p.rebin > 1:
            d = rebindata3d(d, prebin)
            print >>sys.stderr, __name__, ": rebin shape(data):", d.shape
    #
    # return
    #
    return d


def read_darkmean(filename, p):
    """
    Read the dark file and create the mean frame
    """
    #
    # Read the data
    #
    dark = read_data(filename, p)
    #
    # get the mean along axis=0
    #
    dark = dark.astype(float)
    d = np.mean(dark, dtype=np.float64, axis=0)
    print >>sys.stderr, __name__, ": shape(dark_mean):", d.shape
    #
    # save darkmean for reference and return
    #
    #np.save("darkmean.npy", d)
    return d


def aphidc_dark(datafilelist, darkfile, p):
    """
    """
    #
    dark = read_darkmean(darkfile, p)
    #
    # Set up data dictionary
    #
    d = {}
    d['delt'] = p.delt
    d['freq'] = p.freq
    d['nOfData'] = p.nOfData
    d['absErr'] = p.absErr
    d['relErr'] = p.relErr
    d['maxIter'] = p.maxIter
    d['A'] = np.zeros_like(dark)
    d['phi'] = np.zeros_like(dark)
    d['DC'] = np.zeros_like(dark)
    d['delA'] = np.zeros_like(dark)
    d['delPhi'] = np.zeros_like(dark)
    d['delDC'] = np.zeros_like(dark)
    d['iters'] = np.zeros_like(dark, dtype=np.uint)
    d['status'] = np.zeros_like(dark, dtype=np.uint)
    d['chi2ByDOF'] = np.zeros_like(dark, dtype=np.uint)
    #
    # Loop over data file list
    #
    for f in datafilelist:
        print >>sys.stderr, ""
        print >>sys.stderr, "--------------------------"
        print >>sys.stderr, "Data file:", f
        basename = os.path.basename(f)
        basename = os.path.splitext(basename)[0] 
        data = read_data(f, p)
        # Substract dark
        data = data - dark
        data = data.astype(float)
        # start fitting
        g3f.fitAPhiDC(data, d)
        np.save(basename + "_A", d['A'])
        np.save(basename + "_phi", d['phi'])
        np.save(basename + "_DC", d['DC'])
        np.save(basename + "_iters", d['iters'])
        np.save(basename + "_status", d['status'])
        np.save(basename + "_delA", d['delA'])
        np.save(basename + "_delPhi", d['delPhi'])
        np.save(basename + "_delDC", d['delDC'])
        np.save(basename + "_chi2ByDOF", d['chi2ByDOF'])
        print >>sys.stderr, "--------------------------"

def aphidc_flat(datafilelist, darkfile, flatfile, p):
    """
    """
    #
    dark = read_darkmean(darkfile, p)
    #
    # Set up data dictionary
    #
    d = {}
    d['delt'] = p.delt
    d['freq'] = p.freq
    d['nOfData'] = p.nOfData
    d['absErr'] = p.absErr
    d['relErr'] = p.relErr
    d['maxIter'] = p.maxIter
    d['A'] = np.zeros_like(dark)
    d['phi'] = np.zeros_like(dark)
    d['DC'] = np.zeros_like(dark)
    d['delA'] = np.zeros_like(dark)
    d['delPhi'] = np.zeros_like(dark)
    d['delDC'] = np.zeros_like(dark)
    d['iters'] = np.zeros_like(dark, dtype=np.uint)
    d['status'] = np.zeros_like(dark, dtype=np.uint)
    d['chi2ByDOF'] = np.zeros_like(dark, dtype=np.uint)

    #
    # flat file
    #
    flat = np.load(flatfile)
    if p.rebin is not None:
        if p.rebin > 0:
            flat = rebindata2d(flat, p)
            print >>sys.stderr, __name__, ": rebin shape(flat):", flat.shape
    #
    # Loop over data file list
    #
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
        data = data/flat
        # start fitting
        data = data.astype(float)   # needed, because c routines are for standard double data
        g3f.fitAPhiDC(data, d)
        np.save(basename + "_A", d['A'])
        np.save(basename + "_phi", d['phi'])
        np.save(basename + "_DC", d['DC'])
        np.save(basename + "_iters", d['iters'])
        np.save(basename + "_status", d['status'])
        np.save(basename + "_delA", d['delA'])
        np.save(basename + "_delPhi", d['delPhi'])
        np.save(basename + "_delDC", d['delDC'])
        np.save(basename + "_chi2ByDOF", d['chi2ByDOF'])
        print >>sys.stderr, "--------------------------"

def aphidc_flat_andro(datafilelist, darkfile, cfgfile, p, flatfile=None, flatfiledark=None):
    """
    """
    #
    dark = read_darkmean(darkfile, p)
    np.save("darkmean", dark)

    pl = bp.utils.ParamList(cfgfile)
    #
    # Set up data dictionary
    #
    d = {}
    d['delt'] = p.delt
    d['freq'] = p.freq
    d['nOfData'] = p.nOfData
    d['absErr'] = p.absErr
    d['relErr'] = p.relErr
    d['maxIter'] = p.maxIter
    d['A'] = np.zeros_like(dark)
    d['phi'] = np.zeros_like(dark)
    d['DC'] = np.zeros_like(dark)
    d['delA'] = np.zeros_like(dark)
    d['delPhi'] = np.zeros_like(dark)
    d['delDC'] = np.zeros_like(dark)
    d['iters'] = np.zeros_like(dark, dtype=np.uint)
    d['status'] = np.zeros_like(dark, dtype=np.uint)
    d['chi2ByDOF'] = np.zeros_like(dark, dtype=np.uint)

    #
    # flat file
    #
    if not flatfile is None:
        flat = read_darkmean(flatfile, p)
        flat_dark = read_darkmean(flatfiledark, p)
        np.save("flatdarkmean", flat_dark)
        flat = flat - flat_dark
    #
    # Loop over data file list
    #
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
        g3f.fitAPhiDC(data, d)
        for m in p.outputs: 
            output = bp.utils.data_reorient2(d[m], pl)
            np.save(basename + "_" + m, output)

def aphidc_flat_andro_udrot(datafilelist, darkfile, cfgfile, p, flatfile=None, flatfiledark=None):
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
    # Set up data dictionary
    #
    d = {}
    d['delt'] = p.delt
    d['freq'] = p.freq
    d['nOfData'] = p.nOfData
    d['absErr'] = p.absErr
    d['relErr'] = p.relErr
    d['maxIter'] = p.maxIter
    d['A'] = np.zeros_like(dark)
    d['phi'] = np.zeros_like(dark)
    d['DC'] = np.zeros_like(dark)
    d['delA'] = np.zeros_like(dark)
    d['delPhi'] = np.zeros_like(dark)
    d['delDC'] = np.zeros_like(dark)
    d['iters'] = np.zeros_like(dark, dtype=np.uint)
    d['status'] = np.zeros_like(dark, dtype=np.uint)
    d['chi2ByDOF'] = np.zeros_like(dark, dtype=np.uint)

    #
    # flat file
    #
    if not flatfile is None:
        flat = read_darkmean(flatfile, p)
        flat_dark = read_darkmean(flatfiledark, p)
        np.save("flatdarkmean", flat_dark)
        flat = flat - flat_dark
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
        g3f.fitAPhiDC(data, d)
        for m in p.outputs: 
            output = bp.utils.data_reorient_udrot(d[m], pl)
            np.save(basename + "_" + m, output)

def alloc_d_forAPhiDC2(phi_init, shape, delt, freq, nOfData, absErr, relErr, maxIter):
    d = {}
    d['delt'] = delt
    d['freq'] = freq
    d['phi_init'] = phi_init
    d['nOfData'] = nOfData
    d['absErr'] = absErr
    d['relErr'] = relErr
    d['maxIter'] = maxIter

    d['A'] = np.zeros(shape)
    d['phi'] = np.zeros(shape)
    d['DC'] = np.zeros(shape)
    d['delA'] = np.zeros(shape)
    d['delPhi'] = np.zeros(shape)
    d['delDC'] = np.zeros(shape)
    d['iters'] = np.zeros(shape, dtype=np.uint)
    d['status'] = np.zeros(shape, dtype=np.uint)
    d['chi2ByDOF'] = np.zeros(shape)
    d['pcts'] = np.zeros(shape)
    d['rse'] = np.zeros(shape)
    d['nrse'] = np.zeros(shape)
    return d

def getAPhiDC2(data, phi_init, shape=None, delt=0.1, freq=1.0, nOfData=17, absErr=0.0, relErr=1.0e-5, maxIter=200):
    if shape is None:
        shape = data.shape[-2:]

    d = {}
    d['delt'] = delt
    d['freq'] = freq
    d['phi_init'] = phi_init
    d['nOfData'] = nOfData
    d['absErr'] = absErr
    d['relErr'] = relErr
    d['maxIter'] = maxIter

    d['A'] = np.zeros(shape)
    d['phi'] = np.zeros(shape)
    d['DC'] = np.zeros(shape)
    d['delA'] = np.zeros(shape)
    d['delPhi'] = np.zeros(shape)
    d['delDC'] = np.zeros(shape)
    d['iters'] = np.zeros(shape, dtype=np.uint)
    d['status'] = np.zeros(shape, dtype=np.uint)
    d['chi2ByDOF'] = np.zeros(shape)
    d['pcts'] = np.zeros(shape)
    d['rse'] = np.zeros(shape)
    d['nrse'] = np.zeros(shape)
    
    g3f.fitAPhiDC2(data, d)
    return d

def aphidc_flat_andro_udrot2(datafilelist, darkfile, cfgfile, p, phi_init=None, flatfile=None, flatfiledark=None):
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

    #
    # Get d
    #
    d = alloc_d_forAPhiDC2(0., dark.shape, p.delt, p.freq, p.nOfData, p.absErr, p.relErr, p.maxIter)
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
        if phi_init is None:
            if p.pickoff_mask_rows_upto is not None:
                print >>sys.stderr, 'p.pickoff_mask_rows_upto:', p.pickoff_mask_rows_upto
                pinit = get_bkmax_phi(data[:, :p.pickoff_mask_rows_upto, :], p.delt, p.freq, 0) 
            if p.pickoff_mask_cols_from is not None:
                print >>sys.stderr, 'p.pickoff_mask_cols_from:', p.pickoff_mask_cols_from 
                pinit = get_bkmax_phi(data[:, :, p.pickoff_mask_cols_from:], p.delt, p.freq, 0, p.pickoff_mask_cols_from) 
            else:
                pinit = get_bkmax_phi(data, p.delt, p.freq, 0)
            #r,c = rc_bkmax(data)
            #print r, c
            #fit = pxfit.pixel_fftfit(data[:,r, c], p.delt, p.freq)
            #phi_init = fit[2]
        else:
            pinit = phi_init
        d['phi_init'] = pinit
        print >>sys.stderr, 'phi_init:', pinit
        #d = getAPhiDC2(data, phi_init, shape=None, delt=p.delt, freq=p.freq, nOfData=p.nOfData, absErr=p.absErr, relErr=p.relErr, maxIter=p.maxIter)
        g3f.fitAPhiDC2(data, d)
        for m in p.outputs:
            output = bp.utils.data_reorient_udrot(d[m], pl)
            np.save(basename + "_" + m, output)


def getPixelAPhiDC2(data, phi_init, shape=None, row=0, col=0, delt=0.1, freq=1.0, nOfData=17, absErr=0.0, relErr=1.0e-5, maxIter=200):
    d = {}
    d['delt'] = delt
    d['freq'] = freq
    d['phi_init'] = phi_init
    d['nOfData'] = nOfData
    d['absErr'] = absErr
    d['relErr'] = relErr
    d['maxIter'] = maxIter
    d['row'] = row
    d['col'] = col

    if shape is None:
        shape = data.shape[-2:]
    
    d['A'] = 0.
    d['phi'] = 0.
    d['DC'] = 0.
    d['delA'] = 0.
    d['delPhi'] = 0.
    d['delDC'] = 0.
    d['iters'] = 0.
    d['status'] = 0.
    d['chi2ByDOF'] = 0.
    
    g3f.fitPixelAPhiDC2(data, d)
    return d

def getPixelAPhiDCTr(data, phi_init, shape=None, row=0, col=0, delt=0.1, freq=1.0, nOfData=17, absErr=0.0, relErr=1.0e-5, maxIter=200):
    d = {}
    d['delt'] = delt
    d['freq'] = freq
    d['phi_init'] = phi_init
    d['nOfData'] = nOfData
    d['absErr'] = absErr
    d['relErr'] = relErr
    d['maxIter'] = maxIter
    d['row'] = row
    d['col'] = col

    if shape is None:
        shape = data.shape[-2:]

    d['A'] = 0.
    d['phi'] = 0.
    d['DC'] = 0.
    d['delA'] = 0.
    d['delPhi'] = 0.
    d['delDC'] = 0.
    d['iters'] = 0.
    d['status'] = 0.
    d['chi2ByDOF'] = 0.

    g3f.fitPixelAPhiDCTr(data, d)
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
    dark = dark[p.pickoff[0]:p.pickoff[1], p.pickoff[2]:p.pickoff[3]]
    #
    # Get d
    #
    d = alloc_d_forAPhiDC2(0., dark.shape, p.delt, p.freq, p.nOfData, p.absErr, p.relErr, p.maxIter)
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
        data = data[:,p.pickoff[0]:p.pickoff[1], p.pickoff[2]:p.pickoff[3]]
        # dark field correction
        data = data - dark
        # flat field correction
        if not flatfile is None:
            data = data/flat
        # start fitting
        data = data.astype(float)   # needed, because c routines are for standard double data
        if phi_init is None:
            pinit = get_bkmax_phi(data, p.delt, p.freq, 0)
            #r,c = rc_bkmax(data)
            #print r, c
            #fit = pxfit.pixel_fftfit(data[:,r, c], p.delt, p.freq)
            #phi_init = fit[2]
        else:
            pinit = phi_init
        d['phi_init'] = pinit
        print >>sys.stderr, 'phi_init:', pinit
        #d = getAPhiDC2(data, phi_init, shape=None, delt=p.delt, freq=p.freq, nOfData=p.nOfData, absErr=p.absErr, relErr=p.relErr, maxIter=p.maxIter)
        g3f.fitAPhiDC2(data, d)
        for m in p.outputs:
            output = bp.utils.data_reorient_udrot(d[m], pl)
            np.save(basename + "_" + m, output)

def alloc_d_forLinAPhiDC(phi_init, shape, delt, freq, nOfData, absErr, relErr, maxIter):
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
    return d

