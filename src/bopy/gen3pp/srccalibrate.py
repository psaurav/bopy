import os
import sys
import numpy as np
import numpy.ma as ma
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import interp1d
#import matplotlib.pyplot as plt
from bopy.utils import ParamList, struct

def get_mask(root, study, meas, w, c):
    f = os.path.join(root, meas, '{0}_{1}_wl{2}_calib{3}_{4}.npy'.format(study, meas, w, 0, 'A'))
    d0 = np.load(f)
    d0max = np.max(d0)
    print >>sys.stderr, d0max
    mask = d0 < d0max*c
    return mask

def get_calibration(root, study, meas, w, X, c):
    f = os.path.join(root, meas, '{0}_{1}_wl{2}_calib{3}_{4}.npy'.format(study, meas, w, 0, X))
    d0 = np.load(f)
    mask = get_mask(root, study, meas, w, c)
    v = ma.compressed(ma.masked_array(d0, mask=mask))
    print >>sys.stderr, len(v)
    m = np.mean(v)
    calib = [m,]
    for i in range(1, 12):
        f = os.path.join(root, meas, '{0}_{1}_wl{2}_calib{3}_{4}.npy'.format(study, meas, w, i, X))
        d = np.load(f)
        v = ma.compressed(ma.masked_array(d, mask=mask))
        calib.append(np.mean(v))
    return np.array(calib)

def get_calibration_singlepixel(ddir, root, meas, sigma, w, X, s, (r, c)):
    f = os.path.join(ddir, '{0}_{1}_wl{2}_calib{3}_{4}.npy'.format(root, meas, w, s, X))
    d0 = np.load(f)
    d0 = gaussian_filter(d0, sigma=sigma)
    return d0[r, c] 

def create_calibration_vector(calib, interpolkind='linear'):
    xt = np.arange(12).astype(float)*19
    print >>sys.stderr, xt
    print >>sys.stderr, len(xt), len(calib)
    f = interp1d(xt, calib, kind=interpolkind)
    src_index = np.arange(209).astype(float)
    return f(src_index)

def create_calibration_vectors(root, study, meast, measr, w, X, c):
    calib = get_calibration(root, study, measr, w, X, c)
    calibr = create_calibration_vector(calib, interpolkind='linear')
    calib = get_calibration(root, study, meast, w, X, c)
    calibt = create_calibration_vector(calib, interpolkind='linear')
    norm = calibr[0]
    if X == 'phi':
        calibr = calibr - norm
        calibt = calibt - norm
    else:
        calibr = calibr/norm
        calibt = calibt/norm
    return (calibt, calibr)

def read_src_association(ddir, fname):
    fname = os.path.join(ddir, fname)
    f = open(fname)
    mp = {}
    for l in f:
        ss = [int(v) for v in l.strip().split()]
        mp[ss[0]] = ss[1]
    return mp

def read_csrc_values(ddir, fname):
    fname = os.path.join(ddir, fname)
    f = open(fname, 'r')
    d = {}
    for l in f:
        v = [w for w in l.strip().split()]
        d.setdefault(v[0], {}).setdefault(int(v[1]), {}).setdefault(int(v[2]), {})
        d[v[0]][int(v[1])][int(v[2])] = float(v[3])
    return d

def get_maxpixels(ddir, root, meas, po_rowmax, w):
    row_cutoff = 110
    srcs = np.arange(209)+1
    for s in srcs:
        fname = os.path.join(ddir, '{0}_{1}_wl{2}_s{3}_A.npy'.format(root, meas, w, s))
        a = np.load(fname)
        a = a[:po_rowmax, :]
        r = np.argmax(np.sum(a, 1))
        c = np.argmax(np.sum(a, 0))
        print s, r, c

def read_maxpixels(ddir, fname):
    fname = os.path.join(ddir, fname)
    f = open(fname, 'r')
    d = {}
    for l in f:
        a = [int(v) for v in l.strip().split()]
        d[a[0]] = (a[1], a[2])
    return d

def recalibrate_phi(data, calib, i):
    return data - calib[i]

def recalibrate_ADC(data, calib, i):
    return data/calib[i]

def get_calbsrc_params(plfname):
    pl = ParamList(plfname)
    p = struct()
    p.wls = [int(w) for w in pl.get_val('WAVELENGTHS').strip().split()]
    p.sigma = float(pl.get_val('SIGMA'))
    p.calibpixels = tuple([int(s) for s in pl.get_val('CALIB_PIXELS').strip().split()])
    p.dkinds = [v for v in pl.get_val('DKINDS').strip().split()]
    p.calibsrcs = np.array([int(v) for v in pl.get_val('CALIB_SRCS').strip().split()])
    p.root = pl.get_val('ROOT')

    return p

def create_calibssrc(p):
    files = {}
    for X in p.dkinds:
        files[X] = open('srccalib-{0}.txt'.format(X), 'w')

    print p.calibpixels

    for ddir, meas in zip(p.ddirs, p.exs):
        for w in p.wls:
            for s in p.calibsrcs:
                for X in p.dkinds:
                    cval = get_calibration_singlepixel(ddir, root, meas, sigma, w, X, s, p.calibpixels)
                    print >>files[X], meas, w, s, cval

    for X in dkinds:
        files[X].close()

    srcs = np.arange(209)+1
    f = open('src-calibsrc.txt', 'w')
    for s in srcs:
        print >>f, s, p.calibsrcs[np.argmin(np.abs(s - p.calibsrcs))]
    f.close()

