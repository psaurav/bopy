#!/usr/bin/python

import os
import sys
import numpy as np
import bopy.gen3pp.srccalibrate as scalib
import bopy.gen3pp.pickoff as po
from bopy.utils import ParamList, struct

def parse_cfg(cfgfile):
    pl = ParamList(cfgfile)
    p = struct()
    p.trg = pl.get_val('TRG_NAME')
    p.ref = pl.get_val('REF_NAME')
    p.sigma = pl.get_val('SIGMA')
    p.new_dir_ext = 'all-gsmooth-{0}'.format(p.sigma)
    p.old_dir_ext = 'gsmooth-{0}'.format(p.sigma)
    p.root = pl.get_val('ROOT')
    p.srcs = np.array([int(v) for v in pl.get_val('SOURCES').strip().split()])
    p.wls_calib_a = [int(v) for v in pl.get_val('WLS_CALIB_AMPLITUDE').strip().split()] or []
    #p.wls_calib_p = [int(v) for v in pl.get_val('WLS_CALIB_PHASE').strip().split()] or []
    p.wls_pickoff_a = [int(v) for v in pl.get_val('WLS_PICKOFF_AMPLITUDE').strip().split()] or []
    p.wls_pickoff_p = [int(v) for v in pl.get_val('WLS_PICKOFF_PHASE').strip().split()] or []
    return p

#ref = 'REF2'
#trg = 'x2x2'
#new_dir_ext = 'all-gsmooth-8.0'
#old_dir_ext = 'gsmooth-8.0'
#root = '20140724_G3089_2Targ'

#srcs = np.arange(209)+1

## calibration source correction
try:
    cfgfile = sys.argv[1]
except:
    print >>sys.stderr, 'Usage: {0} <cfgfile>'.format(os.path.basename(sys.argv[0]))
    sys.exit(1)

print "cfgfile:", cfgfile
p = parse_cfg(cfgfile)

wls = p.wls_calib_a
calibsrc_dir = '.'
map_file = 'src-calibsrc.txt'
mp = scalib.read_src_association(calibsrc_dir, map_file)

maxpixels = scalib.read_maxpixels('.', 'maxpixels.txt')
print maxpixels

calibvals_a = scalib.read_csrc_values(calibsrc_dir, 'srccalib-A.txt')
calibvals_phi = scalib.read_csrc_values(calibsrc_dir, 'srccalib-phi.txt')
for w in wls:
    a0 = calibvals_a[p.ref][w][0]
    phi0 = calibvals_phi[p.ref][w][0]
    print a0, phi0
    for meas in [p.ref, p.trg]:
        newdirname = '{0}-{1}'.format(meas, p.new_dir_ext)
        olddirname = '{0}-{1}'.format(meas, p.old_dir_ext)
        try:
            os.mkdir(newdirname)
        except:
            pass
        for s in p.srcs:
            # amplitude
            ifile = os.path.join(olddirname, '{0}_{1}_wl{2}_s{3}_A.npy'.format(p.root, meas, w, s))
            a = np.load(ifile)
            cs = mp[s]
            a = a*calibvals_a[meas][w][0]/calibvals_a[meas][w][cs]
            ofile = os.path.join(newdirname, '{0}_{1}_wl{2}_s{3}_A.npy'.format(p.root, meas, w, s))
            np.save(ofile, a)

# pickoff
#wls_phi = [660, 690, 785, 808, 830]
wls_phi = [660, 690, 785, 808, 830]
pofile = 'pickoff.h5'
pofile = os.path.join(calibsrc_dir, pofile)
podata = po.load_alldkinds(pofile)
for w in wls_phi:
    a0 = podata['A'][p.ref][w][0]
    phi0 = podata['phi'][p.ref][w][0]
    d0 = podata['DC'][p.ref][w][0]
    print a0, phi0
    for meas in [p.ref, p.trg]:
        newdirname = '{0}-{1}'.format(meas, p.new_dir_ext)
        olddirname = '{0}-{1}'.format(meas, p.old_dir_ext)
        try:
            os.mkdir(newdirname)
        except:
            pass
        for s in p.srcs:
            # phase
            ifile = os.path.join(olddirname, '{0}_{1}_wl{2}_s{3}_phi.npy'.format(p.root, meas, w, s))
            phi = np.load(ifile)
            corr = podata['phi'][meas][w][s-1] - podata['phi'][meas][w][0]
            if corr > np.pi:
                corr -= np.pi
            elif corr < -np.pi:
                corr += np.pi
            phi = phi - corr
            ofile = os.path.join(newdirname, '{0}_{1}_wl{2}_s{3}_phi.npy'.format(p.root, meas, w, s))
            np.save(ofile, phi)
            # DC (by default)
            ifile = os.path.join(olddirname, '{0}_{1}_wl{2}_s{3}_DC.npy'.format(p.root, meas, w, s))
            dc = np.load(ifile)
            corr = podata['DC'][meas][w][0]/podata['DC'][meas][w][s-1]
            dc = np.load(ifile)
            dc = dc*corr
            ofile = os.path.join(newdirname, '{0}_{1}_wl{2}_s{3}_DC.npy'.format(p.root, meas, w, s))
            np.save(ofile, dc)

#wls_a = [660, 690, 808]
wls_a = p.wls_pickoff_a 
for w in wls_a:
    a0 = podata['A'][p.ref][w][0]
    phi0 = podata['phi'][p.ref][w][0]
    print a0, phi0
    for meas in [p.ref, p.trg]:
        newdirname = '{0}-{1}'.format(meas, p.new_dir_ext)
        olddirname = '{0}-{1}'.format(meas, p.old_dir_ext)
        try:
            os.mkdir(newdirname)
        except:
            pass
        for s in p.srcs:
            # amplitude
            ifile = os.path.join(olddirname, '{0}_{1}_wl{2}_s{3}_A.npy'.format(p.root, meas, w, s))
            phi = np.load(ifile)
            corr = podata['A'][meas][w][0]/podata['A'][meas][w][s-1]
            a = np.load(ifile)
            a = a*corr
            ofile = os.path.join(newdirname, '{0}_{1}_wl{2}_s{3}_A.npy'.format(p.root, meas, w, s))
            np.save(ofile, a)

#print calibvals_a['REF2'][785]
