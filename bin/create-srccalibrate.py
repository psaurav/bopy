#!/usr/bin/python

import os
import sys
import numpy as np
import bopy.gen3pp.srccalibrate as scalib

try:
    cfgfile = sys.argv[1]
    trg = sys.argv[2]
    ref = sys.argv[3]
except:
    print >>sys.stderr, 'Usage: {0} <cfgfile> <trg> <ref>'.format(os.path.basename(sys.argv[0]))
    sys.exit(1)

exs = [trg, ref]
trg = './calibsrc-{0}/'.format(trg)
ref = './calibsrc-{0}/'.format(ref)
ddirs = [trg, ref]

p = scalib.get_calbsrc_params(cfgfile)

print p.calibpixels[0], p.calibpixels[1]

files = {}
for X in p.dkinds:
    files[X] = open('srccalib-{0}.txt'.format(X), 'w')

for ddir, meas in zip(ddirs, exs):
    for w in p.wls:
        for s in p.calibsrcs:
            for X in p.dkinds:
                cval = scalib.get_calibration_singlepixel(ddir, p.root, meas, p.sigma, w, X, s, p.calibpixels)
                print >>files[X], meas, w, s, cval

for X in p.dkinds:
    files[X].close()

srcs = np.arange(209)+1
f = open('src-calibsrc.txt', 'w')
for s in srcs:
    print >>f, s, p.calibsrcs[np.argmin(np.abs(s - p.calibsrcs))]
f.close()
