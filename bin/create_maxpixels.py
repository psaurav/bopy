#!/usr/bin/python

import os
import sys
import numpy as np

try:
    row_cutoff = int(sys.argv[1])
    w = sys.argv[2]
    ddir = sys.argv[3]
    meas = sys.argv[4]
    root= sys.argv[5]
except:
    print >>sys.stderr, 'Usage: {0} <row_cutoff> <w> <ddir> <meas> <root>'.format(os.path.basename(sys.argv[0]))
    sys.exit(1)

print "Maxpixels: got in", row_cutoff

srcs = np.arange(209)+1
f = open('maxpixels.txt', 'w')
for s in srcs:
    fname = os.path.join(ddir, '{0}_{1}_wl{2}_s{3}_A.npy'.format(root, meas, w, s))
    a = np.load(fname)
    a = a[:row_cutoff, :]
    r = np.argmax(np.sum(a, 1))
    c = np.argmax(np.sum(a, 0))
    print >>f, s, r, c
f.close()

