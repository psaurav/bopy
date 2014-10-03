#!/usr/bin/python

import os
import sys

import numpy as np
import scipy.io as sio

try:
    infile = sys.argv[1]
    outfile = os.path.basename(infile).replace('.npy', '.mat')
except:
    sys.exit("Usage: %s <input.npy>"%(os.path.basename(sys.argv[0])))

if not os.path.exists(infile):
    sys.exit("ERROR: infile does not exist - %s"%(infile))

if infile == outfile:
    sys.exit("ERROR: infile and outfile same - %s"%(infile))

d = np.load(infile)
name = os.path.basename(infile)[9:][:-4]
print >>sys.stderr, "d.shape:", d.shape, name, outfile
sio.savemat(outfile, {name: d})
    
