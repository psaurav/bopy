#!/usr/bin/python

import os
import sys
import bopy.gen3pp.phaseunwrap as punwrap
import numpy as np

try:
    tdir = sys.argv[1]
    troot = sys.argv[2]
    rdir = sys.argv[3]
    rroot = sys.argv[4]
    #wls = sys.argv[6]
except:
    print 'Usage: {0} <tdir> <troot> <rdir> <rroot> <mdir>'.format(os.path.basename(sys.argv[0]))
    sys.exit(1)

#wls = [int(w) for w in wls.strip().split()]
    

mdir = '../masks'
#troot = '20140724_G3089_2Targ_x2x2'
#rdir = '../REF2-all-gsmooth-8.0'
#rroot = '20140724_G3089_2Targ_REF2'
#mdir = '../masks03-x2x2-REF2/masks'
mroot = 'mask'
wls = [660, 690, 785, 808, 830]
srcs = np.arange(209)+1
punwrap.trgvsref(troot, rdir, rroot, mdir, mroot, wls, srcs, phasetag='mean')
