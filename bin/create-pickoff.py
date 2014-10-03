#!/usr/bin/python

import os
import sys
import numpy as np
import bopy.gen3pp.pickoff as pickoff

try:
    cfgfile = sys.argv[1]
    trg = sys.argv[2]
    ref = sys.argv[3]
except:
    print >>sys.stderr, 'Usage: {0} <trg> <ref>'.format(os.path.basename(sys.argv[0]))
    sys.exit(0)

exs = [trg, ref]
trg = './pickoff-{0}/'.format(trg)
ref = './pickoff-{0}/'.format(ref)
ddirs = [trg, ref]

p = pickoff.get_pickoff_params(cfgfile)
 
fname='pickoff.h5'

all_dkinds = pickoff.get_alldkinds(ddirs, p.root, exs, p.sigma, p.pickoffpixels, p.wls, p.dkinds)
pickoff.save_alldkinds(fname, all_dkinds)

