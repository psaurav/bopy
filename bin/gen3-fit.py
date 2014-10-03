#!/usr/bin/python

import sys
import os
import numpy as np
import bopy as bp
import bopy.utils as butils
import bopy.gen3fit as bg3f

try:
    cfgfilename = sys.argv[1]
except:
    sys.exit("Usage: "+os.path.basename(sys.argv[0])+" <cfgfilename>")

if not os.path.exists(cfgfilename):
    sys.exit("ERROR: "+os.path.basename(sys.argv[0])+" file does not exist: "+cfgfilename)

pl = bp.utils.ParamList()
pl.read_from_file(cfgfilename)
pl.read_from_file(pl.get_val('GEN3FIT_PREV_CFGFILE'))
srcplate_cfg = pl.get_val('GEN3FIT_PREV_CFGFILE')

ddir = pl.get_val('GEN3FIT_DDIR')
root = pl.get_val('GEN3FIT_ROOT')
meas = pl.get_val('GEN3FIT_MEASUREMENT')
sources = [int(v) for v in pl.get_val('GEN3FIT_SOURCES').split()]
if len(sources) == 1 and sources[0] == 0:
    sources = np.arange(209)+1
wavelengths = [int(v) for v in pl.get_val('GEN3FIT_WAVELENGTHS').split()]
fmeas = pl.get_val('GEN3FIT_FLATFILEDIR')
flatfile = pl.get_val('GEN3FIT_FLATFILENAME')
flatfiledark = pl.get_val('GEN3FIT_FLATFILEDARKNAME')

datalist = []
for s in sources:
    for w in wavelengths:
        #print s, w
        fname = "%s_%s_wl%d_s%d.fits" % (root, meas, w, s)
        fname = os.path.join(ddir, root, meas, fname)
        datalist.append(fname)
        #print fname
        if not os.path.exists(fname):
            print >>sys.stderr, "File:", fname, "does not exist"
            sys.exit(1)

print datalist
 
darkfile = "%s_%s_dark.fits" % (root, meas)
darkfile = os.path.join(ddir, root, meas, darkfile)
if not os.path.exists(darkfile):
    print >>sys.stderr, "File:", darkfile, "does not exist"
    sys.exit(1)

p = butils.struct()

p.rebin = int(pl.get_val('REBIN'))
p.delt = float(pl.get_val('GEN3FIT_SAMPLING_INTERVAL_SEC'))
p.freq = float(pl.get_val('GEN3FIT_BEAT_FREQ_HZ'))
p.nOfData = int (pl.get_val('GEN3FIT_N_DATAPOINTS'))
p.absErr = float(pl.get_val('GEN3FIT_ABS_ERROR'))
p.relErr = float(pl.get_val('GEN3FIT_REL_ERROR'))
p.maxIter = int (pl.get_val('GEN3FIT_MAX_ITERATIONS'))
p.outputs = [v for v in pl.get_val('GEN3FIT_OUTPUTS').split()]


#bg3f.aphidc_flat_andro(datalist, darkfile, cfgfile, p, flatfile, flatfiledark)
bg3f.aphidc_flat_andro(datalist, darkfile, srcplate_cfg, p)

