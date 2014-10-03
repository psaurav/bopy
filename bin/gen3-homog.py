#!/usr/bin/python

import sys
import os
import bopy as bp

try:
    cfgfilename = sys.argv[1]
except:
    print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), "<cfgfilename>"
    sys.exit()

print >>sys.stderr, "cfgfilename:", cfgfilename

homog = bp.spectra.HomogMuaMusp(cfgfilename)
homog.get_mua()
homog.get_musp()
homog.output_values()
