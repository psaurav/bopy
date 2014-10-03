#!/usr/bin/python

import os
import sys
import numpy as np
import bopy as bp

if __name__ == '__main__':
    try:
        cfgfilename = sys.argv[1]
        filename = sys.argv[2]
        name = sys.argv[3]
        column = int(sys.argv[4])
    except:
        sys.exit('Usage: %s <cfgfile> <file> <meas_name> <column>'%(os.path.basename(sys.argv[0])))

    pl = bp.utils.ParamList(cfgfilename)
    wl_nm = np.array([int(v) for v in pl.get_val('WAVELENGTHS').strip().split()])

    val = bp.spectra.lbs_spectra(filename, column, wl_nm)/10.0
    for w, val in zip(wl_nm, val):
        print '%s_%d = %f'%(name, w, val)
