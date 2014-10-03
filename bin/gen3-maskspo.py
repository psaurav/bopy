#!/usr/bin/python

import os
import sys
import numpy as np
from bopy.utils.utilclasses import ParamList
import bopy.utils.masks as masks

if __name__ == '__main__':

    try:
        cfgfile1 = sys.argv[1]
        cfgfile2 = sys.argv[2]
    except:
        sys.exit('Usage: {0} <cfg_file1> <cfg_file2>'.format(os.path.basename(sys.argv[0])))

    pl = ParamList()
    print >>sys.stderr, 'Reading from file:', cfgfile1
    pl.read_from_file(cfgfile1)
    print >>sys.stderr, 'Reading from file:', cfgfile2
    pl.read_from_file(cfgfile2)

    wls = [int(w) for w in pl.get_val('WAVELENGTHS').strip().split()]
    ppdir = pl.get_val('PREPROCESS_DIRECTORY')
    trgdir = pl.get_val('TARGET_DIRECTORY')
    refdir = pl.get_val('REFERENCE_DIRECTORY')
    trg = pl.get_val('MASK_TARGET')
    ref = pl.get_val('MASK_REFERENCE')
    root = pl.get_val('ROOT')

    ncol = int(pl.get_val('RF_DATA_NCOL'))
    nrow = int(pl.get_val('RF_DATA_NROW'))
    shape = (nrow, ncol)

    rmin = int(pl.get_val('MASK_ROWMIN'))
    rmax = int(pl.get_val('MASK_ROWMAX'))
    cmin = int(pl.get_val('MASK_COLMIN'))
    cmax = int(pl.get_val('MASK_COLMAX'))

    boundary = (rmin, rmax, cmin, cmax)
    
    rpo = int(pl.get_val('PICKOFF_TOPRIGHT_ROW'))
    cpo = int(pl.get_val('PICKOFF_TOPRIGHT_COL'))
    cutoff = float(pl.get_val('MASK_CUTOFF'))
    
    masks.create_gen3maskspo(wls, ppdir, trgdir, trg, refdir, ref, root, shape, boundary, rpo, cpo, cutoff)
