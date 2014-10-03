#!/usr/bin/python

import os
import sys
import numpy as np
import bopy as bp

def print_config_example():
    print "ROOT = 20110909_G3_ref"
    print "SOURCES = 0"
    print "WAVELENGTHS = 785"
    print "AMP_FILE_FORMAT_STRING = %s_s%d_wl%s_A_ro.gim"
    print "PHA_FILE_FORMAT_STRING = %s_s%d_wl%s_phi_ro.gim"
    print "MASK_FILE_FORMAT_STRING = %s_s%d_wl%s_mask.npy"
    print "FILE_FORMAT_ORDER = source_wavelength"
    print "PP_DATA_DIR = ."
    print "MASK_CUTOFF = 0.25"
    print "TOP_ROWS_MASK = 14"
    print "MEASUREMENT_REF = REF"
    print "MEASUREMENT_TRG = LB"
    print "N_IN = 1.33"
    print "N_OUT = 1.0"
    print "SLABFIT_FUNCTION = CONTINI"
    print "OPTIMIZATION = l"
    print "SLAB_THICKNESS = 60.0"
    print "RF_FREQ_MHZ = 70"
    print "CALCULATED_MUA_785 = 0.0021429"
    print "CALCULATED_MUSP_785 = 0.8"
    print "PRINT_DIAGNOSTICS = True"

def set_values(pl):
    v = bp.utils.struct()
    v.root = pl.get_val('ROOT')
    v.srcs = np.array([int(val) for val in pl.get_val('SOURCES').split()])
    if len(v.srcs) <= 1:
        if v.srcs[0] == 0:
            v.srcs = np.arange(1, 210) 
    v.wls = [val for val in pl.get_val('WAVELENGTHS').split()]
    v.seriessum = int(pl.get_val('SLAB_SERIES_SUM'))
    v.amp_filestring = pl.get_val('AMP_FILE_FORMAT_STRING')
    v.pha_filestring = pl.get_val('PHA_FILE_FORMAT_STRING')
    v.mask_filestring = pl.get_val('MASK_FILE_FORMAT_STRING')
    v.roffset = int(pl.get_val('SAMPLE_ROWOFFSET'))
    v.coffset = int(pl.get_val('SAMPLE_COLOFFSET'))
    v.rstride = int(pl.get_val('SAMPLE_ROWSTRIDE'))
    v.cstride = int(pl.get_val('SAMPLE_COLSTRIDE'))
    v.file_format_order = pl.get_val('FILE_FORMAT_ORDER')
    v.ppdata_dir = pl.get_val('PP_DATA_DIR')
    v.mask_cutoff = float(pl.get_val('SIGNAL_MASK_CUTOFF'))
    v.signal_mask_output = bool(pl.get_val('SIGNAL_MASK_OUTPUT'))
    v.toprowsmask = float(pl.get_val('TOP_ROWS_MASK'))
    v.meas_ref = pl.get_val('MEASUREMENT_REF')
    v.meas_trg = pl.get_val('MEASUREMENT_TRG')
    v.det_ncol = int(pl.get_val('RF_DATA_NCOL'))
    v.det_nrow = int(pl.get_val('RF_DATA_NROW'))
    v.slab_l = float(pl.get_val('SLAB_THICKNESS'))
    v.diag = bool(pl.get_val('PRINT_DIAGNOSTICS'))
    return v

try:
    cfgfile = sys.argv[1]
except:
    print >>sys.stderr, "Usage:", os.path.basename(sys.argv[0]), "<cfgfile>"
    print 
    print_config_example()
    sys.exit()

pl = bp.utils.ParamList(cfgfile)
v = set_values(pl)
# get the source and detector positions
spos = bp.io.read_sdpos('sources.txt', path=v.ppdata_dir)
dets = bp.io.read_sdpos('detectors.txt', path=v.ppdata_dir)
    
masks_dir = 'masks'
if os.path.exists(masks_dir):
    if not os.path.isdir(masks_dir):
        os.remove(masks_dir)
try:
    os.makedirs(masks_dir)
except:
    print >>sys.stderr, "Directory", masks_dir, "found.  Using it"

# Init master mask to all False
detshape = (v.det_nrow, v.det_ncol)
master = np.zeros(detshape, dtype=bool)      # init to all False
grid_sample = bp.utils.sample_mask(detshape, v.roffset, v.rstride, v.coffset, v.cstride)

print v.ppdata_dir, v.meas_ref, v.meas_trg
print v.wls

for s in v.srcs:
    print >>sys.stderr, 'source:', '{0}\r'.format(s),

    signal_sample = np.ones(detshape, dtype=bool) # all true

    for w in v.wls:
        Aref = os.path.join(v.ppdata_dir, v.meas_ref, v.amp_filestring%(v.root, v.meas_ref, w, s))         
        Pref = os.path.join(v.ppdata_dir, v.meas_ref, v.pha_filestring%(v.root, v.meas_ref, w, s))         
        Atrg = os.path.join(v.ppdata_dir, v.meas_trg, v.amp_filestring%(v.root, v.meas_trg, w, s))         
        Ptrg = os.path.join(v.ppdata_dir, v.meas_trg, v.pha_filestring%(v.root, v.meas_trg, w, s))         
        signal = ~bp.utils.mask2_maxcutoff(bp.io.read_frame(Aref),
                                           bp.io.read_frame(Atrg),
                                           v.mask_cutoff,
                                           toprowsmask=v.toprowsmask)
        if v.signal_mask_output:
            np.save(os.path.join(masks_dir, 'wl%s_s%d.npy'%(w s)), signal)
            
        signal_sample = signal_sample & signal
    master = master | signal_sample
    np.save(os.path.join(masks_dir, v.mask_filestring%(v.root, s)), signal_sample)

np.save(os.path.join(masks_dir, v.root+"_master.npy"), master)
master_grid = master & grid_sample
np.save(os.path.join(masks_dir, v.root+"_master_grid.npy"), master_grid)

