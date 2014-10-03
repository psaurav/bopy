#!/usr/bin/python

import os
import sys
import numpy as np

from bopy.utils.utilclasses import ParamList, struct
import bopy.spectra.spectras as bspect
from bopy.toast.utils import read_mesh_bbox

def get_params(args):
    pl = ParamList()
    for f in args:
        pl.read_from_file(f)

    p = struct() 
    p.shapefile = pl.get_val('BREAST_SHAPEFILE')
    p.meshfile = pl.get_val('BREAST_MESH_FILE')
    p.hbo2 = float(pl.get_val('BREAST_BULK_CONC_HBO2_MICROMOLAR'))
    p.hb = float(pl.get_val('BREAST_BULK_CONC_HB_MICROMOLAR'))
    p.ph2o = float(pl.get_val('BREAST_BULK_PERCENTAGE_WATER'))
    p.plipid = float(pl.get_val('BREAST_BULK_PERCENTAGE_LIPID'))
    p.A = float(pl.get_val('BREAST_BULK_SCATTERING_PREFACTOR_A'))
    p.b = float(pl.get_val('BREAST_BULK_SCATTERING_POWER_B'))
    p.spect_h2o = pl.get_val('BREAST_SPECT_WATER')
    p.spect_lipid = pl.get_val('BREAST_SPECT_LIPID')
    p.spect_hbo2hb = pl.get_val('BREAST_SPECT_HBO2HB')
    p.chestwallzpos = float(pl.get_val('BREAST_CHESTWALLPOS_Z_MM'))
    p.wls = np.array(sorted([int(v) for v in pl.get_val('WAVELENGTHS').strip().split()]))
    p.homog_mua = []
    for w in p.wls:
        x = pl.get_val('CALCULATED_MUA_MM_{0}'.format(w))
        p.homog_mua.append(x)
    p.homog_musp = []
    for w in p.wls:
        x = pl.get_val('CALCULATED_MUSP_MM_{0}'.format(w))
        p.homog_musp.append(x)

    p.breast_ref= float(pl.get_val('BREAST_REFINDEX'))
    p.homog_ref= float(pl.get_val('HOMOG_REFINDEX'))
    p.grid_delta= float(pl.get_val('BREAST_GRID_DELTA_MM'))
    
    return p

if __name__ == '__main__':

    if len(sys.argv) == 1:
        print >>sys.stderr, 'Usage: {0} <one or more cfg files>'.format(os.path.basename(sys.argv[0]))
        sys.exit(1)

    p = get_params(sys.argv[1:])

    ofile='breastvaluesout.txt'

    f = open(ofile, 'w')

    print >>f, 'wavel = [{0}];'.format(';'.join([str(v) for v in p.wls]))

    ext_hb = bspect.get_eHB_mm_uM(p.wls, spectra=p.spect_hbo2hb)
    ext_hbo2 = bspect.get_eHBO2_mm_uM(p.wls, spectra=p.spect_hbo2hb)
    #ext_hbo2 = bspect.get_eHBO2_mm_uM(p.wls, spectra=p.spect_hbo2hb)
    #print ext_hb
    #ext = bspect.get_emat_tissue_mm(p.wls, blood_spec=p.spect_hbo2hb, h2o_spec=p.spect_h2o, fat_spec=p.spect_lipid)
    #print ext
    print >>f, 'extinct = [ ...'
    print >>f, '[{0}]; ...'.format(', '.join([str(v) for v in ext_hbo2]))
    print >>f, '[{0}]; ...'.format(', '.join([str(v) for v in ext_hb]))
    print >>f, '];'
    
    bbox = read_mesh_bbox(p.meshfile)
    xn = int((bbox[1,0]-bbox[0,0])/p.grid_delta)
    yn = int((bbox[1,1]-bbox[0,1])/p.grid_delta)
    zn = int((bbox[1,2]-bbox[0,2])/p.grid_delta)
    grd = [xn, yn, zn]
    print >>f, 'grd = [{0}];'.format(', '.join([str(v) for v in grd]))
    
    print >>f, 'mua_reference = [{0}];'.format(', '.join(p.homog_mua))
    print >>f, 'mus_reference = [{0}];'.format(', '.join(p.homog_musp))
    ref = np.ones(len(p.wls))*p.breast_ref
    homog_ref = np.ones(len(p.wls))*p.homog_ref
    print >>f, 'ref_reference = [{0}];'.format(', '.join([str(v) for v in homog_ref]))

    breast_mua_h20 = bspect.get_mua_fh2o_mm(p.wls, p.ph2o/100., spectra=p.spect_h2o)
    breast_mua_lipid = bspect.get_mua_ffat_mm(p.wls, p.plipid/100., spectra=p.spect_lipid)
    breast_bkg = breast_mua_h20 + breast_mua_lipid
    print >>f, 'mua_bkg = [{0}];'.format(', '.join([str(v) for v in breast_bkg]))
    print >>f, 'ref_bkg = {0}; % we set this to homogeneous for now (actual breast value would be 1.44)'.format(p.homog_ref)
    print >>f, 'C_init = [{0}];'.format(', '.join([str(v) for v in [p.hbo2, p.hb]]))
    print >>f, 'A_init = ', p.A
    print >>f, 'b_init = ', p.b
    f.close()
    print "Breast values output to file:", ofile
