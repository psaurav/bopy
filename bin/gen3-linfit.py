#!/usr/bin/python

import sys
import os
import numpy as np
import bopy as bp
import bopy.utils as butils
import bopy.gen3fit.linfit as blf

def listdatalist_noprepends(filelistname):
    f = open(filelistname, 'r')
    listfiles = f.readlines()
    datalist = []
    for fl in listfiles:
        fname = os.path.join(fl.rstrip())
        if not os.path.isfile(fname):
            print >>sys.stderr, "File:", fname, "does not exist"
            sys.exit(1)
        datalist.append(fname)
    darkfile = datalist.pop(0)
    return darkfile, datalist
        
def listdatalist(ddir, studyname, meas, filelistname):
    f = open(filelistname, 'r')
    listfiles = f.readlines()
    datalist = []
    for fl in listfiles:
        fname = os.path.join(ddir, studyname, meas, fl.rstrip())
        if not os.path.isfile(fname):
            print >>sys.stderr, "File:", fname, "does not exist"
            sys.exit(1)
        datalist.append(fname)
    darkfile = datalist.pop(0)
    return darkfile, datalist

def listdatalist_nodark(ddir, studyname, meas, filelistname):
    f = open(filelistname, 'r')
    listfiles = f.readlines()
    datalist = []
    for fl in listfiles:
        fname = os.path.join(ddir, studyname, meas, fl.rstrip())
        if not os.path.isfile(fname):
            print >>sys.stderr, "File:", fname, "does not exist"
            sys.exit(1)
        datalist.append(fname)
    return None, datalist

def gen3darkfile(ddir, studyname, root, meas):
    darkfile = "%s_%s_dark.fits" % (root, meas)
    darkfile = os.path.join(ddir, studyname, meas, darkfile)
    if not os.path.exists(darkfile):
        print >>sys.stderr, "File:", darkfile, "does not exist"
        sys.exit(1)
    return darkfile


def gen3datalist(ddir, studyname, root, meas, wavelengths, sources):
    datalist = []
    for s in sources:
        for w in wavelengths:
            #print s, w
            fname = "%s_%s_wl%d_s%d.fits" % (root, meas, w, s)
            fname = os.path.join(ddir, studyname, meas, fname)
            datalist.append(fname)
            #print fname
            if not os.path.exists(fname):
                print >>sys.stderr, "File:", fname, "does not exist"
                sys.exit(1)
    return datalist

def gen3linfit(cfgfilename):
    pl = bp.utils.ParamList()
    pl.read_from_file(cfgfilename)
    pl.read_from_file(pl.get_val('GEN3FIT_PREV_CFGFILE'))
    srcplate_cfg = pl.get_val('GEN3FIT_PREV_CFGFILE')

    ddir = pl.get_val('DATA_DIRECTORY')
    studyname = pl.get_val('STUDY_NAME')
    if pl.get_val('ROOT') is not None:
        root = pl.get_val('ROOT')
    else:
        root = studyname
    meas = pl.get_val('GEN3FIT_MEASUREMENT')
    fmeas = pl.get_val('GEN3FIT_FLATFILEDIR')
    flatfile = pl.get_val('GEN3FIT_FLATFILENAME')
    flatfiledark = pl.get_val('GEN3FIT_FLATFILEDARKNAME')
    
    if pl.get_val('GEN3FIT_FILENAME_ENCODING') == 'GEN3':
        wavelengths = [int(v) for v in pl.get_val('WAVELENGTHS').split()]
        sources = [int(v) for v in pl.get_val('SOURCES').split()]
        if len(sources) == 1 and sources[0] == 0:
            sources = np.arange(209)+1
        datalist = gen3datalist(ddir, studyname, root, meas, wavelengths, sources)
        darkfile = gen3darkfile(ddir, studyname, root, meas)
    elif pl.get_val('GEN3FIT_FILENAME_ENCODING') == 'LIST':
        filelistname = 'datalist-firstisdark.txt'
        if not os.path.exists(filelistname):
            print >>sys.stderr, "File", filelistname, "does not exist"
            sys.exit(1)
        darkfile, datalist = listdatalist(ddir, studyname, meas, filelistname)
    elif pl.get_val('GEN3FIT_FILENAME_ENCODING') == 'LIST_NODARK':
        filelistname = 'datalist-nodark.txt'
        if not os.path.exists(filelistname):
            print >>sys.stderr, "File", filelistname, "does not exist"
            sys.exit(1)
        darkfile, datalist = listdatalist_nodark(ddir, studyname, meas, filelistname)
    elif pl.get_val('GEN3FIT_FILENAME_ENCODING') == 'LIST_NOPREPENDS':
        filelistname = pl.get_val('GEN3FIT_FILENAME_DATALIST')
        print >>sys.stderr, "file list name:", filelistname
        if not os.path.exists(filelistname):
            print >>sys.stderr, "File", filelistname, "does not exist"
            sys.exit(1)
        darkfile, datalist = listdatalist_noprepends(filelistname)

    p = butils.struct()

    p.rebin = int(pl.get_val('REBIN'))
    p.delt = float(pl.get_val('GEN3FIT_SAMPLING_INTERVAL_SEC'))
    p.freq = float(pl.get_val('GEN3FIT_BEAT_FREQ_HZ'))
    p.nOfData = int (pl.get_val('GEN3FIT_N_DATAPOINTS'))
    p.absErr = float(pl.get_val('GEN3FIT_ABS_ERROR'))
    p.relErr = float(pl.get_val('GEN3FIT_REL_ERROR'))
    p.maxIter = int (pl.get_val('GEN3FIT_MAX_ITERATIONS'))
    p.outputs = [v for v in pl.get_val('GEN3FIT_OUTPUTS').split()]
    p.drop_first_frames = int(pl.get_val('GEN3FIT_DROP_FIRST_FRAMES'))

    try:
        p.pickoff_mask_rows_upto = int(pl.get_val('GEN3FIT_PICKOFF_MASK_ROWS_UPTO'))
        p.pickoff_mask_cols_from = None
    except:
        p.pickoff_mask_rows_upto = None
    try:
        p.pickoff_mask_cols_from = int(pl.get_val('GEN3FIT_PICKOFF_MASK_COLS_FROM'))
        p.pickoff_mask_rows_upto = None
    except:
        p.pickoff_mask_cols_from = None

    #bg3f.aphidc_flat_andro(datalist, darkfile, cfgfile, p, flatfile, flatfiledark)
    call = pl.get_val('GEN3FIT_CALL')
    if call == 'linfit_flat_andro_udrot2':
        blf.aphidc_flat_andro_udrot2_pickoff(datalist, darkfile, srcplate_cfg, p)
    elif call == 'linfit_flat_andro_udrot2_pickoff':
        p.pickoff = [int(v) for v in pl.get_val('GEN3FIT_PICKOFF').split()]
        blf.aphidc_flat_andro_udrot2_pickoff(datalist, darkfile, srcplate_cfg, p)
    else:
        print >>sys.stderr, 'Nothing to do.  Bye.'

if __name__ == '__main__':
    try:
        cfgfilename = sys.argv[1]
    except:
        sys.exit("Usage: "+os.path.basename(sys.argv[0])+" <cfgfilename>")
    if not os.path.exists(cfgfilename):
        sys.exit("ERROR: "+os.path.basename(sys.argv[0])+" file does not exist: "+cfgfilename)
    gen3linfit(cfgfilename)

