#!/usr/bin/python

import sys
import os
import numpy as np
import bopy as bp
import bopy.utils as butils
import bopy.gen3fit as bg3f
        
def listdatalist(ddir, root, meas, filelistname):
    f = open(filelistname, 'r')
    listfiles = f.readlines()
    datalist = []
    for fl in listfiles:
        fname = os.path.join(ddir, root, meas, fl.rstrip())
        if not os.path.isfile(fname):
            print >>sys.stderr, "File:", fname, "does not exist"
            sys.exit(1)
        datalist.append(fname)
    darkfile = datalist.pop(0)
    return darkfile, datalist

def gen3darkfile(ddir, root, meas):
    darkfile = "%s_%s_dark.fits" % (root, meas)
    darkfile = os.path.join(ddir, root, meas, darkfile)
    if not os.path.exists(darkfile):
        print >>sys.stderr, "File:", darkfile, "does not exist"
        sys.exit(1)
    return darkfile


def gen3datalist(ddir, root, meas, wavelengths, sources):
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
    return datalist

def gen3fit(cfgfilename):
    pl = bp.utils.ParamList()
    pl.read_from_file(cfgfilename)
    pl.read_from_file(pl.get_val('GEN3FIT_PREV_CFGFILE'))
    srcplate_cfg = pl.get_val('GEN3FIT_PREV_CFGFILE')

    ddir = pl.get_val('DATA_DIRECTORY')
    root = pl.get_val('STUDY_NAME')
    meas = pl.get_val('GEN3FIT_MEASUREMENT')
    sources = [int(v) for v in pl.get_val('SOURCES').split()]
    if len(sources) == 1 and sources[0] == 0:
        sources = np.arange(209)+1
    wavelengths = [int(v) for v in pl.get_val('WAVELENGTHS').split()]
    fmeas = pl.get_val('GEN3FIT_FLATFILEDIR')
    flatfile = pl.get_val('GEN3FIT_FLATFILENAME')
    flatfiledark = pl.get_val('GEN3FIT_FLATFILEDARKNAME')
    
    if pl.get_val('GEN3FIT_FILENAME_ENCODING') == 'GEN3':
        datalist = gen3datalist(ddir, root, meas, wavelengths, sources)
        darkfile = gen3darkfile(ddir, root, meas)
    elif pl.get_val('GEN3FIT_FILENAME_ENCODING') == 'LIST':
        filelistname = 'datalist-firstisdark.txt'
        if not os.path.exists(filelistname):
            print >>sys.stderr, "File", filelistname, "does not exist"
            sys.exit(1)
        darkfile, datalist = listdatalist(ddir, root, meas, filelistname)

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

    #bg3f.aphidc_flat_andro(datalist, darkfile, cfgfile, p, flatfile, flatfiledark)
    bg3f.aphidc_flat_andro(datalist, darkfile, srcplate_cfg, p)

if __name__ == '__main__':
    try:
        cfgfilename = sys.argv[1]
    except:
        sys.exit("Usage: "+os.path.basename(sys.argv[0])+" <cfgfilename>")
    if not os.path.exists(cfgfilename):
        sys.exit("ERROR: "+os.path.basename(sys.argv[0])+" file does not exist: "+cfgfilename)
    gen3fit(cfgfilename)

