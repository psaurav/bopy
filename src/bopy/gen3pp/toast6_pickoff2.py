#!/usr/bin/python

import sys
import os
import numpy as np
import bopy as bp
import multiprocessing as mp
import numpy.ma as ma
#from bopy.utils.rfuncs import loess2d

# changes:
#       * ver 6
#           * no phase normalization
#       * loessed ahead of time, so, just read.
#       * bypass logdcorr etc.  if absent, create arrays of zero
#       * sname is not taken from ROOT

# This file will generate masks for the gim files produced by init-fftw2.py.
# The masks are defined as containing all points with an intensity of 
# greater than some pre-determined percentage of the maximum. 

def get_toast_phi((dfile, mfile, span, rpath, rfile)):
    """ Reads phi files, smooths and then returns the unmasked pixels
        In addition, performs normalization on the target files, if rfile is 
        not None
        
        Args:
        dfile - data file
        mfile - mask file
        span - span parameter for LOESS
        rpath - the path to ref files (None if not for target phase normalization)
        rfile - reference data file
    """
    d = np.load(dfile)
    m = np.load(mfile)
    if rfile is not None:
        from bopy.gen3pp.gen3pp_utils import phase_normalization
        r = np.load(os.path.join(rpath, rfile))
        #d = phase_normalization(d, r, m)
    #d = loess2d(d, span)
    d = ma.compressed(ma.masked_array(d, mask=m))
    return -d   # -d because toast phase decreases with r

def get_toast_lnX((dfile, mfile, span)):
    d = np.load(dfile)
    #d = loess2d(d, span)
    m = np.load(mfile)
    d = ma.compressed(ma.masked_array(d, mask=m))
    return np.log(d)

def write_toast_lnX(ddir, dmeas, meas, sname, w, srcs, dtype, fstring, 
                    msmaskdir, msmaskstr, span, corr, outputfile):
    args = ((os.path.join(ddir, dmeas, fstring.format(sname, meas, w, s, dtype)), 
            os.path.join(msmaskdir, msmaskstr.format(sname, s)),
            span) for s in srcs)
    pool = mp.Pool(mp.cpu_count() - 1)
    res = pool.map_async(get_toast_lnX, args)
    res = res.get()
    # get all the arrays in a single line
    res = np.hstack(res)
    res = res - corr
    bp.io.write_fem(outputfile, res)

def write_toast_phi(ddir, dmeas, meas, sname, w, srcs, dtype, fstring, 
                    msmaskdir, msmaskstr, span, corr, outputfile, rpath=None, rfile=None):
    if rpath is None:
        args = ((os.path.join(ddir, dmeas, fstring.format(sname, meas, w, s, dtype)), 
                os.path.join(msmaskdir, msmaskstr.format(sname, s)),
                span, rpath, rfile) for s in srcs)
    else:
        args = ((os.path.join(ddir, dmeas, fstring.format(sname, meas, w, s, dtype)), 
                os.path.join(msmaskdir, msmaskstr.format(sname, s)),
                span, rpath, rfile.format(w, s)) for s in srcs)
    pool = mp.Pool(mp.cpu_count() - 1)
    res = pool.map_async(get_toast_phi, args)
    res = res.get()
    # get all the arrays in a single line
    res = np.hstack(res)
    res = res - corr
    bp.io.write_fem(outputfile, res)

def write_qm_header(qmfile):
    print >>qmfile, "QM file 3D"
    print >>qmfile, "Dimension 3"
    print >>qmfile, ""

def sample_data(ddir, tdir, rdir, tmeas, rmeas, sname, s_nrow, s_ncol, wls, nrow, ncol, 
                 ssmaskdir, musp_corr, span, ro=0, co=0, rs=5, cs=5,
                 sro=0, sco=0, srs=1, scs=1, logdcorr=None, logacorr=None, phicorr=None):
    """ Creates the toast files
    Args:
        ddir: the data director
        tdir: the target subdirector (e.g LB)
        rdir: the reference subdirectory (e.g REF)
        tmeas: target measurement (e.g LB)
        rmeas: reference measurement (e.g REF)
        sname: study name (eg 20130410_G3033_Patient)
        #fstring: the file name convention (e.g {0}_{1}_wl{3}, s{4})
        #dext: the filename extension for data files (e.g npy or gim)
        srcs: list of source indexes
        wls: the list of wavelengths (int, nm)
        nrow: the number of rows of detector data matrix
        ncol: the number of cols of detector data matrix
        ssmaskdir: single spectra mask dir
        ro: row offset
        rs: row stride
        co: column offset
        cs: column stride
        musp_corr: musp for point source correction
    """

    print >>sys.stderr, 'ddir =', ddir
    print >>sys.stderr, 'tdir =', tdir
    print >>sys.stderr, 'rdir =', rdir
    print >>sys.stderr, 'tmeas =', tmeas
    print >>sys.stderr, 'rmeas =', rmeas
    print >>sys.stderr, 'sname =', sname

    # some variables
    df_string = '{0}_{1}_wl{2}_s{3}_{4}.npy'
    ssmask_string = 'mask_s{0}.npy'
    msmask_dir = 'masks_ms'
    msmask_string = '{0}_s{1}_msmask.npy'
    msmask_master_string = '{0}_msmask_master.npy'

    if not os.path.exists(msmask_dir):
        os.makedirs(msmask_dir)

    src_master_mask = np.ones((s_nrow, s_ncol), dtype=bool)  # all True


    # create regular grid sample mask for sources
    (r, c) = np.indices(src_master_mask.shape)
    ir = r % srs != sro   # row band
    ic = c % scs != sco   # column band
    g_src = ir | ic         # intersection of row and column bands
    index_src = r*s_ncol + c
    srcs = ma.compressed(ma.masked_array(index_src, mask=g_src))
    srcs = srcs + 1     # unit-offset indexing

    print "srcs:", srcs

    # Init mask
    master_mask = np.ones((nrow, ncol), dtype=bool)  # all True
    # create regular grid sample mask
    (r, c) = np.indices(master_mask.shape)
    ir = r % rs != ro   # row band
    ic = c % cs != co   # column band
    g = ir | ic         # intersection of row and column bands

    # Loop over sources
    for s in srcs:
        # create common mask
        src_mask = np.zeros((nrow, ncol), dtype=bool) # all False
        for w in wls:
            mask = np.load(os.path.join(ssmaskdir, ssmask_string.format(s)))
            src_mask |= mask 
        
        # pick out only those on the sample grid and save it to disk
        src_mask |= g
        np.save(os.path.join(msmask_dir, msmask_string.format(sname, s)), src_mask) 
        # save them on the master mask
        master_mask &= src_mask
    # save master grid on disk
    np.save(os.path.join(msmask_dir, msmask_master_string.format(sname)), master_mask)

    ## DEBUG ###
    #sys.exit()

    # the detector indices
    oldindex = r*ncol + c
    # create a 1D vector of indices of the non-masked elements
    unmasked_indices = ma.compressed(ma.masked_array(oldindex, mask=master_mask))
    # create a 1D vector of new continuous indices for the non-masked elements
    newindex = np.arange(len(unmasked_indices))
    # create a map from the old index to new index
    old2new = dict(zip(unmasked_indices, newindex))

    print >>sys.stderr, 'nDetectors:', newindex.shape

    #
    # QM file portion
    #
    
    # sample source positions
    spos = bp.io.read_sdpos(ddir, 'sources.txt')    
    srcs_idx = srcs - 1
    spos = spos[srcs_idx]
    # correct sy for point source
    spos[:,1] = spos[:,1] + 1.0/musp_corr

    # sample detector positions
    dpos = bp.io.read_sdpos(ddir, 'detectors.txt')    
    dpos = dpos[unmasked_indices]
    
    # create toast data directory
    toastddir = 'toastdata_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(sro, sco, srs, scs, ro, co, rs, cs)
    if not os.path.exists(toastddir):
        os.makedirs(toastddir)
    
    # QM file
    qmfile = os.path.join(toastddir, 
                          '{0}_{1}_{2}_{3}_{4}.qm'.format(sname, ro, co, rs, cs))
    qmfile = open(qmfile, 'w')
    
    # write QM file header
    write_qm_header(qmfile)
    
    # write source positions
    print >>qmfile, "SourceList", len(spos), "fixed"
    np.savetxt(qmfile, spos)
    print >>qmfile, ""

    # write detector positions
    print >>qmfile, "MeasurementList", len(dpos)
    np.savetxt(qmfile, dpos)
    print >>qmfile, ""

    # write linked list
    ll = []
    print >>qmfile, "LinkList"
    for s in srcs:
        m = np.load(os.path.join(msmask_dir, msmask_string.format(sname, s)))
        sdets = ma.compressed(ma.masked_array(oldindex, mask=m))
        snew = np.array([old2new[i] for i in sdets])
        n = len(snew)
        print >>qmfile, '%d:'%(n), ' '.join(str(i) for i in snew)
    qmfile.close()

    #
    # Data files portion
    # 
    # loop over
    femstring = '{0}_{1}_{2}.fem'
    dtype = {'A': 'lnA', 'D': 'lnD', 'phi': 'phi'}
    mtype = {tmeas: 'trg', rmeas: 'ref'}
    fem_str = '{0}_{1}_{2}.fem'
    for w, dcorr, acorr, pcorr in zip(wls, logdcorr, logacorr, phicorr):
        print >>sys.stderr, w
        #for mt, mname in zip(mtype):
        # D trg
        print >>sys.stderr, '  D...',
        outputfile = os.path.join(toastddir, fem_str.format('trg', 'lnD', w))
        write_toast_lnX(ddir, tdir, tmeas, sname, w, srcs, 'DC', 
                            df_string, msmask_dir, msmask_string, span, dcorr, outputfile)
        # D ref
        outputfile = os.path.join(toastddir, fem_str.format('ref', 'lnD', w))
        write_toast_lnX(ddir, rdir, rmeas, sname, w, srcs, 'DC', 
                            df_string, msmask_dir, msmask_string, span, 0., outputfile)
        print >>sys.stderr, 'done'

        # A trg
        print >>sys.stderr, '  A...',
        outputfile = os.path.join(toastddir, fem_str.format('trg', 'lnA', w))
        write_toast_lnX(ddir, tdir, tmeas, sname, w, srcs, 'A', 
                            df_string, msmask_dir, msmask_string, span, acorr, outputfile)
        # A ref
        outputfile = os.path.join(toastddir, fem_str.format('ref', 'lnA', w))
        write_toast_lnX(ddir, rdir, rmeas, sname, w, srcs, 'A', 
                            df_string, msmask_dir, msmask_string, span, 0., outputfile)
        print >>sys.stderr, 'done'

        # phi ref
        outputfile = os.path.join(toastddir, fem_str.format('ref', 'phi', w))
        write_toast_phi(ddir, rdir, rmeas, sname, w, srcs, 'phi', 
                            df_string, msmask_dir, msmask_string, span, 0., outputfile)
        # phi trg
        print >>sys.stderr, '  phi...',
        outputfile = os.path.join(toastddir, fem_str.format('trg', 'phi', w))
        rpath = os.path.join(ddir, rdir)  
        rfile = sname + '_' + rmeas + '_wl{0}_s{1}_phi.npy'
        write_toast_phi(ddir, tdir, tmeas, sname, w, srcs, 'phi', 
                            df_string, msmask_dir, msmask_string, span, pcorr, outputfile, 
                            rpath, rfile)
        print >>sys.stderr, 'done'

def driver_toast_data(*args):
    pl = bp.utils.ParamList()
    for a in args:
        pl.read_from_file(a)

    # sources
    s_nrow = int(pl.get_val('SRC_NROW'))
    s_ncol = int(pl.get_val('SRC_NCOL'))
    sro = int(pl.get_val('SRC_ROW_OFFSET'))
    sco = int(pl.get_val('SRC_COL_OFFSET'))
    srs = int(pl.get_val('SRC_ROW_STRIDE'))
    scs = int(pl.get_val('SRC_COL_STRIDE'))

    # Data
    ddir = pl.get_val('PP_DATA_DIR')
    tdir = pl.get_val('TRG_DIR_NAME')
    rdir = pl.get_val('REF_DIR_NAME')
    tmeas = pl.get_val('TRG_NAME')
    rmeas = pl.get_val('REF_NAME') 
    sname = pl.get_val('ROOT')
    wls = np.array([int(w) for w in pl.get_val('WAVELENGTHS').split()])
    nrow = int(pl.get_val('RF_DATA_NROW'))
    ncol = int(pl.get_val('RF_DATA_NCOL')) 
    ssmaskdir = pl.get_val('SINGLE_SPECT_MASK_DIR')
    ro = int(pl.get_val('DATA_ROW_OFFSET'))
    co = int(pl.get_val('DATA_COL_OFFSET'))
    rs = int(pl.get_val('DATA_ROW_STRIDE'))
    cs = int(pl.get_val('DATA_COL_STRIDE'))
    if pl.get_val('LOGDCORR') is not None:
        logdcorr = np.array([float(w) for w in pl.get_val('LOGDCORR').split()])
        logacorr = np.array([float(w) for w in pl.get_val('LOGACORR').split()])
        phicorr = np.array([float(w) for w in pl.get_val('PHASECORR').split()])
    else:
        logdcorr = np.array([0. for w in wls])
        logacorr = np.array([0. for w in wls])
        phicorr = np.array([0. for w in wls])
    
    musp_corr = float(pl.get_val('MUSP_CORR_SOURCE_POSITION'))
    span = float(pl.get_val('LOESS_SPAN'))
    sample_data(ddir, tdir, rdir, tmeas, rmeas, sname, s_nrow, s_ncol, wls, nrow, ncol, 
                ssmaskdir, musp_corr, span, ro=ro, co=ro, rs=rs, cs=cs,
                sro=sro, sco=sco, srs=srs, scs=scs,
                logdcorr=logdcorr, logacorr=logacorr, phicorr=phicorr)
    
