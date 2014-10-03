import os
import numpy as np
import numpy.ma as ma
import bopy as bp
from bopy.utils.rfuncs import loess2d
from scipy.ndimage.filters import gaussian_filter

def save_loessed_signal((infilename, ofilename, span)):
    d = bp.io.read_frame(infilename)
    d = loess2d(d, span=span)
    np.save(ofilename, d)
    return

def save_gaussed_signal((infilename, ofilename, sigma)):
    d = bp.io.read_frame(infilename)
    d = gaussian_filter(d, sigma)
    np.save(ofilename, d)
    return

def mp_save_loessed_signal(nprocs, infilenames, span):
    from multiprocessing import Pool
    ofilenames = map(os.path.basename, infilenames) 
    for i, o in zip(infilenames, ofilenames):
        print i, o
    args = ((infile, ofile, span) for infile, ofile in zip(infilenames, ofilenames))
    pool = Pool(processes=nprocs)
    result = pool.map_async(save_loessed_signal, args)
    return len(result.get())

def sp_save_loessed_signal(nprocs, infilenames, span):
    ofilenames = map(os.path.basename, infilenames) 
    for infile, ofile in zip(infilenames, ofilenames):
        print infile, ofile
        args = (infile, ofile, span)
        save_loessed_signal(args)
    return 0

def sp_save_gaussed_signal(nprocs, infilenames, sigma):
    ofilenames = map(os.path.basename, infilenames) 
    for infile, ofile in zip(infilenames, ofilenames):
        #print infile, ofile
        args = (infile, ofile, sigma)
        save_gaussed_signal(args)
    return 0


def get_signal_mask(ddir, meas1, meas2, studyname, filestring, s, w, cutoff, ptype='A'):
    r = bp.io.read_ppfile(ddir, meas1, studyname, filestring, s, w, ptype)
    t = bp.io.read_ppfile(ddir, meas2, studyname, filestring, s, w, ptype)
    rm = bp.utils.mask_maxcutoff_gt(r, cutoff)
    tm = bp.utils.mask_maxcutoff_gt(t, cutoff)
    return rm & tm

def get_loessed_signal_mask2(ddir, meas1dir, meas2dir, meas1, meas2, studyname, filestring, s, w, cutoff, toprowsmask=0, span=0.02, ptype='A'):
    """Returns the mask, given two amplitude images.  
    """
    r = bp.io.read_ppfile2(ddir, meas1dir, meas1, studyname, filestring, s, w, ptype)
    t = bp.io.read_ppfile2(ddir, meas2dir, meas2, studyname, filestring, s, w, ptype)
    if not span == 0.0:
        r = loess2d(r, span=span)
        t = loess2d(t, span=span)
    rm = bp.utils.mask_maxcutoff(r, cutoff)
    tm = bp.utils.mask_maxcutoff(t, cutoff)
    if toprowsmask > 0:
        return bp.utils.mask_toprows(rm | tm, toprowsmask)
    else:
        return rm | tm

def save_loessed_signal_mask2((ddir, meas1dir, meas2dir, meas1, meas2, studyname, filestring, s, w, cutoff, toprowsmask, span, ptype)):
    mask = get_loessed_signal_mask2(ddir, meas1dir, meas2dir, meas1, meas2, studyname, filestring, s, w, cutoff, toprowsmask, span, ptype)
    dirname = 'masks'
    fname = '{0}_wl{1}_s{2}_mask.npy'.format(studyname, w, s)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    f = os.path.join(dirname, fname)
    np.save(f, mask)
    return 1

def mp_save_loessed_signal_mask2(nproc, ddir, meas1dir, meas2dir, meas1, meas2, studyname, filestring, srcs, wls, cutoff, toprowsmask, span, ptype):
    from multiprocessing import Pool
    args = ((ddir, meas1dir, meas2dir, meas1, meas2, studyname, filestring, s, w, cutoff, toprowsmask, span, ptype) for s in srcs for w in wls)
    pool = Pool(processes=nproc)
    result = pool.map_async(save_loessed_signal_mask2, args)
    return len(result.get())

def driver_mp_save_loessed_signal_mask2(pl):
    ddir = '.'
    meas1dir = pl.get_val('TRG_NAME')
    meas1 = pl.get_val('TRG_NAME')
    meas2dir = pl.get_val('REF_NAME')
    meas2 = pl.get_val('REF_NAME')
    studyname = pl.get_val('STUDY_NAME')
    filestring = pl.get_val('FILESTRING')
    srcs = [int(s) for s in pl.get_val('SOURCES').split()]
    wls = [int(w) for w in pl.get_val('WAVELENGTHS').split()]
    #
    # new
    #
    cutoff = float(pl.get_val('SIGNAL_MAX_CUTOFF'))
    toprowsmask = int(pl.get_val('TOPROWSMASK'))
    span = float(pl.get_val('LOESS_SPAN'))
    ptype = 'A'
    res = mp_save_loessed_signal_mask2(nproc, ddir, meas1dir, meas2dir, 
                                       meas1, meas2, studyname, filestring, 
                                       srcs, wls, cutoff, toprowsmask, span, ptype)
    
    

def masked_phase_diff_mean(pt, pr, mask):
    return np.mean(ma.masked_array(pt - pr, mask=mask)) 

def phase_normalization(pt, pr, mask):
    diff_mean = masked_phase_diff_mean(pt, pr, mask)
    if diff_mean > np.pi:
        return pt - 2.0*np.pi
    elif diff_mean < -np.pi:
        return pt + 2.0*np.pi
    else:
        return pt 

def loessmin_diff((file1, file2, filex1, filex2, cutoff, toprowsmask, span)):
    from bopy.io import read_frame
    from bopy.utils import mask2_maxcutoff 
    from bopy.utils.rfuncs import loess2d
    d1 = read_frame(file1)
    d2 = read_frame(file2)
    dmask = mask2_maxcutoff(read_frame(filex1),
            read_frame(filex2),
            cutoff, toprowsmask=toprowsmask)
    d1 = loess2d(d1, span=span)
    d1 = np.min(ma.compressed(ma.masked_array(d1, mask=dmask)))
    d2 = loess2d(d2, span=span)
    d2 = np.min(ma.compressed(ma.masked_array(d2, mask=dmask)))
    return d1 - d2

def loess_diff((file1, file2, filex1, filex2, cutoff, toprowsmask, span)):
    from bopy.io import read_frame
    from bopy.utils import mask2_maxcutoff 
    from bopy.utils.rfuncs import loess2d
    d1 = read_frame(file1)
    d2 = read_frame(file2)
    dmask = mask2_maxcutoff(read_frame(filex1),
            read_frame(filex2),
            cutoff, toprowsmask=toprowsmask)
    d1 = loess2d(d1, span=span)
    d2 = loess2d(d2, span=span)
    return ma.masked_array(d1 - d2, mask=dmask)

def loess_logdiv((file1, file2, filex1, filex2, cutoff, toprowsmask, span)):
    from bopy.io import read_frame
    from bopy.utils import mask2_maxcutoff 
    from bopy.utils.rfuncs import loess2d
    d1 = read_frame(file1)
    d2 = read_frame(file2)
    dmask = mask2_maxcutoff(read_frame(filex1),
            read_frame(filex2),
            cutoff, toprowsmask=toprowsmask)
    d1 = loess2d(d1, span=span)
    d2 = loess2d(d2, span=span)
    return ma.masked_array(np.log(d1/d2), mask=dmask)

def togauss(meas, sigma=2.):
    import glob
    nprocs = 8
    dirname = '{0}-gsmooth-{1}'.format(meas, sigma)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    files = glob.glob('../{0}/*[A|DC|phi].npy'.format(meas))
    files = sorted(files)
    f = open('togauss.log', 'w')
    print >>f, 'MEAS =', meas
    print >>f, 'GAUSS_SIGMA =', sigma
    f.close()
    sp_save_gaussed_signal(nprocs, files, sigma)
    os.chdir('..')
