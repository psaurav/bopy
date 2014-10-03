import sys
import numpy as np

def data_mean(d, axis=0):
    """
    """
    d = np.mean(d, dtype=np.float64, axis=axis)
    return d.astype(float)

def data_reorient(d, p):
    """
    """
    print >>sys.stderr, "Reorienting...",
    if (p.corr_flipud == 1):
        d = np.flipud(d)
    if (p.corr_fliplr == 1):
        d = np.fliplr(d)
    while (p.corr_rotation > 0):
        d = np.rot90(d)
        corr_rotation = corr_rotation - 1
    print >>sys.stderr, "done."
    return d

def data_reorient2(d, pl):
    """
    """
    print >>sys.stderr, "Reorienting...",
    lr = int(pl.get_val('REORIENT_CCD_FLIPLR'))
    rot90 = int(pl.get_val('REORIENT_CCD_ROT90'))
    if rot90 > 0:
        d = np.rot90(d, rot90)
    if lr > 0:
        d = np.fliplr(d)
    print >>sys.stderr, "done."
    return d

def data_reorient_udrot(d, pl, ud=None, rot90=None):
    """
    """
    print >>sys.stderr, "Reorienting...",
    if (not ud is None) and (not rot90 is None):
        if ud > 0:
            d = np.flipud(d)
        if rot90 > 0:
            d = np.rot90(d, rot90)
        print >>sys.stderr, "done."
        return d
        
    ud = pl.get_val('REORIENT_CCD_FLIPUD')
    if ud == None:
        ud = pl.get_val('CORRECTION_FLIPUD')
    ud = int(ud)
    rot90 = pl.get_val('REORIENT_CCD_ROT90')
    if rot90 == None:
        rot90 = pl.get_val('CORRECTION_ROTATION')
    rot90 = int(rot90)
    if ud > 0:
        d = np.flipud(d)
    if rot90 > 0:
        d = np.rot90(d, rot90)
    print >>sys.stderr, "done."
    return d
