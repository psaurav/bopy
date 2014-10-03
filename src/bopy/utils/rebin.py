import sys

def rebin2d(a, shape):
    """
    Rebin 2D numpy array according to shape

    Keyword parameters:
        a -- numpy array
        shape -- new shape
    """
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def rebin3d(a, shape):
    """
    Rebin 3D numpy array according to shape

    Parameters:
        a -- numpy array
        shape -- new shape
    """
    sh = shape[0],shape[1],a.shape[1]//shape[1],shape[2],a.shape[2]//shape[2]
    return a.reshape(sh).mean(-1).mean(2)

def rebindata2d(d, n):
    """
    """
    print >>sys.stderr, "shape(d2d):", d.shape
    remainder = d.shape[0] % n
    if remainder == 0:
        nr = d.shape[0]
    else:
        nr = d.shape[0] - remainder

    remainder = d.shape[1] % n
    if remainder == 0:
        nc = d.shape[1]
    else:
        nc = d.shape[1] - remainder
    print >>sys.stderr, "new: nr =", nr, "nc =", nc
    d = d[0:nr,0:nc]
    print >>sys.stderr, " new shape(d2d):", d.shape
    nr = nr/n
    nc = nc/n
    return rebin2d(d, (nr,nc))

def rebindata3d(d, n):
    """
    """
    print >>sys.stderr, "shape(d3d):", d.shape
    nf = d.shape[0]
    remainder = d.shape[1] % n
    if remainder == 0:
        nr = d.shape[1]
    else:
        nr = d.shape[1] - remainder

    remainder = d.shape[2] % n
    if remainder == 0:
        nc = d.shape[2]
    else:
        nc = d.shape[2] - remainder
    print >>sys.stderr, "new: nr =", nr, "nc =", nc
    d = d[:,0:nr,0:nc]
    print >>sys.stderr, " new shape(d3d):", d.shape
    nr = nr/n
    nc = nc/n
    return rebin3d(d, (nf,nr,nc))
