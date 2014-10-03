from bopy.utils import ParamList, struct

def get_calbsrc_params(plfname):
    pl = ParamList(plfname)
    p = struct()
    p.wls = [int(w) for w in pl.get_val('WAVELENGTHS').strip().split()]
    p.sigma = float(pl.get_val('SIGMA'))
    p.calib_pixels = (int(s) for s in pl.get_val('CALIB_PIXELS').strip().split())
    p.dkinds = [v for v in pl.get_val('DKINDS').strip().split()]
    p.calibsrcs = np.array([int(v) for v in pl.get_val('CALIB_SRCS').strip().split()])

    return p

def create_calibssrc(p):
    files = {}
    for X in p.dkinds:
        files[X] = open('srccalib-{0}.txt'.format(X), 'w')

    for ddir, meas in zip(p.ddirs, p.exs):
        for w in p.wls:
            for s in p.calibsrcs:
                for X in p.dkinds:
                    cval = scalib.get_calibration_singlepixel(ddir, root, meas, sigma, w, X, s, calibpixels)
                    print >>files[X], meas, w, s, cval

    for X in dkinds:
        files[X].close()

    srcs = np.arange(209)+1
    f = open('src-calibsrc.txt', 'w')
    for s in srcs:
        print >>f, s, calibsrcs[np.argmin(np.abs(s - calibsrcs))]
    f.close()

