import numpy as np
import multiprocessing as mp
from bopy.utils import rc_bkmax
from bopy.gen3fit.pixelfit import pixel_cosfit3, get_wt, pixel_fftfit


def gen3fit_simple(data, dt, freq, phi_init=None, r=None, c=None):

    if phi_init is None:
        if r is None or c is None:
            r, c = rc_bkmax(data)
        s = data[:, r, c]
        d, a, phi_init = pixel_fftfit(s, dt, freq)

    shape = data.shape
    rs = shape[1:3]
    wt = get_wt(len(s), dt, freq)

    rvec, cvec = np.indices(data[0,:,:].shape)
    rvec = rvec.ravel()
    cvec = cvec.ravel()
    args = ((rr*shape[2]+cc, data[:, rr, cc], wt, phi_init) for rr, cc in zip(rvec, cvec))
    pool = mp.Pool(mp.cpu_count() - 1)
    res = pool.map_async(pixel_cosfit3, args)
    res = res.get()
    pool.close()
    pool.join()
    res = np.array(res)
    return res[:,1].reshape(rs), res[:,2].reshape(rs), res[:,3].reshape(rs) 
