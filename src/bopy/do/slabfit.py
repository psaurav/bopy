import numpy as np
import numpy.ma as ma

def read_muamusp(wls, dims=(11,19)):
    mua = {}
    musp = {}
    for w in wls:
        res = np.loadtxt('results_{0}.txt'.format(w))
        mua[w] = res[:,3].reshape(dims)
        musp[w] = res[:,4].reshape(dims)
    return mua, musp

def masks(wls, mu, mn=0.0, mx=0.05, dims=(11,19)):
    m = {}
    for w in wls:
        m[w] = (mu[w] < mn) | (mu[w] > mx)
    return m

def get_mean(wls, mu, mask, (rn, rx, cn, cx)):
    mean = []
    for w in wls:
        x = ma.mean(ma.masked_array(mu[w], mask=mask[w])[rn:rx, cn:cx])
        mean.append(x)
    return np.array(mean)
