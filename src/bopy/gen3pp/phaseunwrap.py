import os
import numpy as np
import numpy.ma as ma

def trgvsref(troot, rdir, rroot, mdir, mroot, wls, srcs, phasetag='median'):

    if phasetag is 'mean':
        print 'Phasetag:', phasetag
        phasetag = np.mean
    else:
        print 'Phasetag:', phasetag
        phasetag = np.median
        
    for s in srcs:
        mname = os.path.join(mdir, '{0}_s{1}.npy'.format(mroot, s))
        m = np.load(mname)
        for w in wls:
            tname = '{0}_wl{1}_s{2}_phi.npy'.format(troot, w, s)
            t = np.load(tname)
            tm = ma.masked_array(t, mask=m)
            rname = os.path.join(rdir, '{0}_wl{1}_s{2}_phi.npy'.format(rroot, w, s))
            r = np.load(rname)
            rm = ma.masked_array(r, mask=m)
            d = phasetag(tm-rm)
            if d > np.pi:
                print 'w:', w, ' s:', s, '-2PI'
                t = t - 2.0*np.pi
                os.rename(tname, '{0}_wrapped'.format(tname))
                np.save(tname, t)
            elif d < -np.pi: 
                print 'w:', w, ' s:', s, '+2PI'
                t = t + 2.0*np.pi
                os.rename(tname, '{0}_wrapped'.format(tname))
                np.save(tname, t)
