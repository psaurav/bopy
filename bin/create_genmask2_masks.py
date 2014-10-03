#!/usr/bin/python

import os
import sys
import numpy as np
import bopy.utils.masks as masks

if __name__ == '__main__':
    
    try:
        fname = sys.argv[1]
        rmin = int(sys.argv[2])
        rmax = int(sys.argv[3])
        cmin = int(sys.argv[4])
        cmax = int(sys.argv[5])
        rpo = int(sys.argv[6])
        cpo = int(sys.argv[7])
        wls = sys.argv[8]
        ddir = sys.argv[9]
        refdir = sys.argv[10]
        trgdir = sys.argv[11]
        root = sys.argv[12]
        ref = sys.argv[13]
        trg = sys.argv[14]
        cutoff = float(sys.argv[15])
    except:
        sys.exit('Usage: {0} <fname.npy> <rmin> <rmax> \
                 <cmin> <cmax> <rpo> <cpo> <"wls"> \
                 <ddir> <refdir> <trgdir> <root> <ref> <trg> \
                 <cutoff>'.format(os.path.basename(sys.argv[0])))

    # Creating and saving genmask
    shape = np.load(fname).shape
    boundary = (rmin, rmax, cmin, cmax)
    righttop = (rpo, cpo)
    
    genmask = masks.genmask(shape, boundary, righttop)
    np.save('genmask2.npy', genmask)
    print 'saved genmask2.npy'

    # Creating and saving masks for each source
    if not os.path.exists('masks'):
        print 'masks does not exist'
        os.makedirs('masks')
    else:
        print 'masks exists'
        if os.path.isfile('masks'):
            print 'masks is file'
            os.remove('masks')
            os.makedirs('mask')

    print '  wls:', wls
    wls = [int(v) for v in wls.strip().split()]
    print '  wls:', wls
    #sys.exit()

    srcs = np.arange(209)+1
    for measdir, meas in zip([refdir, trgdir], [ref, trg]):
        print measdir, meas

    master = masks.create(shape)
    master = np.logical_not(master)
    for s in srcs:
        m = masks.create(shape)
        for w in wls:
            for measdir, meas in zip([refdir, trgdir], [ref, trg]):
            #for measdir, meas in zip([refdir, ], [ref, ]):
                fname = '{0}_{1}_wl{2}_s{3}_A.npy'.format(root, meas, w, s)
                fname = os.path.join(ddir, measdir, fname)
                if not os.path.isfile(fname):
                    print >>sys.stderr, 'File does not exist:', fname
                    sys.exit
                d = np.load(fname)
                m1 = masks.bkmaxcutoff_bottommasked(d, cutoff, toprbb=rpo)
                m = m | m1
        m = m | genmask
        ofile = 'mask_s{0}.npy'.format(s)
        np.save(os.path.join('masks', ofile), m)
        master = master & m
    ofile = 'mask_master.npy'
    np.save(os.path.join('masks', ofile), master)
