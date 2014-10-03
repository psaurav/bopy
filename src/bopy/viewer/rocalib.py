import os
import sys
import bopy as bp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

def _init_p():
    p = bp.utils.struct()
    #
    # Set default parameters
    #
    p.rebin = 1
    p.rot90 = 0
    p.lr = -1
    p.dpixel = 0.
    p.x0 = 0
    p.z0 = 0
    p.cos_theta = 0.
    p.sin_theta = 0. 
    p.lastx = 0
    p.lasty = 0
    p.currx = 0
    p.curry = 0

    ncols = bp.gen3fit.Gen3.ncols
    sspacings = bp.gen3fit.Gen3.sspacings
    p.dist = float((ncols-1)*sspacings) 

    return p

def _write2file(p):
    """
    """
    global pl
    global cfgname
    
    #
    # write cfg file
    #
    pl.set_val('SOURCE_NCOLUMNS', bp.gen3fit.Gen3.ncols)
    pl.set_val('SOURCE_NROWS', bp.gen3fit.Gen3.nrows)
    pl.set_val('SOURCE_SPACINGS', bp.gen3fit.Gen3.sspacings)
    pl.set_val('TANK_THICKNESS', p.tank_thickness)
    pl.set_val('REBIN', p.rebin)
    pl.set_val('CALIBRATED_CCD_X0', p.x0)
    pl.set_val('CALIBRATED_CCD_Z0', p.z0)
    pl.set_val('CALIBRATED_CCD_DPIXEL', p.dpixel)
    pl.set_val('CALIBRATED_CCD_COS_THETA', p.cos_theta)
    pl.set_val('CALIBRATED_CCD_SIN_THETA', p.sin_theta)
    pl.set_val('REORIENT_CCD_FLIPLR', p.lr)
    pl.set_val('REORIENT_CCD_ROT90', p.rot90)

    pl.write_to_file(cfgname) 

    #
    # write sources file
    #
    ncol = bp.gen3fit.Gen3.ncols
    nrow = bp.gen3fit.Gen3.nrows
    delS = bp.gen3fit.Gen3.sspacings

    c0 = (ncol - 1) / 2
    r0 = (nrow - 1) / 2
    count = 1
    srcfile = open('sources.txt', 'w')
    for r in xrange(0, nrow):
        for c in xrange(0, ncol):
            x = delS*(c-c0)

            y = 0.0
            z = delS*(r-r0)
            line = '%d %d %d %f %f %f'%(count, c, r, x, y, z)
            srcfile.write(line + '\n')
            count = count + 1
    srcfile.close()
    print >>sys.stderr, "Wrote sources.txt"

    #
    # write detector file
    #
    x0 = float(pl.get_val("CALIBRATED_CCD_X0"))
    z0 = float(pl.get_val("CALIBRATED_CCD_Z0"))
    L = float(pl.get_val("TANK_THICKNESS"))
    dpixel = float(pl.get_val("CALIBRATED_CCD_DPIXEL"))
    ct = float(pl.get_val("CALIBRATED_CCD_COS_THETA"))
    st = float(pl.get_val("CALIBRATED_CCD_SIN_THETA"))
    nrow = int(pl.get_val("DATA_CCD_NROWS"))
    ncol = int(pl.get_val("DATA_CCD_NCOLS"))
    count = 1
    detfile = open('detectors.txt', 'w')
    for r in xrange(0, nrow):
        for c in xrange(0, ncol):
            # translate
            z = r-z0
            x = c-x0
            # rotate
            zp =  z*ct + x*st
            xp = -z*st + x*ct
            # scale
            zp = dpixel*zp
            xp = dpixel*xp
            y = L
            line = '%d %d %d %f %f %f'%(count, c, r, xp, y, zp)
            detfile.write(line + '\n')
            count = count + 1
    detfile.close()
    print >>sys.stderr, "Wrote detectors.txt"
    print >>sys.stderr, "Done"    

def _onclick(event):
    #print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
    #    event.button, event.x, event.y, event.xdata, event.ydata)
    global d
    global p
    if (event.button == 3):
        p.lastx = p.currx
        p.lasty = p.curry
        p.currx = float(int(event.xdata + 0.5))
        p.curry = float(int(event.ydata + 0.5))
        if (p.currx < p.lastx): # awkward--should just be abs(lasty - curry) or something
            delx = p.lastx-p.currx
            dely = p.lasty-p.curry
        else:
            delx = p.currx-p.lastx
            dely = p.curry-p.lasty
        d = math.sqrt(delx**2 + dely**2)
        p.cos_theta = delx/d
        # sin theta is negative because 
        # z increases down
        p.sin_theta = -dely/d
    
        p.dpixel = float(p.dist)/d
        print >>sys.stderr, "lastx, lasty, currx, curry:", p.lastx, p.lasty, p.currx, p.curry
        print >>sys.stderr, "cos_theta:", p.cos_theta
        print >>sys.stderr, "sin_theta:", p.sin_theta
        print >>sys.stderr, "dpixel:", p.dpixel

    if (event.button == 2):
        p.x0 = float(int(event.xdata + 0.5))
        p.z0 = float(int(event.ydata + 0.5))
        print >>sys.stderr, "x0:", p.x0, "  z0:", p.z0
        

def _onkey(event):
    global p
    global d
    print >>sys.stderr, event.key
    if (event.key == 'w'):
        _write2file(p)
        return
    if (event.key == 'l'):
        d = np.fliplr(d)
        p.lr = -1*p.lr
        plt.imshow(d, interpolation='nearest', cmap=cm.Greys_r)
    if (event.key == 'r'):
        if p.lr > 0:    # no rotation allowed for flipped data
            print >>sys.stderr, "WARNING: Cannot rotate after flips.  Flip back to rotate again."
            return
        d = np.rot90(d)
        p.rot90 = p.rot90 + 1
        if (p.rot90 >= 4):
            p.rot90 = p.rot90 - 4
        plt.imshow(d, interpolation='nearest', cmap=cm.Greys_r)
    plt.draw()

def rocalib(filename, thickness=60.0, rebin=1):
    """
    """
    global d
    d = bp.io.read_frame(filename)
    # take the first frame
    if len(d.shape) == 3:
        d = d[0]
    d = np.array(d)
    print d.shape
    #
    # set up parameter
    #
    global p
    p = _init_p()
    p.tank_thickness = thickness
    #
    # Rebin data if needed
    #
    if rebin > 1:
        d = bp.utils.rebindata2d(d, rebin)
        p.rebin = rebin
    global cfgname
    global pl
    cfgname = "%s_%dx%d.cfg" % (os.path.splitext(os.path.basename(filename))[0], rebin, rebin)
    pl = bp.utils.ParamList()
    pl.set_val('DATA_CCD_NROWS', np.shape(d)[0])
    pl.set_val('DATA_CCD_NCOLS', np.shape(d)[1])
    #
    # Display
    #
    global fig
    global ax
    fig, ax = plt.subplots()
    plt.imshow(d, interpolation='nearest', cmap=cm.Greys_r)
    #plt.axis('off')
    #fig = plt.figure(1)
    #fig.axes.set_yscale('linear')
    ax.set_yscale('linear')
    ax.set_xscale('linear')


    kid = fig.canvas.mpl_connect('key_press_event', _onkey)
    cid = fig.canvas.mpl_connect('button_press_event', _onclick)

    #plt.show()
    plt.show()
