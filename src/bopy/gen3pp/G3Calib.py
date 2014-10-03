#!/usr/bin/python

import sys
import os
import math
import numpy as np
import pylab 
import readline
import bopy as bp
from bopy.utils import struct, ParamList

class G3Calib():

    def __init__(self, initfilename, outputfilename):

        self.pl = ParamList(initfilename)
        self.outputfilename = outputfilename
        self.guide()
        self.initvals()

    def initvals(self):
        self.flipUD = -1
        self.rotated = 0
        self.currx = 0
        self.curry = 0
 
        self.src_ncols = int(self.pl.get_val("SOURCE_NCOLUMNS"))
        self.src_col_length = float(self.pl.get_val("SOURCE_SPACING"))
        self.dist = (self.src_ncols -1)*self.src_col_length
        
    def guide(self):
        print >>sys.stderr, " u                 : flip up-down"
        print >>sys.stderr, " r                 : rotate anti-clockwise 90-degrees"
        print >>sys.stderr, " w                 : write to file"
        print >>sys.stderr, " button-2 (middle) : set the origin"
        print >>sys.stderr, " button-2 (right)  : set point for calibration"
        print >>sys.stderr, " w                 : write to file"
        print >>sys.stderr, " s                 : show sources"
        return

    def set_data(self, filename):
        # fits file
        self.data = bp.io.read_fits(filename)
        # record the data dimension
        #nz, ny, nx = np.shape(scidata)
        nd = np.shape(self.data)
        if len(nd) == 3:
            self.data = self.data[0]
        print >>sys.stderr, "max-min:", np.max(self.data), np.min(self.data)
        self.pl.set_val('ORIG_RF_DATA_NCOL', nd[-1])
        self.pl.set_val('ORIG_RF_DATA_NROW', nd[-2])
        self.pl.set_val('REBIN', 1)
        print >>sys.stderr, "nx:", nd[-1], "  ny:", nd[-2]
        return

    def updown(self):
        print >>sys.stderr, "Pressed u"
        if (self.rotated != 0):
            print >>sys.stderr, "Cannot perform flipUD: rotated = ", self.rotated
            return
        self.data = np.flipud(self.data)
        self.flipUD = -1*self.flipUD
        print >>sys.stderr, "Current flipUD: ", self.flipUD

    def rotate(self):
        print >>sys.stderr, "Pressed r"
        self.data = np.rot90(self.data)
        self.rotated = self.rotated + 1
        if (self.rotated >= 4):
            self.rotated = self.rotated - 4
        print >>sys.stderr, "Current rotation: ", self.rotated

    def set_x0z0(self, x0, z0):
        self.x0 = x0
        self.z0 = z0
        print >>sys.stderr, "x0:", x0, "z0:", z0

    def set_xz(self, x, z):
        
        self.lastx = self.currx
        self.lasty = self.curry
        self.currx = x
        self.curry = z
        if (self.currx < self.lastx): # awkward--should just be abs(lasty - curry) or something
            delx = self.lastx - self.currx
            dely = self.lasty - self.curry
        else:
            delx = self.currx - self.lastx
            dely = self.curry - self.lasty
        
        d = math.sqrt(delx**2 + dely**2)
        self.cos_theta = delx/d
        # sin theta is negative because 
        # z increases down
        self.sin_theta = -dely/d
        self.dpixel = float(self.dist)/d

        print 'd=%f, last=(%f,%f), current=(%f,%f)'%(d, self.lastx, self.lasty, self.currx, self.curry)
        print '  dpixel =', self.dpixel
        print '  costheta =', self.cos_theta
        print '  sintheta =', self.sin_theta

    def write(self):

        # Orientation
        self.pl.set_val('CORRECTION_ROTATION', self.rotated)
        self.pl.set_val('CORRECTION_FLIPUD', self.flipUD)

        # Calibration
        self.pl.set_val('CALIBRATED_CCD_DPIXEL', self.dpixel)
        self.pl.set_val('CALIBRATED_CCD_COS_THETA', self.cos_theta)
        self.pl.set_val('CALIBRATED_CCD_SIN_THETA', self.sin_theta)
        self.pl.set_val('CALIBRATED_CCD_X0', self.x0)
        self.pl.set_val('CALIBRATED_CCD_Z0', self.z0)

        nd = np.shape(self.data)
        self.pl.set_val('RF_DATA_NCOL', nd[-1])
        self.pl.set_val('RF_DATA_NROW', nd[-2])
         
        self.pl.write_to_file(self.outputfilename)

        # sources
        ncol = int(self.pl.get_val("SOURCE_NCOLUMNS"))
        nrow = int(self.pl.get_val("SOURCE_NROWS"))
        delS = float(self.pl.get_val("SOURCE_SPACING"))
        c0 = (ncol - 1) / 2
        r0 = (nrow - 1) / 2
        #count is the source index used by Han
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
        # write the detector file
        #
        x0 = self.x0
        z0 = self.z0
        L = float(self.pl.get_val("BREAST_TANK_THICKNESS_MM"))
        ct = self.cos_theta
        st = self.sin_theta
        nrow, ncol = np.shape(self.data)
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
                zp = self.dpixel*zp
                xp = self.dpixel*xp
                y = L
                line = '%d %d %d %f %f %f'%(count, c, r, xp, y, zp)
                detfile.write(line + '\n')
                count = count + 1
        detfile.close()
        print >>sys.stderr, "Wrote detectors.txt"
        print >>sys.stderr, "Done"

    def get_calib_fitsfile(self):
        #SourcePlate/20130711_G3040_TwoTarget_IndiaNigro_SourcePlate_picture.fits
        ddir = self.pl.get_val('DATA_DIRECTORY')
        sname = self.pl.get_val('STUDY_NAME')
        root = self.pl.get_val('ROOT')
        fname = "{0}_SourcePlate_picture.fits".format(root)
        return os.path.join(ddir, sname, 'SourcePlate', fname)

def run_g3calib(initfilename, outputfilename):
    g3c = G3Calib(initfilename, outputfilename) 

def mult(x, A):
    return np.dot(A, x)

def xz2cr(xz, pl):
    dp = float(pl.get_val("CALIBRATED_CCD_DPIXEL"))
    ct = float(pl.get_val("CALIBRATED_CCD_COS_THETA"))
    st = float(pl.get_val("CALIBRATED_CCD_SIN_THETA"))
    x0 = float(pl.get_val("CALIBRATED_CCD_X0"))
    z0 = float(pl.get_val("CALIBRATED_CCD_Z0"))
    v0 = np.array([x0, z0])

    R = np.array([ct, st, -st, ct]).reshape((2,2))

    #xz = np.loadtxt(srcfile)
    #xz = xz[:,3:6:2]
    ns = len(xz)

    xz = xz/dp              # scale
    cr = np.zeros(ns*2).reshape((ns, 2))
    for i in xrange(ns):
        cr[i] = mult(xz[i], R)
    cr = cr + v0
    return cr
