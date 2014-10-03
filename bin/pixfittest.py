#!/usr/bin/python

import pyfits 
import pylab as plt
import sys
import numpy as np
import scipy

import re

import gen3py as g3
import bopy.io as bio


def read_gen3pp_txt(filename, x, y, nx, ny):
    myFile=file(filename,mode='rt')
    data=[]
    for line in myFile:
        data.append([float(value) for value in line.split()])

    d = np.reshape(data[0], (ny, nx))    
    
    return d[y,x]

# 
# main
#
try:
    filename = sys.argv[1]
    darkfile = sys.argv[2]
    afile = sys.argv[3]
    phifile = sys.argv[4]
    dfile = sys.argv[5]
    x = int(sys.argv[6])
    y = int(sys.argv[7])
    delt = float(sys.argv[8])
    f = float(sys.argv[9])
except:
    print "Usage:", sys.argv[0], "<fitsfile> <darkfile> <A-file> <phi-file> <D-file> <x1> <y1> <delt> <freq>"
    sys.exit()

#
# read fits data
#
fdata = bio.read_frame(filename)[:, y, x]
A = bio.read_frame(afile)[:, y, x]
p = bio.read_frame(phifile)[:, y, x]
D = bio.read_frame(dfile)[:, y, x]

#
# read fitting parameters
#
pl          = ParamList('configure.cfg')
delT        = float(pl.get_val('RF_DEL_TIME'))
f           = float(pl.get_val('RF_MODULATED_FREQ'))
doffset     = float(pl.get_val('DOFFSET'))
if doffset:
    fdata = fdata + doffset
    fdata2 = fdata2 + doffset
w           = 2.0*f*scipy.pi

# create timestamps
tmax        = delT*nz
t           = plt.arange(0.0, tmax, delT)
tf          = plt.arange(0.0, tmax, 0.001)
    
print t

#
# read fitted parameters
#
A           = g3.read_frame("A.gim")
A1 = A[y,x] 
A2 = A[y2,x2] 
phi         = g3.read_frame("phi.gim")
phi1= phi[y,x]
phi2= phi[y2,x2]
DC          = g3.read_frame("DC.gim")
DC1= DC[y,x]
DC2= DC[y2,x2]
dark        = g3.read_frame("darkmean3.gim")
dark1 = dark[y,x]
dark2 = dark[y2,x2]

fdata_corr  = fdata - dark1
fdata2_corr  = fdata2 - dark2

print "Fitted params"
print "  DC1    =", DC1
print "  DC2    =", DC2
print "  A1     =", A1
print "  A2     =", A2
print "  phi1   =", phi1
print "  phi2   =", phi2
print "  dark1  =", dark1
print "  dark2  =", dark2

func_fitted     = DC1 + A1*plt.cos(w*tf + phi1)
func_fitted2    = DC2 + A2*plt.cos(w*tf + phi2)

plt.plot(t, fdata_corr, 'ro')
plt.plot(tf, func_fitted)
plt.plot(t, fdata2_corr, 'ro')
plt.plot(tf, func_fitted2)
plt.legend(('Cordata1', 'param1', 'Cordata2', 'param2'),
           'upper right', shadow=True, fancybox=True)
plt.xlabel('time [sec]')
plt.ylabel('corrected CCD reading')
plt.title("Fitting pixel(%d, %d) (%d, %d)"%(x, y, x2, y2))
plt.show()


