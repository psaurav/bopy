#!/usr/bin/python

import os
import sys
import numpy as np
import pylab 
from bopy.utils import ParamList
from bopy.gen3pp import G3Calib

def onkey(event):
    if (event.key == 'u'):
        g3c.updown()
        pylab.imshow(g3c.data, interpolation='nearest')
        pylab.draw()
        return
    if (event.key == 'r'):
        g3c.rotate()
        pylab.imshow(g3c.data, interpolation='nearest')
        pylab.draw()
        return
    if (event.key == 'w'):
        g3c.write()
        return

def onclick(event):
    if (event.button == 3):
        x = float(int(event.xdata + 0.5))
        z = float(int(event.ydata + 0.5))
        g3c.set_xz(x, z)
        return
    if (event.button == 2):
        x0 = float(int(event.xdata + 0.5))
        z0 = float(int(event.ydata + 0.5))
        g3c.set_x0z0(x0, z0)
        return

if __name__ == '__main__':

    (path, prog) = os.path.split(sys.argv[0])
    try:
        initfilename = sys.argv[1]
        #filename = sys.argv[2]
        print >>sys.stderr, ""
        print >>sys.stderr, prog, "- get the orientation/calibration"
    except:
        #print >>sys.stderr, "Usage:", prog, "<init filename> <fits filename>"
        print >>sys.stderr, "Usage:", prog, "<init filename>"
        sys.exit()

    outputfilename = 'configure_calibration.cfg'
    g3c = G3Calib(initfilename, outputfilename)
    filename = g3c.get_calib_fitsfile()
    g3c.set_data(filename)

    fig = pylab.figure(1)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    kid = fig.canvas.mpl_connect('key_press_event', onkey)
    
    pylab.imshow(g3c.data, interpolation='nearest')
    pylab.axis('off')
    #pylab.colorbar()
    pylab.show()
