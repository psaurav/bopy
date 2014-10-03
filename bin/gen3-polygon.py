#!/usr/bin/python

import os
import sys
import re

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.nxutils import points_inside_poly as pip
from matplotlib.patches import Polygon

import bopy as bp

class PolygonMask:
    
    def __init__(self, filename):
        self.set_up(filename)
        fig = plt.figure()
        kid = fig.canvas.mpl_connect('key_press_event', self.onkey)
        plt.imshow(self.data, interpolation='nearest')
        plt.show()

    def onkey(self, event):
        print event.key
        if event.key == 'm':
            self.get_polygon('polygon-mask.npy')
    
    def get_polygon(self, maskfile):
        polygon = plt.ginput(0, 0)
        plt.axes().add_patch(Polygon(polygon, closed=True, fill=False, linestyle='dotted', linewidth=1.0))
        plt.draw()
        polygon = np.array(polygon)
        polygon = (polygon + 0.5).astype(int).astype(float)

        idx = np.indices(self.data.shape)
        x = idx[0].flatten()
        y = idx[1].flatten()
        xyindices = np.vstack((y,x)).T
        mask = pip(xyindices, polygon)
        mask = mask.reshape(self.data.shape)

        np.save(maskfile, mask)
        print >>sys.stderr, "Saved to ", maskfile


    def set_up(self, filename):
        
        data = bp.io.read_frame(filename)
        if len(data.shape) == 3:
            data = data[0,:,:]
        self.data = np.copy(data)
        
if __name__ == '__main__':
    (path, prog) = os.path.split(sys.argv[0])
    try:
        filename = sys.argv[1]
    except:
        print >>sys.stderr, "Usage:", prog, "<fits filename>"
        sys.exit()

    pm = PolygonMask(filename)
