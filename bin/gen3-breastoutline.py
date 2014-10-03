#!/usr/bin/python

import os
import sys
import re

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
import pyfits

#import bopy as bp

def read_fits(filename):
    """
    Read fits file

    Arguments:
        filename (str) -- name of FITS file

    Returns:
        d (numpy array) -- data
    """
    if not os.path.exists(filename):
        sys.exit("IOError: "+__name__+": File does not exist: "+filename)

    hdulist = pyfits.open(filename)
    n = len(hdulist) - 1
    d = hdulist[n].data.astype(float)
    if 'DOFFSET' in hdulist[n].header:
        doffset = float(hdulist[n].header['DOFFSET'])
        d = d + doffset
    return d

class BreastOutline:
    
    def __init__(self, filename):
        self.set_up(filename)
        #self.points = np.array([])
        fig = plt.figure()
        kid = fig.canvas.mpl_connect('key_press_event', self.onkey)
        cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        plt.imshow(self.data, interpolation='nearest')
        plt.show()

    def redraw(self, data, points, vmin, vmax):
        plt.clf()
        plt.imshow(data, interpolation='nearest', vmin=vmin, vmax=vmax)
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        try:
            plt.plot(points[:,1], points[:,0], color='c')
        except:
            pass
        plt.xlim((xmin, xmax))
        plt.ylim((ymin, ymax))
        plt.draw()

    def savepoints(self, points):
        f = open('breast-outline.txt', 'w')
        print >>f, '#', 'row', 'col'
        for p in points:
            print >>f, p[0], p[1]

    def onkey(self, event):
        print event.key
        if event.key == 'w':
            self.savepoints(self.points)
            return
        if event.key == 'down':
            self.vmax = self.vmax - (self.vmax-self.vmin)*0.1
            self.redraw(self.data, self.points, self.vmin, self.vmax)
            return
        if event.key == 'j':
            self.vmax = self.vmax - (self.vmax-self.vmin)*0.1
            self.redraw(self.data, self.points, self.vmin, self.vmax)
            return
        if event.key == 'k':
            self.vmax = self.vmax + (self.vmax-self.vmin)*0.1
            self.redraw(self.data, self.points, self.vmin, self.vmax)
            return
        if (event.key == 'u'):
            #g3c.updown()
            #self.redraw(self.data, self.points, self.vmin, self.vmax)
            return
        if (event.key == 'r'):
            #g3c.rotate()
            #self.redraw(self.data, self.points, self.vmin, self.vmax)
            return
    
    def onclick(self, event):
        print event.button
        if (event.button == 1):
            x = int(event.xdata + 0.5)
            z = int(event.ydata + 0.5)
            try:
                self.points = np.append(self.points, np.array([[z, x]]), axis=0)
            except:
                self.points = np.array([[z, x]])
            self.redraw(self.data, self.points, self.vmin, self.vmax)
            return
        if (event.button == 3):
            self.points = np.delete(self.points, -1, 0)
            self.redraw(self.data, self.points, self.vmin, self.vmax)
        if event.button == 2:
            self.vmax = self.vmax - (self.vmax-self.vmin)*0.1
            self.redraw(self.data, self.points, self.vmin, self.vmax)
            return

    def set_up(self, filename):
        #data = bp.io.read_frame(filename)
        data = read_fits(filename)
        if len(data.shape) == 3:
            data = data[0,:,:]
        self.data = np.copy(data)
        self.vmax = np.max(self.data)
        self.vmin = np.min(self.data)
        self.points = None

        
if __name__ == '__main__':
    (path, prog) = os.path.split(sys.argv[0])
    try:
        filename = sys.argv[1]
    except:
        print >>sys.stderr, "Usage:", prog, "<fits filename>"
        sys.exit()

    bol = BreastOutline(filename)
