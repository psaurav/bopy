#!/usr/bin/python

import os
import sys
import re

import numpy as np
import matplotlib.pyplot as plt 

import bopy.io as bio
import bopy.utils as butils
import bopy.gen3pp.utils as g3pputils

class Points:
    
    def __init__(self, filename, cr=None):
        self.set_up(filename)
        #self.points = np.array([])
        fig = plt.figure()
        kid = fig.canvas.mpl_connect('key_press_event', self.onkey)
        cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        plt.imshow(self.data, interpolation='nearest', cmap='gray')
        self.cr = cr
        if self.cr is not None:
            plt.axes().autoscale(False)
            plt.plot(self.cr[:, 0], self.cr[:, 1], 'y.')
        plt.show()

    def redraw(self, data, points, vmin, vmax):
        plt.clf()
        plt.imshow(data, interpolation='nearest', vmin=vmin, vmax=vmax, cmap='gray')
        if self.cr is not None:
            plt.axes().autoscale(False)
            plt.plot(self.cr[:, 0], self.cr[:, 1], 'y.')
        try:
            plt.axes().autoscale(False)
            plt.plot(points[:,1], points[:,0], color='c')
        except:
            pass
        plt.draw()

    def savepoints(self, points):
        f = open('points.txt', 'w')
        print >>f, '#', 'row', 'col'
        for p in points:
            print >>f, p[0], p[1]

    def remove_srcs(self, points):
        rows, cols = self.cr.shape
        for i in range(rows):
            s = i + 1
            c = self.cr[i,0]
            r = self.cr[i,1] 
            if g3pputils.point_in_poly(r, c, points):
                try:
                    self.picked_srcs.append(s)
                except:
                    self.picked_srcs = [s, ]

        print self.picked_srcs
        np.savetxt('src-choice.txt', self.picked_srcs, fmt='%.0f')

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
            if self.vmax < self.vmin:
                self.vmax = self.vmin
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
        if (event.key == 'x'):
            self.remove_srcs(self.points)
            return
    
    def onclick(self, event):
        #print event.button
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
        data = bio.read_frame(filename)
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
        print >>sys.stderr, "Usage:", prog, "<fits filename> [<cfgfile> <sources.txt>]"
        sys.exit()

    try:
        cfgfile = sys.argv[2]
        srcfile = sys.argv[3]
        print >>sys.stderr, 'Read cfg, src filenames'
    except:
        pass
        
    try:
        pl = butils.ParamList(cfgfile)
        src = bio.read_sdpos(srcfile)
        src = np.delete(src, 1, 1)  # delete the y column
        cr = g3pputils.xz2cr(src, pl)
        print >>sys.stderr, ' ** cr.shape: ', cr.shape
        bol = Points(filename, cr=cr)
    except:
        bol = Points(filename)

