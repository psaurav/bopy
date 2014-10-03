#!/usr/bin/python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import bopy as bp

def points(filename):
    pts = plt.ginput(0,0)
    x = map(lambda x: x[0], pts) # map applies the function passed as 
    y = map(lambda x: x[1], pts) # first parameter to each element of pts
    plt.plot(x, y, '-o')
    plt.draw()
    pts = np.array(pts)
    np.savetxt(filename, pts)
    

def onkey(event):
    print event.key
    if event.key == 'c':
        points('contact-line.txt')
    if event.key == 'f':
        points('free-line.txt')

def onclick(event):
    print event.key

def main(filename, cfgfilename):
    d = bp.io.read_frame(filename)
    if len(d.shape) > 2:
        d = d[0,:,:]
    pl = bp.utils.ParamList(cfgfilename)
    d = bp.utils.data_reorient_udrot(d, pl)
    fig = plt.figure()
    kid = fig.canvas.mpl_connect('key_press_event', onkey)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.imshow(d, interpolation='nearest', cmap = cm.Greys_r)
    plt.show()

if __name__ == '__main__':

    try:
        filename = sys.argv[1]
        cfgfilename = sys.argv[2]
    except:
        print >>sys.stderr, 'Usage:', os.path.basename(sys.argv[0]), '<breastfile> <cfgfile>'
        sys.exit()
    
    main(filename, cfgfilename)
