#!/usr/bin/python

import re
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import bopy as bp

def view2d_udlr(d):
    """
    """
    dmax = np.max(d)
    dmin = np.min(d)
    N = 100
    deld = (dmax-dmin)/N

    global dxi
    global dmi
    dxi = N
    dmi = 0

    fig = plt.figure()
    ax = fig.add_subplot(111)

    im1 = ax.imshow(d)
    fig.colorbar(im1)

    def press(event):
        #print 'pressed:', event.key
        global dxi
        global dmi
        if event.key=='t':
            s = raw_input("Title: ")
            ax.set_title(s)
            fig.canvas.draw()
        if event.key=='x':
            s = raw_input("XLabel: ")
            ax.set_xlabel(s)
            fig.canvas.draw()
        if event.key=='y':
            s = raw_input("YLabel: ")
            ax.set_ylabel(s)
            fig.canvas.draw()
        if event.key=='w':
            s = raw_input("Save as filename: ")
            fig.savefig(s, bbox_inches='tight')
            print "Saved as x.png" 
        if event.key=='up':
            dxi = dxi + 1
            if dxi > N:
                dxi = N
            im1.set_clim([dmin+deld*dmi,dmin+deld*dxi])
            fig.canvas.draw()
        if event.key=='down':
            dxi = dxi - 1
            if dxi < 0:
                dxi = 0
            im1.set_clim([dmin+deld*dmi,dmin+deld*dxi])
            fig.canvas.draw()
        if event.key=='left':
            dmi = dmi - 1
            if dmi < 0:
                dmi = 0
            im1.set_clim([dmin+deld*dmi,dmin+deld*dxi])
            fig.canvas.draw()
        if event.key=='right':
            dmi = dmi + 1
            if dmi > N:
                dmi = N
            im1.set_clim([dmin+deld*dmi,dmin+deld*dxi])
            fig.canvas.draw()

    fig.canvas.mpl_connect('key_press_event', press)

    plt.show()

if __name__ == '__main__':
    try:
        filename = sys.argv[1]
    except:
        sys.exit("Usage: "+os.path.basename(sys.argv[0])+" <filename>")

    d = bp.io.read_frame(filename)
    view2d_udlr(d) 
