#!/usr/bin/env python

import sys
import os

from matplotlib.widgets import Cursor
import numpy as np
import matplotlib.pyplot as plt

#import bopy.io as bio
#from bopy.utils import ParamList, data_reorient_udrot

from bopy.gen3pp.pickoff import RowMax

try:
    ifile = sys.argv[1]
    ofile = sys.argv[2]
except:
    sys.exit('Usage: {0} <in-cfg-filename> <out-cfg-filename>'.format(os.path.basename(sys.argv[0])))

rowmax = RowMax(ifile, ofile)
data = rowmax.read_data() 

def onclick(event):
    if (event.button == 2):
        z0 = int(event.ydata + 0.5)
        print >>sys.stderr, 'z0:', z0
        rowmax.set_rowmax(z0) 
        return
    return

def onkey(event):
    if (event.key == 'w'):
        rowmax.write_cfg()
        print >>sys.stderr, 'Wrote to file:', ofile
        return
    return

fig = plt.figure()
ax = fig.add_subplot(111)
cid = fig.canvas.mpl_connect('button_press_event', onclick)
kid = fig.canvas.mpl_connect('key_press_event', onkey)

ax.imshow(data, interpolation='nearest')
cursor = Cursor(ax, useblit=True, color='red', linewidth=2 )
plt.show()
