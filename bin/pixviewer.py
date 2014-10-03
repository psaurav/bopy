#!/usr/bin/python

import sys
import os
import bopy.io as bio
import matplotlib.pyplot as plt
import matplotlib.cm as cm

if __name__ == '__main__':
    try:
        filename = sys.argv[1]
        col = int(sys.argv[2])
        row = int(sys.argv[3])
    except:
        sys.exit("Usage: "+os.path.basename(sys.argv[0])+" <fitsfilename> <col> <row>")

    if os.path.exists(filename):
        d = bio.read_frame(filename)
    else:
        sys.exit("IOError: "+os.path.basename(sys.argv[0])+": file does not exist: "+filename)

    d = d[:,col,row] 

    plt.plot(d, '.')
    plt.show()
