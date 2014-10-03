#!/usr/bin/python

import sys
import os
import bopy.io as bio
import bopy as bp

if __name__ == '__main__':
    try:
        filename = sys.argv[1]
    except:
        sys.exit("Usage: "+os.path.basename(sys.argv[0])+" <fitsfilename>")

    if os.path.exists(filename):
        d = bio.read_frame(filename)
    else:
        sys.exit("IOError: "+os.path.basename(sys.argv[0])+": file does not exist: "+filename)

    v = bp.viewer.View3D(d)
