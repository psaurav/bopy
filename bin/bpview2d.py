#!/usr/bin/python

import sys
import bopy as bp
import numpy as np

if __name__ == '__main__':
    try:
        filename = sys.argv[1]
    except:
        sys.exit("Usage: "+os.path.basename(sys.argv[0])+" <filename>")

    d = bp.io.read_frame(filename)
    bp.viewer.view2d_udlr(d)
