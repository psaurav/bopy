#!/usr/bin/python

import sys
import numpy as np
import bopy as bp

if __name__=='__main__':
    #Build some strange looking data:
    try:
        filename = sys.argv[1]
    except:
        print >>sys.stderr, "Usage:", sys.argv[0], "<filename>"
        sys.exit(1)

    d = bp.io.read_frame(filename)
    fig_v=bp.viewer.Viewer2d_xsection(d)
