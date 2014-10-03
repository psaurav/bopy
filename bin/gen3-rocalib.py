#!/usr/bin/python

import sys
import os
import bopy as bp

if __name__ == '__main__':
    
    try:
        filename = sys.argv[1]
        thickness = float(sys.argv[2])
    except:
        sys.exit("Usage: "+os.path.basename(sys.argv[0])+" <filename> <tankthickness(mm)> [<rebin>]")

    try:
        rebin = int(sys.argv[3])
    except:
        rebin = 1

    bp.viewer.rocalib(filename, thickness, rebin)
