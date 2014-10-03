#!/usr/bin/python

import sys
import bopy as bp

if __name__ == "__main__":

    try:
        filename = sys.argv[1]
        print sys.argv[0], filename
        bp.gen3pp.run_g3study(filename)
    except:
        bp.gen3pp.run_g3study()
