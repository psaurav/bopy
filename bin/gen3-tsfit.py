#!/usr/bin/python

import os
import sys
import bopy.gen3pp.toaststudies as tstudy

if __name__ == '__main__':

    try:
        ddir = sys.argv[1]
        meas = sys.argv[2]
        cfgfile = sys.argv[3]
        w = int(sys.argv[4])
        kind = sys.argv[5]
    except:
        print >>sys.stderr, 'Usage: {0} <ddir> <meas> <cfgfile> <w> <kind>'.format(os.path.basename(sys.argv[0]))
        sys.exit(0)

    print meas, w, kind
    tstudy.tsfit(ddir, meas, cfgfile, w, kind)
