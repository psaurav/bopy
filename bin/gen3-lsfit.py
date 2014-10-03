#!/usr/bin/python

import os
import sys
import bopy.gen3pp.longseries as ls

if __name__ == '__main__':

    try:
        ddir = sys.argv[1]
        cfgfile = sys.argv[2]
        w = int(sys.argv[3])
        kind = sys.argv[4]
    except:
        print >>sys.stderr, 'Usage: {0} <ddir> <cfgfile> <w> <kind>'.format(os.path.basename(sys.argv[0]))
        sys.exit(0)

    print w, kind
    ls.lsfit(ddir, cfgfile, w, kind)
