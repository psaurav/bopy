#!/usr/bin/python

import os
import sys
import bopy.gen3pp.longseries as ls

try:
    plfname = sys.argv[1]
except:
    print >>sys.stderr, 'Usage: {0} <param filename>'.format(os.path.basename(sys.argv[0]))

if not os.path.isfile(plfname):
    ls.print_defaults(plfname=plfname)
    sys.exit(0)

ls.driver(plfname)
ls.driver(plfname, correction='calib')
ls.driver(plfname, correction='pickoff')
