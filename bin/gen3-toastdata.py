#!/usr/bin/python

import os
import sys

if __name__ == '__main__':

    try:
        cfgfile1 = sys.argv[1]
        cfgfile2 = sys.argv[2]
    except:
        print >>sys.stderr, 'Usage: {0} <configure_calibration.cfg> <toast.cfg>'.format(os.path.basename(sys.argv[0]))
        sys.exit(1)

    import bopy.gen3pp.toast5_pickoff as toast
    toast.driver_toast_data(os.path.join('..', cfgfile1), os.path.join('..', cfgfile2))
