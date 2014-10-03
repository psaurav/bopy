#!/usr/bin/python

import os
import sys
from bopy.gen3pp.gen3pp_utils import togauss

if __name__ == '__main__':
    try:
        meas = sys.argv[1]
        sigma = float(sys.argv[2])
    except:
        sys.exit('Usage: {0} <meas> <sigma>'.format(os.path.basename(sys.argv[0])))
    togauss(meas, sigma)
        
