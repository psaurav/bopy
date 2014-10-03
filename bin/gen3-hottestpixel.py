#!/usr/bin/python

import os
import sys
import bopy.gen3pp.utils as bppu

if __name__ == '__main__':
    
    try:
        ddir = sys.argv[1]
        root = sys.argv[2]
        exs = sys.argv[3]
        wls = sys.argv[4]
        nsrc = int(sys.argv[5])
    except:
        sys.exit('Usage: {0} <ddir> <root> <"meas"> <"wls"> <nsrc>'.format(os.path.basename(sys.argv[0])))

    try:
        zupto = int(sys.argv[6])
    except:
        zupto = None

    exs = [v for v in exs.strip().split()]
    wls = [int(v) for v in wls.strip().split()]
    print exs
    print wls
    print zupto

    bppu.light_levels(ddir, exs, root, wls, zupto=zupto, ns=nsrc)

    
    

    
        
        
