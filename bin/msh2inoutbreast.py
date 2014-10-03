#!/usr/bin/python

import os
import sys
import numpy as np

import bopy.io as bio
from bopy.breastshape.breast import BreastResult
from bopy.toast.utils import inbreast2nim

if __name__ == '__main__':
    
    try:
        mshfilename = sys.argv[1]
        brfilename = sys.argv[2]
        zchestwall = float(sys.argv[3])
    except:
        sys.exit('Usage: {0} <mesh filename> <breast result filename> <zchestwall>'.format(os.path.basename(sys.argv[0]))) 

    nodes = bio.read_mesh_nodes(mshfilename)
    br = BreastResult(brfilename)
    inbreast = []
    for p in nodes:
        ib = br.isInsideBreast(p[0], p[1], p[2], zchestwall=zchestwall)
        inbreast.append(ib)
    inbreast = np.array(inbreast)
    np.save('inbreast.npy', inbreast)
    inbreast2nim('inbreast.npy', 'in1out0.nim', 1, 0)
    inbreast2nim('inbreast.npy', 'refindex.nim', 1.44, 1.33)
