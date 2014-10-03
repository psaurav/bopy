import os
import sys
import numpy as np

def node_in_breast(nodes, bs, method='cubic'):
    """
    Returns an bool array 
    
    Args:
        nodes - mesh nodes np.array (nn, 3)
        bs - breast and chest shape 
        method - interpolation method [default: cubic]
    """
    from scipy.interpolate import griddata
    zn = griddata(bs[:,0:2], bs[:,2], nodes[:,0:2], method=method)
    return self.nodes[:,2] < zn
