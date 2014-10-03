#import sys
#import os
import numpy as np
#import math
import rpy2.robjects as ro
from rpy2.rinterface import FloatSexpVector as FVector
#from rpy2.robjects.vectors import FloatVector as FVector
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
import gc

def loess2d(d, span=0.6):
    s = d.shape
    x, y = np.indices(s)
    x = x.ravel()
    y = y.ravel()
    d = d.ravel()

    #r = ro.r
    stats = importr("stats")
    f = ro.Formula('d ~ x+y')
    env = f.environment
    env['x'] = FVector(x)
    env['y'] = FVector(y)
    env['d'] = FVector(d)
    l = stats.loess(f, span=span)
    p = stats.predict(l)
    gc.collect()
    return np.reshape(np.array(p), s)

def loess1d(d, span=0.6):
    s = d.shape
    if len(s) > 1:
        print >>sys.stderr, 'Loess1d: data shape', s, 'not acceptable'
        return None
    x = np.arange(s[0])
    stats = importr("stats")
    f = ro.Formula('d ~ x')
    env = f.environment
    env['x'] = FVector(x)
    env['d'] = FVector(d)
    l = stats.loess(f, span=span)
    p = stats.predict(l)
    gc.collect()
    return np.array(p)
