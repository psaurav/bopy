import os
import sys
import re
import numpy as np
import bopy as bp
import bopy.io as bio
from bopy.utils import struct

class Fit:

    def __init__(self, d, dt=0.1, f=1.0):
        self.d = np.array(d, type=float)
        self.t = np.arange(len(self.d))
        self.t = dt*self.t
        self.w = np.pi*2.0*f

    def fitAPD(self):
        D = np.mean(self.d)
        A = (np.max(self.d) - np.min(self.d))/2.0
        p = 0.0

        out = opt.leastsq(self.of_APD, [D, A, p], factor = 100, maxfev = 2000, xtol=1e-12, ftol=1e-12)
        x = out[0]
        return np.array([x[0], x[1], x[2]])
        

    def of_APD(self, x):
        return (self.d - (x[0]+x[1]*np.cos(self.w*self.t+x[2])) )
