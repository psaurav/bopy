import sys
import numpy as np
from scipy.stats import chisquare 

import rpy2.robjects as ro
rqchisq = ro.r['qchisq']

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import bopy as bp

class View3D:
    def __init__(self, d):
        #self.d = d
        self.__dict__['d'] = d
        self.nmax = len(d[:,0,0])
        self.n = 0
        fig = plt.figure()
        kid = fig.canvas.mpl_connect('key_press_event', self.onkey)
        cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        plt.imshow(self.d[self.n,:,:], interpolation='nearest', cmap = cm.Greys_r)
        #plt.colorbar()
        plt.show()

    def onkey(self, event):
        if (event.key == 'n'):
            self.n = self.n + 1
            if self.n >= self.nmax:
                self.n = 0
            print >>sys.stderr, self.n, sys.getrecursionlimit()
            plt.imshow(self.d[self.n,:,:], interpolation='nearest', cmap = cm.Greys_r)
            #plt.colorbar()
            plt.draw()
        if (event.key == 'p'):
            self.n = self.n - 1
            if self.n < 0:
                self.n = self.nmax -1
            print >>sys.stderr, self.n, sys.getrecursionlimit()
            plt.imshow(self.d[self.n,:,:], interpolation='nearest', cmap = cm.Greys_r)
            #plt.colorbar()
            plt.draw()

    def onclick(self, event):
        if (event.button == 2):
            c = int(event.xdata + 0.5)
            r = int(event.ydata + 0.5)
            print >>sys.stderr, r, c
            plt.figure(2)
            plt.plot(self.d[:,r,c])
            plt.show()
        if (event.button == 3):
            
            # Get the pixel
            c = int(event.xdata + 0.5)
            r = int(event.ydata + 0.5)
            print >>sys.stderr, r, c

            # Extract time series for this pixel.
            s = self.d[:,r,c]
            n = len(s)
            delt_t = 0.1            # sampling rate fixed.  Might change in the future
            freq = 1.0              # x-correlation freq fixed. Might change
            omega = 2.*np.pi*freq
            wt = bp.gen3fit.get_wt(n, delt_t, freq)

            D1, A1, P1 = bp.gen3fit.pixel_fftfit(s, delt_t, freq=freq)
            #D2, C2, P2 = bp.gen3fit.pixel_cosfit(s, P1, delt_t, freq)
            D2, A2, P2 = bp.gen3fit.pixel_cosfit2(s, wt, P1)
            print >>sys.stderr, "FFT:", D1, A1, P1
            #print >>sys.stderr, "COS:", D2, C2*C2, P2
            print >>sys.stderr, "COS:", D2, A2, P2
            tt = np.arange(n)*delt_t
            t = np.arange(0., np.max(tt), 0.01)
            s_fft = D1 + A1*np.cos(omega*t + P1)
            #s_cos = D2 + C2*C2*np.cos(omega*t + P2)
            s_cos = D2 + A2*np.cos(omega*t + P2)
            s_exp = D2 + A2*np.cos(wt + P2)
            chsq, p = chisquare(s, s_exp, len(s)-3)

            pcts = bp.gen3fit.get_pcts(s, s_exp)
            #pval = bp.gen3fit.get_pvalue(s, s_exp, 3)

            # Now plot.
            fig = plt.figure(2)
            plt.show()
            fig.clf()
            plt.plot(tt, self.d[:,r,c], 'go', label='Data')
            plt.plot(t, s_fft, 'r-', label='FFT fit')
            plt.plot(t, s_cos, 'b-', label='COS fit')
            #plt.title('csq={0:.4f}, p={1:.4f}, pcts={2:.4f}, p={3:.4f}'.format(chsq, p, pcts, pval))
            plt.title('csq={0:.4f}, p={1:.4f}, pcts={2:.4f}'.format(chsq, p, pcts))
            plt.xlabel('t [sec]')
            plt.ylabel('CCD')
            plt.legend()
            plt.draw()
