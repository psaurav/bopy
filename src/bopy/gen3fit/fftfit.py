import sys
import numpy as np
import bopy as bp

class FFTfit:
    def __init__(self, freq, delt, n):
        if n >= 20:
            self.n = 20
        elif n < 20:
            if n >= 10:
                self.n = 10
            else:
                self.n = None
        self.nf = 1.0/delt/2.0            # Nyquist frequency
        self.freqs = np.linspace(0, self.nf, self.n/2+1)
        self.idx = np.argmin(np.abs(self.freqs-freq))
        print >>sys.stderr, '-----------------'
        print >>sys.stderr, 'FFTfit'
        print >>sys.stderr, '-----------------'
        print >>sys.stderr, 'Nyquist freq:', self.nf
        print >>sys.stderr, 'n           :', self.n
        print >>sys.stderr, 'idx         :', self.idx
        print >>sys.stderr, 'freqs       :', self.freqs 

    def fftfit_dap(self, d):
        d = bp.io.read_frame(d)
        s_fft = np.fft.rfft(d, n=self.n, axis=0)
        a = np.abs(s_fft)/(self.n/2.)      # amplitude normalization
        a[0] = a[0]/2.          # amplitude (f=0) additional normalization
        p = np.angle(s_fft)             # phase
        #p[p<0.0] += np.pi*2.0
        return a[0], a[self.idx], p[self.idx]
