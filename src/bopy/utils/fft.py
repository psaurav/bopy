import numpy as np

def get_fft_fap(x, delta_t):
    """
        All FFT stuff done here.
        returns frequencies, amplitudes and phases
    """
    n = len(x)                      # length of data array
    nf = 1.0/delta_t/2.0            # Nyquist frequency
    f = np.linspace(0, nf, n/2+1)   # List of frequencies
    y = np.fft.rfft(x)              # Unnormalized Fourier coefficients (complex)
    a = np.abs(y)/(n/2.)             # Amplitude
    a[0] = a[0]/2.                   # Special normalization to compensate for power
    ph = np.angle(y)                # phase
    return f, a, ph


def get_fft_ap(s):
    s_fft = np.fft.rfft(s)
    a = np.abs(s_fft)/(s.size/2.)   # amplitude
    a[0] = a[0]/2.
    p = np.angle(s_fft)             # phase
    return a, p
    
    
