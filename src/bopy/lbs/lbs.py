
import os
import re
import numpy as np
from scipy.interpolate import interp1d

def interpolate1d(x, y, xs, kind='linear'):
    f = interp1d(x, y, kind=kind)
    return f(xs)

def read_header(filename):
    f = open(filename, 'r')
    header = []
    for l in f:
        header.append(l)
        if 'end of header' in l:
            header = [x.strip() for x in header]
            return header
    return None

def read_lasers(header):
    l = filter(lambda x:re.search(r'Laser names', x), header)
    return np.array([int(x) for x in l[0].split(':')[1].split('*')[0].strip().split()])

def read_source_detector_sep(header):
    l = filter(lambda x:re.search(r'Source-Detector', x), header)
    return float(l[0].split(':')[1].strip())

def get_phantom_properties_mm(name, wavel):
    tt = np.load(name)
    mua = interpolate1d(tt[:,0], tt[:,1], wavel)
    musp = interpolate1d(tt[:,0], tt[:,2], wavel)
    return mua, musp

def get_dcswitch_expected(seminf, rho, mua_calib, musp_calib, wbyv, ac_real=None):
    for mua, musp in zip(mua_calib, musp_calib):
        G = seminf.reflectance(mua, musp, rho, wbyv)
        if not ac_real is None:
            amp = G.real
        else:
            amp = np.abs(G)
        phs = np.angle(G)
        try:
            amp_expected = np.vstack((amp_expected, amp))
            phs_expected = np.vstack((phs_expected, phs))
        except:
            amp_expected = amp
            phs_expected = phs
    return amp_expected.T, phs_expected.T

class DCSwitchCalibration:
    def __init__(self, seminf, rho, mua_calib, musp_calib, wbyv, amp_data, phs_data, ac_real=None):
        for mua, musp in zip(mua_calib, musp_calib):
            G = seminf.reflectance(mua, musp, rho, wbyv)
            if not ac_real is None:
                amp = G.real
            else:
                amp = np.abs(G)
            phs = np.angle(G)
            try:
                amp_expected = np.vstack((amp_expected, amp))
                phs_expected = np.vstack((phs_expected, phs))
            except:
                amp_expected = amp
                phs_expected = phs
        #amp_calib = amp_expected.T
        #phs_calib = phs_expected.T
        #self._amp_calib = amp_calib/amp_data
        #self._phs_calib = phs_data - phs_calib 
        self._amp_calib = amp_data/amp_expected.T
        self._phs_calib = phs_data - phs_expected.T 
        print "Calibration prepared:"
        print '   ', self._amp_calib.shape, self._phs_calib.shape

    def calibrate_all(self, amp, phs):
        amp_c = amp/self._amp_calib
        phs_c = phs - self._phs_calib
        return amp_c, phs_c
    
    def calibrate_column(self, index, amp, phs):
        print amp.shape, phs.shape
        amp_c = amp*self._amp_calib[:,index]
        phs_c = phs - self._phs_calib[:,index]
        return amp_c, phs_c

def read_tis(filename, minwl, maxwl, skiprows=10):
    d = np.loadtxt(filename, skiprows=skiprows)
    d = d[d[:,0]>=minwl]
    d = d[d[:,0]<=maxwl]
    return d

def get_tis_avg(filenames, minwl, maxwl):
    # no check is done to see if the wavelengths are same
    for f in filenames:
        if not os.path.isfile(f):
            print f, 'does not exist'
        d = np.loadtxt(f, skiprows=10)
        d = d[d[:,0]>=minwl]
        d = d[d[:,0]<=maxwl]                        # remove unwanted wavelengths
        try:
            dd = np.vstack((dd, d[np.newaxis,:]))     # stack them together
        except:
            dd = d[np.newaxis,:]
    w = dd[0,:,0]                                     # get the wavelength column
    dd = np.nanmean(dd[:,:,1::], axis=0)              # get the mean of amplitudes 
    return np.hstack((w[:,np.newaxis], dd))           # recreate and return 

def read_dcswitch(filename, minfreq, maxfreq):
    d = np.loadtxt(filename, skiprows=16)
    d = d[d[:,0]>=minfreq]
    d = d[d[:,0]<=maxfreq]
    #d[:,0] = d[:,0]*2*np.pi*1.e6
    return d

def get_dcswitch_amplitude_phase_radians(filename, minfreq, maxfreq):
    d = read_dcswitch(filename, minfreq, maxfreq)
    amp = d[:,2::2]
    phs = d[:,1::2]
    phs = phs*np.pi/180.    # to radians
    return amp, phs

def get_dcswitch_avg(filenames, minfreq, maxfreq):

    # no check is done to see if the frequencies are same
    for f in filenames:
        d = np.loadtxt(f, skiprows=16)
        d = d[d[:,0]>=minfreq]
        d = d[d[:,0]<=maxfreq]                        # remove unwanted frequencies
        try:
            dd = np.vstack((dd, d[np.newaxis,:]))     # stack them together
        except:
            dd = d[np.newaxis,:]
    #w = dd[0,:,0]*2*np.pi*1.e6                       # get the frequency column and convert to omega
    w = dd[0,:,0]                                     # get the frequency column and convert to omega
    dd = np.nanmean(dd[:,:,1::], axis=0)              # get the mean of amplitudes and phases
    return np.hstack((w[:,np.newaxis], dd))           # recreate and return 

import glob

class ReadDCSwitch:
    def __init__(self, (minfreq, maxfreq)):
        self._minfreq = minfreq
        self._maxfreq = maxfreq

    def read(self, ddir, filename):
        filename = os.path.join(ddir, filename)
        d = np.loadtxt(filename, skiprows=16)
        d = d[d[:,0]>=self._minfreq]
        d = d[d[:,0]<=self._maxfreq]
        w = d[:,0]*2.*np.pi*1e6
        phs = d[:,1::2]
        amp = d[:,2::2]
        return w, amp, phs*np.pi/180.

    def read_freqs(self, ddir, filename):
        filename = os.path.join(ddir, filename)
        d = np.loadtxt(filename, skiprows=16)
        d = d[d[:,0]>=self._minfreq]
        d = d[d[:,0]<=self._maxfreq]
        return d[:,0]

    def read_avg(self, ddir, filename):
        filenames = '{0}/{1}'.format(ddir, filename)
        filenames = sorted(glob.glob(filenames))
        print 'DCSwitch files read:'
        print '   ', '\n    '.join(filenames) 
        for f in filenames:
            d = np.loadtxt(f, skiprows=16)
            d = d[d[:,0]>=self._minfreq]
            d = d[d[:,0]<=self._maxfreq]                        # remove unwanted frequencies
            try:
                dd = np.vstack((dd, d[np.newaxis,:]))     # stack them together
            except:
                dd = d[np.newaxis,:]
        w = dd[0,:,0]*2.*np.pi*1e6                                     # get the frequency column and convert to omega
        dd = np.nanmean(dd[:,:,1::], axis=0)              # get the mean of amplitudes and phases
        phs = dd[:,0::2]*np.pi/180. 
        amp = dd[:,1::2]
        return w, amp, phs            

    def read_from_header(self, ddir, file_name):
        filename = os.path.join(ddir, file_name)
        header = read_header(filename)
        lasers = read_lasers(header)
        rho =  read_source_detector_sep(header)
        return rho, lasers

class ReflectanceCorrrection:
    def __init__(self, filenames, minwl, maxwl):
        self._data = get_tis_avg(filenames, minwl, maxwl)

    def get_correction_vector():
        return self._data[0], self._data[:,1]

    def correct(d):
        return d/self._data[:,1]

def display_check_data(freqs, lasers, amp_data, phs_data, amp_expected, phs_expected):
    nrows = amp_data.shape[1]
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    f, axes = plt.subplots(nrows, 2, sharex='col')
    f.set_size_inches(20.5,40.5)
    for i in range(nrows):
        c = i % 2
        axes[i][0].plot(freqs, amp_data[:,i], 'b.', label='{0}'.format(lasers[i]))
        axes[i][0].plot(freqs, amp_expected[:,i], 'r-', label='{0}'.format(lasers[i]))
        axes[i][1].plot(freqs, phs_data[:,i], 'b.', label='{0}'.format(lasers[i]))
        axes[i][1].plot(freqs, phs_expected[:,i], 'r-', label='{0}'.format(lasers[i]))
        axes[i][0].yaxis.set_major_locator(mticker.MaxNLocator(3))
        axes[i][1].yaxis.set_major_locator(mticker.MaxNLocator(3))
        axes[i][0].set_ylabel('{0}'.format(lasers[i]))
        axes[i][0].yaxis.get_major_formatter().set_powerlimits((-1, 1))
        #axes[i][0].legend()
        #axes[i][1].legend()
    axes[0][0].set_title('Amplitude')
    axes[0][1].set_title('Phase')
    f.suptitle(r'Check data and Theory')
    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.show()

def display_data(data_type, freqs, lasers, amp_data, phs_data):
    nrows = amp_data.shape[1]
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    f, axes = plt.subplots(nrows, 2, sharex='col')
    f.set_size_inches(20.5,40.5)
    for i in range(nrows):
        c = i % 2
        axes[i][0].plot(freqs, amp_data[:,i], 'b.', label='{0}'.format(lasers[i]))
        axes[i][1].plot(freqs, phs_data[:,i], 'b.', label='{0}'.format(lasers[i]))
        axes[i][0].yaxis.set_major_locator(mticker.MaxNLocator(3))
        axes[i][1].yaxis.set_major_locator(mticker.MaxNLocator(3))
        axes[i][0].set_ylabel('{0}'.format(lasers[i]))
        axes[i][0].yaxis.get_major_formatter().set_powerlimits((-1, 1))
        #axes[i][0].legend()
        #axes[i][1].legend()
    axes[0][0].set_title('Amplitude')
    axes[0][1].set_title('Phase')
    f.suptitle(r'Data ({0}'.format(data_type))
    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.show()

